#!/usr/bin/env python3
"""
PMD detection faithfully implementing the dnmtools algorithm.

Reference: dnmtools-1.5.0 (Smith et al.)
  - src/analysis/pmd.cpp        — PMD pipeline
  - src/common/TwoStateHMM_PMD.cpp — Baum-Welch BetaBinomial HMM
  - src/common/EmissionDistribution.cpp — BetaBinomial MLE

Algorithm:
  1. Bin CpG-level counts into genomic windows (default 1kb)
  2. Two-state Beta-Binomial HMM with Baum-Welch training
  3. CpG desert handling (break HMM at large inter-bin gaps)
  4. Posterior decoding → classify bins: PMD (LMD) vs HMD
  5. Extract continuous PMD regions, merge nearby domains
  6. Output BED format (compatible with dnmtools naming)

Key improvements over Gaussian HMM:
  - Beta-Binomial emissions correctly model count data
  - Accounts for varying coverage depth per bin
  - CpG desert segmentation prevents cross-desert merging
  - Convergence tolerance-based early stopping
"""

import gzip
import argparse
import sys
import numpy as np
from scipy.special import betaln, digamma, gammaln


# ─── Beta-Binomial Emission Distribution ─────────────────────────────────

class BetaBinomial:
    """
    Beta-Binomial emission distribution for methylation count data.

    Log-likelihood of observing k methylated reads out of N total reads:
      log P(k | N, α, β) = ln C(N,k) + ln B(k+α, N-k+β) - ln B(α, β)

    Parameters α, β control mean (= α/(α+β)) and dispersion.
    Lower α+β → higher variance (overdispersion relative to binomial).
    """

    def __init__(self, alpha=1.0, beta=1.0):
        self.alpha = alpha
        self.beta = beta
        self._lnbeta = betaln(alpha, beta)

    def log_likelihood(self, meth_first, meth_second):
        """Log P(k, N-k | α, β) for bin with k=meth_first, N-k=meth_second."""
        k = int(meth_first)
        n_minus_k = int(meth_second)
        n = k + n_minus_k
        if n == 0:
            return 0.0
        # ln C(n,k) via gammaln
        ln_choose = (gammaln(n + 1) - gammaln(k + 1) -
                     gammaln(n_minus_k + 1))
        return ln_choose + betaln(self.alpha + k, self.beta + n_minus_k) - self._lnbeta

    def __call__(self, meth_pair):
        """Convenience: meth_pair = (meth_count, unmeth_count)."""
        return self.log_likelihood(meth_pair[0], meth_pair[1])

    def fit(self, vals_a, vals_b, weights):
        """
        M-step: estimate α, β from weighted observations.

        Uses method-of-moments in logit space:
          E[log(p)] ≈ ψ(α) - ψ(α+β)  →  alpha_rhs = Σ wᵢ·log(p̂ᵢ)/Σ wᵢ
          E[log(1-p)] ≈ ψ(β) - ψ(α+β) → beta_rhs = Σ wᵢ·log(1-p̂ᵢ)/Σ wᵢ

        Then Newton iteration:
          α_new = invψ(ψ(α_old+β_old) + alpha_rhs)
          β_new = invψ(ψ(α_old+β_old) + beta_rhs)
        """
        vals_a = np.asarray(vals_a)
        vals_b = np.asarray(vals_b)
        weights = np.asarray(weights)

        w_sum = weights.sum()
        if w_sum <= 0:
            return

        alpha_rhs = np.dot(weights, vals_a) / w_sum
        beta_rhs = np.dot(weights, vals_b) / w_sum

        self.alpha = 0.01
        self.beta = 0.01
        tolerance = 1e-10

        for _ in range(100):
            prev_alpha, prev_beta = self.alpha, self.beta
            psi_sum = digamma(prev_alpha + prev_beta)
            self.alpha = self._invpsi(psi_sum + alpha_rhs)
            self.beta = self._invpsi(psi_sum + beta_rhs)

            if (self._movement(self.alpha, prev_alpha) <= tolerance or
                self._movement(self.beta, prev_beta) <= tolerance):
                break

        self._lnbeta = betaln(self.alpha, self.beta)

    @staticmethod
    def _invpsi(x):
        """Inverse digamma function (binary search)."""
        L = 1.0
        Y = np.exp(x)
        while L > 1e-10:
            Y += L * np.sign(x - digamma(Y))
            L /= 2.0
        return Y

    @staticmethod
    def _movement(curr, prev):
        return abs(curr - prev) / max(abs(curr), abs(prev))


# ─── Two-State HMM (BetaBinomial emissions) ──────────────────────────────

class TwoStateHMM:
    """
    Two-state Hidden Markov Model for methylation domain segmentation.

    State 0 = foreground (PMD/LMD): hypomethylated → α_fg small, β_fg large
    State 1 = background (HMD): hypermethylated → α_bg large, β_bg small

    Transitions (log-scale):
      s_f = P(start in FG),  s_b = P(start in BG)
      ff, fb = P(FG→FG), P(FG→BG)
      bf, bb = P(BG→FG), P(BG→BG)
      ft, bt = P(FG→end), P(BG→end)
    """

    def __init__(self, min_prob=1e-10, tolerance=1e-4, max_iter=10,
                 verbose=False):
        self.MIN_PROB = min_prob
        self.tolerance = tolerance
        self.max_iter = max_iter
        self.verbose = verbose

        # Transition initial values (matching dnmtools defaults)
        self.p_sf = 0.5;   self.p_sb = 0.5
        self.p_ff = 0.99;  self.p_fb = 0.01
        self.p_bf = 0.01;  self.p_bb = 0.99
        self.p_ft = 1e-10; self.p_bt = 1e-10

        # Emission initial values (matching dnmtools)
        self.fg = BetaBinomial(0.05, 0.95)  # hypomethylated
        self.bg = BetaBinomial(0.95, 0.05)  # hypermethylated

    @staticmethod
    def _log_sum_log(p, q):
        """log(exp(p) + exp(q)) — numerically stable."""
        if np.isneginf(p): return q
        if np.isneginf(q): return p
        larger = max(p, q)
        smaller = min(p, q)
        return larger + np.log1p(np.exp(smaller - larger))

    def _forward_backward(self, meth_vals, start, end):
        """
        Forward-Backward algorithm on one segment [start, end).

        meth_vals: list of (meth_count, unmeth_count) per bin

        Returns:
          forward: array of (fg_score, bg_score) in log space
          backward: array of (fg_score, bg_score) in log space
          total_score: log-likelihood of the segment
        """
        n = end - start
        forward = np.zeros((n, 2))
        backward = np.zeros((n, 2))

        lp_sf = np.log(self.p_sf)
        lp_sb = np.log(self.p_sb)
        lp_ff = np.log(self.p_ff)
        lp_fb = np.log(self.p_fb)
        lp_ft = np.log(self.p_ft)
        lp_bf = np.log(self.p_bf)
        lp_bb = np.log(self.p_bb)
        lp_bt = np.log(self.p_bt)

        # ── Forward ──
        forward[0, 0] = lp_sf
        forward[0, 1] = lp_sb
        val = meth_vals[start]
        if val[0] + val[1] >= 1:
            forward[0, 0] += self.fg(val)
            forward[0, 1] += self.bg(val)

        for i in range(1, n):
            k = i - 1
            forward[i, 0] = self._log_sum_log(
                forward[k, 0] + lp_ff,
                forward[k, 1] + lp_bf)
            forward[i, 1] = self._log_sum_log(
                forward[k, 0] + lp_fb,
                forward[k, 1] + lp_bb)
            val = meth_vals[start + i]
            if val[0] + val[1] >= 1:
                forward[i, 0] += self.fg(val)
                forward[i, 1] += self.bg(val)

        total = self._log_sum_log(
            forward[n-1, 0] + lp_ft,
            forward[n-1, 1] + lp_bt)

        # ── Backward ──
        backward[n-1, 0] = lp_ft
        backward[n-1, 1] = lp_bt

        for k in range(n-1, 0, -1):
            i = k - 1
            fg_a = backward[k, 0]
            bg_a = backward[k, 1]
            val = meth_vals[start + k]
            if val[0] + val[1] >= 1:
                fg_a += self.fg(val)
                bg_a += self.bg(val)
            backward[i, 0] = self._log_sum_log(fg_a + lp_ff, bg_a + lp_fb)
            backward[i, 1] = self._log_sum_log(fg_a + lp_bf, bg_a + lp_bb)

        return forward, backward, total

    def _single_iteration(self, meth_vals, reset_points):
        """
        One Baum-Welch iteration: E-step + M-step.
        Returns total log-likelihood.
        """
        n_bins = len(meth_vals)

        # Prepare log-ratio data for emission fitting (vectorized)
        meth_arr = np.array(meth_vals, dtype=np.float64)
        totals = meth_arr[:, 0] + meth_arr[:, 1]
        has_data = totals >= 1
        vals_a = np.zeros(n_bins)
        vals_b = np.zeros(n_bins)
        p_vals = np.clip(meth_arr[:, 0] / np.maximum(totals, 1), 0.01, 0.99)
        vals_a[has_data] = np.log(p_vals[has_data])
        vals_b[has_data] = np.log1p(-p_vals[has_data])

        lp_ff = np.log(self.p_ff); lp_fb = np.log(self.p_fb)
        lp_bf = np.log(self.p_bf); lp_bb = np.log(self.p_bb)

        total_score = 0.0

        # Accumulators for transitions
        ff_vals = np.full(n_bins, -np.inf)
        fb_vals = np.full(n_bins, -np.inf)
        bf_vals = np.full(n_bins, -np.inf)
        bb_vals = np.full(n_bins, -np.inf)

        # Weight accumulator per state for emission M-step
        fg_weights = []
        bg_weights = []
        fg_a_vals = []
        fg_b_vals = []
        bg_a_vals = []
        bg_b_vals = []

        for seg_i in range(len(reset_points) - 1):
            seg_start = int(reset_points[seg_i])
            seg_end = int(reset_points[seg_i + 1])

            # Skip empty segments (can happen with adjacent deserts)
            if seg_end <= seg_start:
                continue

            forward, backward, score = self._forward_backward(
                meth_vals, seg_start, seg_end)
            total_score += score

            seg_n = seg_end - seg_start

            # Accumulate weighted observations for this segment
            fw_seg = np.empty(seg_n)
            bw_seg = np.empty(seg_n)
            for i in range(seg_n):
                fg_s = forward[i, 0] + backward[i, 0]
                bg_s = forward[i, 1] + backward[i, 1]
                denom = self._log_sum_log(fg_s, bg_s)
                fw_seg[i] = np.exp(fg_s - denom)
                bw_seg[i] = np.exp(bg_s - denom)

            # Estimate transitions (vectorized where possible)
            fw_k = forward[:-1, 0]  # indices 0..n-2
            bw_k = forward[:-1, 1]
            bw_i_0 = backward[1:, 0].copy()
            bw_i_1 = backward[1:, 1].copy()

            # Add emissions to backward terms for transitions
            seg_meth = meth_vals[seg_start:seg_end]
            for i in range(seg_n - 1):
                val = seg_meth[i + 1]
                if has_data[seg_start + i + 1]:
                    bw_i_0[i] += self.fg(val)
                    bw_i_1[i] += self.bg(val)

                idx = seg_start + i
                ff_vals[idx] = fw_k[i] + lp_ff + bw_i_0[i] - score
                fb_vals[idx] = fw_k[i] + lp_fb + bw_i_1[i] - score
                bf_vals[idx] = bw_k[i] + lp_bf + bw_i_0[i] - score
                bb_vals[idx] = bw_k[i] + lp_bb + bw_i_1[i] - score

            # Collect weighted observations
            seg_mask = has_data[seg_start:seg_end]
            for i in range(seg_n):
                if seg_mask[i]:
                    fg_weights.append(fw_seg[i])
                    bg_weights.append(bw_seg[i])
                    idx = seg_start + i
                    fg_a_vals.append(vals_a[idx])
                    fg_b_vals.append(vals_b[idx])
                    bg_a_vals.append(vals_a[idx])
                    bg_b_vals.append(vals_b[idx])

        # ── M-step: update transitions ──
        n_blocks = len(reset_points) - 1

        # log-sum-exp of transition accumulators
        def log_sum_exp_vec(vec):
            valid = vec[~np.isinf(vec)]
            if len(valid) == 0:
                return -np.inf
            m = valid.max()
            return m + np.log(np.sum(np.exp(valid - m)))

        ff_sum = log_sum_exp_vec(ff_vals[:-1])  # exclude last (no transition from last)
        fb_sum = log_sum_exp_vec(fb_vals[:-1])
        bf_sum = log_sum_exp_vec(bf_vals[:-1])
        bb_sum = log_sum_exp_vec(bb_vals[:-1])

        # -inf initialization already excludes blocks with no transitions,
        # so no (n_blocks-1) correction is needed (unlike C++ which uses 0-init)
        p_ff_new = np.exp(ff_sum)
        p_fb_new = np.exp(fb_sum)
        p_bf_new = np.exp(bf_sum)
        p_bb_new = np.exp(bb_sum)

        denom = p_ff_new + p_fb_new
        if denom > 0:
            self.p_ff = max(p_ff_new / denom - self.p_ft / 2.0, self.MIN_PROB)
            self.p_fb = max(p_fb_new / denom - self.p_ft / 2.0, self.MIN_PROB)

        denom = p_bf_new + p_bb_new
        if denom > 0:
            self.p_bf = max(p_bf_new / denom - self.p_bt / 2.0, self.MIN_PROB)
            self.p_bb = max(p_bb_new / denom - self.p_bt / 2.0, self.MIN_PROB)

        self.p_sb = (self.p_bb + self.p_fb) / 2.0
        self.p_sf = (self.p_bf + self.p_ff) / 2.0

        # ── M-step: update emissions ──
        self.fg.fit(fg_a_vals, fg_b_vals, fg_weights)
        self.bg.fit(bg_a_vals, bg_b_vals, bg_weights)

        return total_score

    def train(self, meth_vals, reset_points):
        """Baum-Welch training across all segments."""
        prev_total = -np.inf

        for it in range(self.max_iter):
            total = self._single_iteration(meth_vals, reset_points)

            if self.verbose:
                fg_mean = self.fg.alpha / (self.fg.alpha + self.fg.beta)
                bg_mean = self.bg.alpha / (self.bg.alpha + self.bg.beta)
                sys.stderr.write(
                    f"  iter {it+1:3d}: ll={total:.2f}  "
                    f"FG a={self.fg.alpha:.3f} b={self.fg.beta:.3f} mu={fg_mean:.3f}  "
                    f"BG a={self.bg.alpha:.3f} b={self.bg.beta:.3f} mu={bg_mean:.3f}  "
                    f"FB_len={1/self.p_fb:.0f} BF_len={1/self.p_bf:.0f}\n"
                )

            if abs(total - prev_total) < self.tolerance:
                if self.verbose:
                    sys.stderr.write(f"  Converged at iter {it+1}\n")
                break
            prev_total = total

        return prev_total

    def posterior_decode(self, meth_vals, reset_points):
        """
        Posterior decoding: classify each bin as FG (PMD) or BG (HMD).
        Returns (classes, scores) where classes[i]=True for FG state.
        """
        n_bins = len(meth_vals)
        forward_all = np.zeros((n_bins, 2))
        backward_all = np.zeros((n_bins, 2))

        for seg_i in range(len(reset_points) - 1):
            seg_start = int(reset_points[seg_i])
            seg_end = int(reset_points[seg_i + 1])
            if seg_end <= seg_start:
                continue
            fw, bw, _ = self._forward_backward(meth_vals, seg_start, seg_end)
            forward_all[seg_start:seg_end] = fw
            backward_all[seg_start:seg_end] = bw

        classes = np.zeros(n_bins, dtype=bool)
        scores = np.zeros(n_bins)

        for i in range(n_bins):
            fg_state = forward_all[i, 0] + backward_all[i, 0]
            bg_state = forward_all[i, 1] + backward_all[i, 1]
            classes[i] = fg_state > bg_state
            scores[i] = np.exp(fg_state - self._log_sum_log(fg_state, bg_state))

        return classes, scores


# ─── Data loading and binning ────────────────────────────────────────────

def bin_methylation_counts(bedmethyl_path, chrom_filter, bin_size=1000,
                           min_cov=1):
    """
    Read bedMethyl, bin CpGs, return (meth_counts, bin_regions, reset_points).

    Each bin: meth_counts[i] = (methylated_reads, unmethylated_reads)

    Desert handling: segments are broken at gaps > desert_size bp
    between consecutive non-empty bins.
    """
    desert_size = max(5000, 5 * bin_size)  # dnmtools default

    # ── First pass: load per-chromosome CpG data ──
    opener = gzip.open if bedmethyl_path.endswith('.gz') else open

    chrom_data = {}
    total_lines = 0

    sys.stderr.write(f"Loading: {bedmethyl_path}\n")
    with opener(bedmethyl_path, 'rt') as f:
        for line in f:
            if line.startswith('track') or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 11:
                continue

            chrom = cols[0]
            if chrom_filter and chrom not in chrom_filter:
                continue

            total_lines += 1
            cov = int(cols[9])
            if cov < min_cov:
                continue

            start_pos = int(cols[1])   # 0-based
            end_pos = int(cols[2])     # 1-based
            # Use exact Nmod count from column 11 rather than rounding
            # from percentage (column 10) to avoid ±1 read errors per CpG.
            # Fall back to rounding if Nmod column is missing or empty.
            if len(cols) > 11 and cols[11].strip():
                n_mod = int(cols[11])
            else:
                meth_pct = float(cols[10])
                n_mod = int(round(meth_pct / 100.0 * cov))

            if chrom not in chrom_data:
                chrom_data[chrom] = {'start_pos': [], 'end_pos': [],
                                     'n_mod': [], 'cov': []}

            chrom_data[chrom]['start_pos'].append(start_pos)
            chrom_data[chrom]['end_pos'].append(end_pos)
            chrom_data[chrom]['n_mod'].append(n_mod)
            chrom_data[chrom]['cov'].append(cov)

            if total_lines % 2_000_000 == 0:
                sys.stderr.write(f"  {total_lines/1e6:.1f}M CpGs loaded...\n")

    sys.stderr.write(f"  Done: {total_lines:,} CpGs across "
                     f"{len(chrom_data)} chromosomes\n")

    # ── Second pass: bin each chromosome ──
    chrom_order = sorted(chrom_data.keys(),
                         key=lambda c: len(chrom_data[c]['end_pos']),
                         reverse=True)

    all_meth_counts = []
    all_bins = []       # (chrom, start, end)
    all_reset_points = [0]

    for chrom in chrom_order:
        n_cpg = len(chrom_data[chrom]['end_pos'])
        sys.stderr.write(f"\n{chrom}: {n_cpg:,} CpGs\n")

        starts = np.array(chrom_data[chrom]['start_pos'])
        ends = np.array(chrom_data[chrom]['end_pos'])
        n_mod_vals = np.array(chrom_data[chrom]['n_mod'], dtype=np.int64)
        cov_vals = np.array(chrom_data[chrom]['cov'])

        # Bin (vectorized with bincount)
        chrom_bins = []
        chrom_meth_counts = []

        bin_idx = ends // bin_size
        bin_min = bin_idx.min()
        bin_shifted = bin_idx - bin_min
        n_bins_total = bin_shifted.max() + 1

        bin_meth = np.bincount(bin_shifted, weights=n_mod_vals,
                               minlength=n_bins_total)
        bin_cov = np.bincount(bin_shifted, weights=cov_vals,
                              minlength=n_bins_total)

        # Track per-bin actual CpG position extents (for boundary refinement)
        # Use 64-bit sentinels that fit within 64-bit range
        SENTINEL_MAX = np.iinfo(np.int64).max
        bin_min_start = np.full(n_bins_total, SENTINEL_MAX, dtype=np.int64)
        np.minimum.at(bin_min_start, bin_shifted, starts.astype(np.int64))
        bin_max_end = np.full(n_bins_total, -1, dtype=np.int64)
        np.maximum.at(bin_max_end, bin_shifted, ends.astype(np.int64))

        nonempty = bin_cov > 0
        nonempty_idx = np.where(nonempty)[0]
        bin_numbers = nonempty_idx + bin_min

        for bi, b in enumerate(bin_numbers):
            n_meth_val = bin_meth[nonempty_idx[bi]]
            n_total_val = bin_cov[nonempty_idx[bi]]
            n_meth_int = int(n_meth_val)
            n_total_int = int(n_total_val)
            n_unmeth = n_total_int - n_meth_int
            bin_start = int(b) * bin_size
            bin_end = bin_start + bin_size
            # Refined boundaries: actual CpG positions within this bin
            refined_start = int(bin_min_start[nonempty_idx[bi]])
            refined_end = int(bin_max_end[nonempty_idx[bi]])
            chrom_bins.append((chrom, bin_start, bin_end,
                               refined_start, refined_end))
            chrom_meth_counts.append((n_meth_int, n_unmeth))

        sys.stderr.write(f"  {len(chrom_bins):,} non-empty bins\n")

        if len(chrom_bins) == 0:
            continue

        # Desert handling within chromosome
        # Insert reset points at gaps > desert_size
        chrom_reset_offset = len(all_meth_counts)
        all_reset_points.append(chrom_reset_offset)

        for i in range(1, len(chrom_bins)):
            gap = chrom_bins[i][1] - chrom_bins[i-1][2]
            if gap > desert_size:
                all_reset_points.append(chrom_reset_offset + i)

        all_meth_counts.extend(chrom_meth_counts)
        all_bins.extend(chrom_bins)

        sys.stderr.write(f"  {len(chrom_meth_counts):,} bins retained\n")

        # Free memory
        del chrom_data[chrom]

    all_reset_points.append(len(all_meth_counts))

    # Remove duplicate reset points (adjacent deserts)
    rp_unique = [all_reset_points[0]]
    for rp in all_reset_points[1:]:
        if rp != rp_unique[-1]:
            rp_unique.append(rp)
    all_reset_points = rp_unique

    sys.stderr.write(f"\nTotal: {len(all_meth_counts):,} bins, "
                     f"{len(all_reset_points)-2} deserts\n")

    return all_meth_counts, all_bins, all_reset_points


# ─── PMD extraction ──────────────────────────────────────────────────────

def extract_pmds(classes, bins, meth_counts, reset_points=None,
                 min_length=100_000, merge_dist=None):
    """
    Extract PMDs from HMM classification.

    Parameters
    ----------
    classes : np.ndarray[bool] — True = FG state (candidate PMD)
    bins : list of (chrom, bin_start, bin_end, refined_start, refined_end)
    meth_counts : list of (meth, unmeth)
    reset_points : list of int — segment boundary indices (deserts/chrom breaks)
    min_length : int — minimum PMD size
    merge_dist : int — merge nearby PMDs (default: 2*bin_size)

    Returns list of (chrom, start, end, mean_meth, n_bins)
    """
    if merge_dist is None:
        merge_dist = 2000

    # Build per-segment extraction with boundary enforcement
    raw_pmds = []
    cross_desert_splits = 0

    if reset_points is None or len(reset_points) < 2:
        reset_points = [0, len(classes)]

    for seg_i in range(len(reset_points) - 1):
        seg_start = int(reset_points[seg_i])
        seg_end = int(reset_points[seg_i + 1])
        if seg_end <= seg_start:
            continue

        in_pmd = False
        pmd_start = None
        pmd_chrom = None
        meth_vals = []

        for i in range(seg_start, seg_end):
            if classes[i]:
                if not in_pmd:
                    pmd_start = bins[i][3]   # refined_start (actual first CpG)
                    pmd_chrom = bins[i][0]
                    meth_vals = []
                meth_vals.append(
                    meth_counts[i][0] / max(1, meth_counts[i][0] + meth_counts[i][1])
                )
                in_pmd = True
            else:
                if in_pmd and meth_vals:
                    end = bins[i-1][4]      # refined_end (actual last CpG)
                    n_bins = len(meth_vals)
                    reg_mean = np.mean(meth_vals)
                    length = end - pmd_start
                    raw_pmds.append({
                        'chrom': pmd_chrom, 'start': pmd_start, 'end': end,
                        'mean_meth': reg_mean, 'length': length, 'n_bins': n_bins
                    })
                    in_pmd = False

        # Close any open PMD at segment boundary
        if in_pmd and meth_vals:
            end = bins[seg_end - 1][4]     # refined_end of last bin in segment
            n_bins = len(meth_vals)
            reg_mean = np.mean(meth_vals)
            length = end - pmd_start
            raw_pmds.append({
                'chrom': pmd_chrom, 'start': pmd_start, 'end': end,
                'mean_meth': reg_mean, 'length': length, 'n_bins': n_bins
            })

    # Filter by minimum length
    pmds = [r for r in raw_pmds if r['length'] >= min_length]

    # Merge nearby PMDs on same chromosome (within same segment only —
    # cross-desert merges are prevented by segment-boundary enforcement above)
    if len(pmds) > 1:
        merged = [pmds[0]]
        for r in pmds[1:]:
            prev = merged[-1]
            if (r['chrom'] == prev['chrom'] and
                r['start'] - prev['end'] <= merge_dist):
                # Weighted merge
                total_len = r['end'] - prev['start']
                w1 = prev['length'] / total_len
                w2 = r['length'] / total_len
                prev['end'] = r['end']
                prev['mean_meth'] = prev['mean_meth'] * w1 + r['mean_meth'] * w2
                prev['length'] = total_len
                prev['n_bins'] += r['n_bins']
            else:
                merged.append(r)
        pmds = merged

    return [(r['chrom'], r['start'], r['end'], r['mean_meth'], r['n_bins'])
            for r in pmds]


# ─── Main ─────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="PMD detection — dnmtools algorithm (Beta-Binomial HMM)")
    parser.add_argument("bedmethyl", help="Input bedMethyl file (gzipped)")
    parser.add_argument("-o", "--output", default="pmd_regions.bed",
                        help="Output BED file")
    parser.add_argument("-b", "--bin-size", type=int, default=1000)
    parser.add_argument("-c", "--min-coverage", type=int, default=1)
    parser.add_argument("-l", "--min-length", type=int, default=100000)
    parser.add_argument("-i", "--max-iter", type=int, default=10,
                        help="Max Baum-Welch iterations (dnmtools default: 10)")
    parser.add_argument("-t", "--tolerance", type=float, default=1e-4,
                        help="Convergence tolerance")
    parser.add_argument("--chroms", type=str, default=None,
                        help="Comma-separated chromosomes")
    parser.add_argument("--include-sex", action="store_true",
                        help="Include sex chromosomes (chrX, chrY) in processing")
    parser.add_argument("--ucsc-names", action="store_true",
                        help="Use UCSC-style chromosome names (chr1..chr22) "
                             "instead of NCBI RefSeq accessions")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--merge-dist", type=int, default=2000,
                        help="Max gap to merge adjacent PMDs (default: 2×bin)")
    args = parser.parse_args()

    # Chromosomes to process — autosomes only by default.
    # Support both NCBI RefSeq (NC_000001.11) and UCSC (chr1) naming conventions.
    NCBI_AUTOSOMES = {
        "NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12",
        "NC_000005.10", "NC_000006.12", "NC_000007.14", "NC_000008.11",
        "NC_000009.12", "NC_000010.11", "NC_000011.10", "NC_000012.12",
        "NC_000013.11", "NC_000014.9", "NC_000015.10", "NC_000016.10",
        "NC_000017.11", "NC_000018.10", "NC_000019.10", "NC_000020.11",
        "NC_000021.9", "NC_000022.11",
    }
    UCSC_AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
    NCBI_SEX = {"NC_000023.11", "NC_000024.10"}
    UCSC_SEX = {"chrX", "chrY"}

    if args.ucsc_names:
        chrom_filter = UCSC_AUTOSOMES.copy()
        if args.include_sex:
            chrom_filter |= UCSC_SEX
            sys.stderr.write("Note: including sex chromosomes (chrX, chrY) — "
                             "results may be unreliable\n")
        sys.stderr.write("Using UCSC-style chromosome names\n")
    else:
        chrom_filter = NCBI_AUTOSOMES.copy()
        if args.include_sex:
            chrom_filter |= NCBI_SEX
            sys.stderr.write("Note: including sex chromosomes (NC_000023.11, "
                             "NC_000024.10) — results may be unreliable\n")

    if args.chroms:
        chrom_filter = set(args.chroms.split(","))

    # ── Step 1: Load & bin ──
    meth_counts, bins, reset_points = bin_methylation_counts(
        args.bedmethyl, chrom_filter, bin_size=args.bin_size,
        min_cov=args.min_coverage)

    if len(meth_counts) < 100:
        sys.stderr.write("ERROR: insufficient data (<100 non-empty bins)\n")
        sys.exit(1)

    # ── Step 2: Train HMM ──
    sys.stderr.write(f"\n{'='*60}\n")
    sys.stderr.write("Training Beta-Binomial HMM (dnmtools algorithm)\n")
    sys.stderr.write(f"  bins={len(meth_counts)}, segments={len(reset_points)-1}\n")
    sys.stderr.write(f"  max_iter={args.max_iter}, tolerance={args.tolerance}\n")

    hmm = TwoStateHMM(max_iter=args.max_iter, tolerance=args.tolerance,
                      verbose=args.verbose)
    hmm.train(meth_counts, reset_points)

    sys.stderr.write(f"\nTrained HMM parameters:\n")
    sys.stderr.write(f"  FG (PMD): a={hmm.fg.alpha:.3f}, b={hmm.fg.beta:.3f}, "
                     f"mu={hmm.fg.alpha/(hmm.fg.alpha+hmm.fg.beta):.4f}\n")
    sys.stderr.write(f"  BG (HMD): a={hmm.bg.alpha:.3f}, b={hmm.bg.beta:.3f}, "
                     f"mu={hmm.bg.alpha/(hmm.bg.alpha+hmm.bg.beta):.4f}\n")
    sys.stderr.write(f"  Expected domain lengths: FG={1/hmm.p_fb:.0f} bins, "
                     f"BG={1/hmm.p_bf:.0f} bins\n")

    # ── Step 3: Posterior decoding ──
    sys.stderr.write(f"\nPosterior decoding...\n")
    classes, scores = hmm.posterior_decode(meth_counts, reset_points)

    n_fg = classes.sum()
    sys.stderr.write(f"  FG (PMD) bins: {n_fg:,} / {len(classes):,} "
                     f"({100*n_fg/len(classes):.1f}%)\n")

    # ── Step 4: Posterior score distribution (Q7) ──
    if args.verbose:
        sys.stderr.write(f"\nPosterior score distribution:\n")
        for pct in [10, 25, 50, 75, 90, 95, 99]:
            sys.stderr.write(f"  {pct}th: {np.percentile(scores, pct):.4f}\n")
        sys.stderr.write(f"  FG mean score: {scores[classes].mean():.4f}\n")
        sys.stderr.write(f"  BG mean score: {scores[~classes].mean():.4f}\n")

    # ── Step 5: Extract PMDs ──
    pmds = extract_pmds(classes, bins, meth_counts,
                        reset_points=reset_points,
                        min_length=args.min_length,
                        merge_dist=args.merge_dist)

    sys.stderr.write(f"  PMDs (≥{args.min_length/1000:.0f}kb): {len(pmds)}\n")

    # ── Step 6: Write output ──
    with open(args.output, 'w') as out:
        out.write("#track name=PMDs description=\"dnmtools Beta-Binomial HMM "
                  "Partially Methylated Domains\"\n")
        out.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t"
                  "thickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\t"
                  "blockStarts\n")
        for i, (chrom, start, end, meth, n_bins) in enumerate(pmds):
            length = end - start
            name = f"PMD{i}:{meth:.4f}:{n_bins}"
            out.write(f"{chrom}\t{start}\t{end}\t{name}\t{n_bins}\t.\t"
                      f"{start}\t{end}\t0,0,255\t1\t{length}\t"
                      f"{int(length/1000)}\n")

    total_bp = sum(p[2] - p[1] for p in pmds)
    sys.stderr.write(f"\n{'='*60}\n")
    sys.stderr.write(f"Total: {len(pmds)} PMDs, {total_bp/1e6:.1f} Mb\n")
    sys.stderr.write(f"Mean PMD size: {total_bp/max(1,len(pmds))/1e3:.1f} kb\n")
    sys.stderr.write(f"Output: {args.output}\n")


if __name__ == "__main__":
    main()
