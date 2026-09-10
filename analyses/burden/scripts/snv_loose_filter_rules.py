#!/usr/bin/env python3
"""Loose SNV candidate filters for pregnancy-loss discovery.

Design goals
------------
1. Use ANNOVAR functional layer + in-silico scores fully (orthogonal axes).
2. ClinVar: keep tiers 1–4; drop only clear Benign (tier 5). Tier 5 is NOT
   synonymous — synonymous is a functional class, Benign is clinical class.
3. Two-layer rarity (loose): population AF gate + own-control absence.

Axes (do not collapse into one score)
-------------------------------------
A. Function   Func / ExonicFunc          → hard gate for discovery set
B. ClinVar    CLNSIG                     → exclude Benign-only; keep 1–4 + blank
C. Scores     SIFT / PolyPhen / CADD / REVEL → flags for ranking, not hard drop
D. Rarity     gnomAD AF + control AC=0   → two-layer candidate gate

ClinVar clinical tiers (this project)
-------------------------------------
  1  Pathogenic
  2  Likely_pathogenic
  3  Uncertain_significance (VUS)
  4  Likely_benign
  5  Benign

Loose ClinVar rule: keep if unannotated OR highest retained tier in {1,2,3,4};
drop only when the label set is Benign-dominated (tier 5 only, or Benign with
no P/LP/VUS token). Likely_benign is kept on purpose (loose).

Recommended primary analysis set
--------------------------------
  function_ok AND clinvar_ok AND pop_af_loose AND control_absent

Then aggregate at gene / panel level (burden), not single-site claims.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable

# ---------------------------------------------------------------------------
# ClinVar
# ---------------------------------------------------------------------------
CLINVAR_TIER = {
    "Pathogenic": 1,
    "Likely_pathogenic": 2,
    "Uncertain_significance": 3,
    "Likely_benign": 4,
    "Benign": 5,
}

# Tokens that are not ACMG 1–5 clinical significance (keep; do not auto-drop)
CLINVAR_NON_ACMG = {
    "Conflicting_interpretations_of_pathogenicity",
    "Conflicting_classifications_of_pathogenicity",
    "not_provided",
    "other",
    "drug_response",
    "risk_factor",
    "association",
    "protective",
    "Affects",
    "association_not_found",
}

KEEP_CLINVAR_TIERS = {1, 2, 3, 4}  # drop only clear Benign (5) under loose rule


def clinvar_tokens(raw: str | None) -> set[str]:
    if raw is None:
        return set()
    s = str(raw).strip()
    if s in {"", ".", "NA", "na", "nan", "-"}:
        return set()
    parts = re.split(r"[|,;/]+", s)
    return {p.strip() for p in parts if p.strip()}


def clinvar_tiers_present(raw: str | None) -> set[int]:
    tiers: set[int] = set()
    for t in clinvar_tokens(raw):
        # normalize common variants
        key = t.replace(" ", "_")
        if key in CLINVAR_TIER:
            tiers.add(CLINVAR_TIER[key])
        elif t in CLINVAR_TIER:
            tiers.add(CLINVAR_TIER[t])
    return tiers


def clinvar_ok_loose(raw: str | None) -> bool:
    """Keep blank / 1–4 / conflicting / non-ACMG; drop Benign-only (tier 5)."""
    toks = clinvar_tokens(raw)
    if not toks:
        return True
    tiers = clinvar_tiers_present(raw)
    if not tiers:
        # e.g. risk_factor only, conflicting without parsed ACMG word
        return True
    if tiers & {1, 2, 3}:
        return True
    if 4 in tiers and 5 not in tiers:
        return True
    if 4 in tiers and 5 in tiers:
        # mixed LB/B → still keep under loose (borderline)
        return True
    # only Benign (and maybe junk)
    if tiers == {5}:
        return False
    if 5 in tiers and not (tiers & KEEP_CLINVAR_TIERS):
        return False
    return True


# ---------------------------------------------------------------------------
# Function layer (hard gate)
# ---------------------------------------------------------------------------
FUNC_KEEP = {"exonic", "splicing", "exonic;splicing"}

LOF_EXONIC = (
    "stopgain",
    "stoploss",
    "frameshift",
    "startloss",
    "nonsense",
)


def consequence_class(func: str | None, exonic: str | None) -> str:
    """Return: splicing | lof | missense | other_nonsyn | synonymous | skip."""
    func = (func or "").strip()
    exonic = (exonic or "").strip().lower()
    if "splicing" in func:
        return "splicing"
    if any(x in exonic for x in LOF_EXONIC):
        return "lof"
    if "missense" in exonic or "nonsynonymous" in exonic:
        return "missense"
    if "synonymous" in exonic:
        return "synonymous"
    if func not in FUNC_KEEP:
        return "skip"
    if exonic and exonic not in {".", "unknown", ""}:
        return "other_nonsyn"
    if func == "splicing" or "splicing" in func:
        return "splicing"
    return "skip"


def function_ok_loose(func: str | None, exonic: str | None) -> bool:
    """Discovery hard gate: non-synonymous coding + splicing; drop synonymous."""
    cons = consequence_class(func, exonic)
    return cons in {"splicing", "lof", "missense", "other_nonsyn"}


# ---------------------------------------------------------------------------
# Prediction scores (soft flags; not hard exclude in loose mode)
# ---------------------------------------------------------------------------
@dataclass
class ScoreFlags:
    sift_damaging: bool = False
    polyphen_damaging: bool = False
    cadd_ge15: bool = False
    cadd_ge20: bool = False
    revel_ge05: bool = False
    revel_ge075: bool = False
    n_damaging_votes: int = 0
    labels: list[str] = field(default_factory=list)


def _parse_float(x) -> float | None:
    if x is None:
        return None
    s = str(x).strip()
    if s in {"", ".", "NA", "na", "nan", "-"}:
        return None
    # ANNOVAR often "0.12,D" or "D(0.12)"
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    if not m:
        return None
    try:
        return float(m.group(0))
    except ValueError:
        return None


def score_flags(
    sift: str | None = None,
    polyphen_hdiv: str | None = None,
    polyphen_hvar: str | None = None,
    cadd: str | None = None,
    revel: str | None = None,
) -> ScoreFlags:
    """Parse common ANNOVAR score columns into soft damaging flags."""
    f = ScoreFlags()
    sift_s = (sift or "").upper()
    if "D" in sift_s.split(",") or sift_s.startswith("D") or ";D" in sift_s:
        f.sift_damaging = True
        f.labels.append("SIFT_D")
    for raw in (polyphen_hdiv, polyphen_hvar):
        u = (raw or "").upper()
        if "D" in u or "P" in u.split(",")[0:1] or re.search(r"\bD\b|\bP\b", u):
            f.polyphen_damaging = True
            f.labels.append("PolyPhen_D/P")
            break
    cadd_v = _parse_float(cadd)
    if cadd_v is not None:
        if cadd_v >= 15:
            f.cadd_ge15 = True
            f.labels.append("CADD>=15")
        if cadd_v >= 20:
            f.cadd_ge20 = True
            f.labels.append("CADD>=20")
    revel_v = _parse_float(revel)
    if revel_v is not None:
        if revel_v >= 0.5:
            f.revel_ge05 = True
            f.labels.append("REVEL>=0.5")
        if revel_v >= 0.75:
            f.revel_ge075 = True
            f.labels.append("REVEL>=0.75")
    votes = [
        f.sift_damaging,
        f.polyphen_damaging,
        f.cadd_ge15,
        f.revel_ge05,
    ]
    f.n_damaging_votes = sum(1 for v in votes if v)
    return f


# ---------------------------------------------------------------------------
# Two-layer rarity (loose)
# ---------------------------------------------------------------------------
# Layer 1 — population database (gnomAD)
POP_AF_LOOSE = 0.01          # keep AF < 1% or missing
POP_AF_STRICT = 0.0          # keep only AF == 0 (missing handled separately)

# Layer 2 — own controls
# loose: zero alt carriers among phenotype==control
# (equivalent to control_AC == 0 for biallelic GT coding used here)


def pop_af_ok_loose(gnomad_af: float | None, threshold: float = POP_AF_LOOSE) -> bool:
    """Layer 1 loose: missing AF kept; observed AF must be < threshold."""
    if gnomad_af is None:
        return True
    return gnomad_af < threshold


def pop_af_ok_strict_absent(gnomad_af: float | None, allow_missing: bool = True) -> bool:
    """Layer 1 strict: AF == 0; missing optional (True = treat missing as absent)."""
    if gnomad_af is None:
        return allow_missing
    return gnomad_af <= POP_AF_STRICT


def control_absent_ok(control_carrier_n: int) -> bool:
    """Layer 2 (original): variant must not appear in any control carrier."""
    return int(control_carrier_n) <= 0


# Alternate Layer 2 — cohort frequency gate (parallel filter, does not replace control==0)
CTRL_CARRIER_RATE_MAX = 0.05
CASE_CARRIER_RATE_MIN = 0.05


def cohort_freq_gate_ok(
    n_control_carriers: int,
    n_control: int,
    n_case_carriers: int,
    n_case: int,
    *,
    ctrl_rate_max: float = CTRL_CARRIER_RATE_MAX,
    case_rate_min: float = CASE_CARRIER_RATE_MIN,
) -> bool:
    """Keep if control rate < 5% and combined-case (abnormal+normal) rate > 5%."""
    if n_control <= 0 or n_case <= 0:
        return False
    return (n_control_carriers / n_control) < ctrl_rate_max and (
        n_case_carriers / n_case
    ) > case_rate_min


# ---------------------------------------------------------------------------
# Combined decision
# ---------------------------------------------------------------------------
@dataclass
class FilterDecision:
    keep: bool
    function_ok: bool
    clinvar_ok: bool
    pop_af_ok: bool
    control_ok: bool
    consequence: str
    clinvar_tiers: set[int]
    score: ScoreFlags
    reasons: list[str]


def decide_loose(
    *,
    func: str | None,
    exonic: str | None,
    clnsig: str | None,
    gnomad_af: float | None,
    control_carrier_n: int,
    sift: str | None = None,
    polyphen_hdiv: str | None = None,
    polyphen_hvar: str | None = None,
    cadd: str | None = None,
    revel: str | None = None,
    require_control_absent: bool = True,
) -> FilterDecision:
    """Primary loose candidate rule for discovery.

    keep = function_ok & clinvar_ok & pop_af_loose & (control_absent if required)

    Scores are attached for ranking / secondary strata; they do not veto keep.
    """
    cons = consequence_class(func, exonic)
    f_ok = function_ok_loose(func, exonic)
    c_ok = clinvar_ok_loose(clnsig)
    p_ok = pop_af_ok_loose(gnomad_af)
    ctrl_ok = control_absent_ok(control_carrier_n) if require_control_absent else True
    scores = score_flags(sift, polyphen_hdiv, polyphen_hvar, cadd, revel)
    tiers = clinvar_tiers_present(clnsig)

    reasons: list[str] = []
    if not f_ok:
        reasons.append(f"drop_function:{cons}")
    if not c_ok:
        reasons.append("drop_clinvar_benign")
    if not p_ok:
        reasons.append(f"drop_pop_af:{gnomad_af}")
    if require_control_absent and not ctrl_ok:
        reasons.append(f"drop_in_control:n={control_carrier_n}")
    if f_ok and c_ok and p_ok and ctrl_ok:
        reasons.append("keep_loose_candidate")
        if scores.n_damaging_votes:
            reasons.append(f"score_votes:{scores.n_damaging_votes}")

    keep = f_ok and c_ok and p_ok and ctrl_ok
    return FilterDecision(
        keep=keep,
        function_ok=f_ok,
        clinvar_ok=c_ok,
        pop_af_ok=p_ok,
        control_ok=ctrl_ok,
        consequence=cons,
        clinvar_tiers=tiers,
        score=scores,
        reasons=reasons,
    )


def stratum_label(decision: FilterDecision) -> str:
    """Optional post-hoc stratum for tables/plots (not a hard filter)."""
    if not decision.keep:
        return "excluded"
    if decision.consequence in {"lof", "splicing"}:
        base = "lof_splice"
    elif decision.score.revel_ge075 or decision.score.n_damaging_votes >= 3:
        base = "missense_high"
    elif decision.score.n_damaging_votes >= 1:
        base = "missense_supported"
    else:
        base = "missense_unscored_or_tol"
    if decision.clinvar_tiers & {1, 2}:
        return base + "|clinvar_PLP"
    if 3 in decision.clinvar_tiers:
        return base + "|clinvar_VUS"
    if 4 in decision.clinvar_tiers:
        return base + "|clinvar_LB"
    return base + "|clinvar_blank"


# ---------------------------------------------------------------------------
# Human-readable rule card (print / README excerpt)
# ---------------------------------------------------------------------------
RULE_CARD = """
SNV loose dual-filter rule card
===============================

Clarify
  - ClinVar tier 5 (Benign) != synonymous SNV.
  - Synonymous = functional class (ExonicFunc); Benign = clinical class (CLNSIG).
  - Many Benign variants are missense; rare synonymous can still be Pathogenic.

A. Function (hard keep)
  Keep: splicing, stopgain/stoploss/frameshift/startloss, nonsynonymous/missense,
        other coding nonsyn (e.g. nonframeshift indel).
  Drop: synonymous; intergenic/intronic/UTR/ncRNA unless splicing.

B. ClinVar (loose)
  Tiers: 1 P, 2 LP, 3 VUS, 4 LB, 5 B.
  Keep: blank (unannotated), tiers 1–4, conflicting / non-ACMG labels.
  Drop: Benign-only (tier 5 without 1–3).

C. Prediction scores (soft; do not hard-drop in loose mode)
  Flag damaging if any of:
    SIFT D; PolyPhen D/P; CADD >= 15 (stronger >= 20); REVEL >= 0.5 (stronger >= 0.75)
  Use votes for ranking / strata: lof_splice > missense_high > missense_supported
                                 > missense_unscored_or_tol

D. Two-layer rarity (loose discovery set)
  Layer 1 population (gnomAD_genome_AF_ALL):
    keep if AF is missing OR AF < 0.01
  Layer 2 own control (original discovery set):
    keep if control carrier count == 0
  Primary set: A & B & Layer1 & Layer2

  Alternate Layer 2 (parallel set; does not replace the original):
    control carrier rate < 5%, and case (abnormal+normal combined) rate > 5%

E. Analysis recommendation
  Candidate variants → gene/panel burden (case vs control), not single-site
  pathogenicity claims. Optional sensitivity: Layer1 strict AF==0 (missing kept
  or dropped as a separate run).
"""


def print_rule_card() -> None:
    print(RULE_CARD.strip())


if __name__ == "__main__":
    print_rule_card()
    # tiny self-check
    d = decide_loose(
        func="exonic",
        exonic="nonsynonymous SNV",
        clnsig="Uncertain_significance",
        gnomad_af=0.002,
        control_carrier_n=0,
        revel="0.8",
        cadd="22",
    )
    assert d.keep
    assert stratum_label(d).startswith("missense_high")
    d2 = decide_loose(
        func="exonic",
        exonic="synonymous SNV",
        clnsig=".",
        gnomad_af=None,
        control_carrier_n=0,
    )
    assert not d2.keep
    d3 = decide_loose(
        func="exonic",
        exonic="nonsynonymous SNV",
        clnsig="Benign",
        gnomad_af=0.0,
        control_carrier_n=0,
    )
    assert not d3.keep
    assert cohort_freq_gate_ok(0, 153, 25, 495)
    assert not cohort_freq_gate_ok(8, 153, 25, 495)
    assert not cohort_freq_gate_ok(0, 153, 24, 495)
    print("self-check OK")
