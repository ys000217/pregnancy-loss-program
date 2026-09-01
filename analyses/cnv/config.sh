#!/usr/bin/env bash
# Shared paths and 10x-paired CNV defaults. Source this from every script.
# Do not write into NGS_Rawdata or ONT raw directories.

set -euo pipefail

# --- cluster paths (5250028_songyang on jhinno) ---
export PROJECT_ID="${PROJECT_ID:-5250028_songyang}"
export SHARE_ROOT="${SHARE_ROOT:-/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang}"
# Analysis outputs on the same share as BAM/FASTA — not $HOME
export WORKDIR="${WORKDIR:-${SHARE_ROOT}/cnv_work}"
export MANIFEST="${MANIFEST:-${WORKDIR}/manifest.from_ont.tsv}"

# ONT BAM trees (00c rglob *merged.bam, then dedup by sample id)
export ONT_BAM_ROOT="${ONT_BAM_ROOT:-${SHARE_ROOT}/ONT_processed_data/modbam_Dorodo}"
export ONT_RAW_MERGED_ROOT="${ONT_RAW_MERGED_ROOT:-/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_Lab_shared_data/ONT_Rawdata}"

# Illumina WGS: raw PE FASTQ (read-only)
export NGS_RAWDATA="${NGS_RAWDATA:-/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_Lab_shared_data/ONT_Rawdata/NGS_Rawdata}"

# Reference: MUST match ONT BAM @SQ (RefSeq accessions NC_*, NT_* on this project)
export REF_FASTA="${REF_FASTA:-${SHARE_ROOT}/GRCh38.p14.fasta}"
# Sniffles TR BED is optional; default GIAB TR beds use chr1/chr2 and do NOT match RefSeq @SQ.
# Leave empty or point to a RefSeq-named TR BED. Empty = Sniffles runs without --tandem-repeats.
export TR_BED="${TR_BED:-}"
export ANNOTSV_GENOME="${ANNOTSV_GENOME:-GRCh38}"

export THREADS="${THREADS:-16}"

# --- 10x operating point ---
export BIN_PRIMARY=100000
export BIN_LARGE=500000
export MIN_SVLEN=50
export MIN_DEPTH_CNV=100000
export MIN_QUAL_SV=20
export MOSDEPTH_MAPQ=20

# Warning thresholds (qc.tsv pass=0). Nominal ~10x often lands at unique ~6–8x after dups.
export MIN_MEDIAN_DEPTH=8
export MIN_BREADTH_5X=0.80
# Hard abort only if coverage is too thin to call 100 kb CNVs. 0 = never abort on QC.
export QC_ABORT="${QC_ABORT:-1}"
export QC_ABORT_MEDIAN="${QC_ABORT_MEDIAN:-3}"
export QC_ABORT_BREADTH="${QC_ABORT_BREADTH:-0.40}"

export MERGE_RO=0.50
# Merge filters (set before mask probe so a failed probe cannot skip these under set -u)
export MAX_CNV_EVENT="${MAX_CNV_EVENT:-10000000}"
export MERGE_MASK_FRAC="${MERGE_MASK_FRAC:-0.50}"
# CNVpytor call QC (disable with CNVPYTOR_QC=0)
export CNVPYTOR_QC="${CNVPYTOR_QC:-1}"
export CNVPYTOR_Q0_MAX="${CNVPYTOR_Q0_MAX:-0.5}"
export CNVPYTOR_PN_MAX="${CNVPYTOR_PN_MAX:-0.5}"
export CNVPYTOR_EVAL_MAX="${CNVPYTOR_EVAL_MAX:-1e-4}"
# Spectre mosdepth window (bp). mosdepth itself has no default --by; Spectre needs ~1 kb.
export SPECTRE_MOSDEPTH_BY="${SPECTRE_MOSDEPTH_BY:-1000}"
# Sex chromosomes dropped in merge by default (autosomes only)
export KEEP_SEX_CHROM="${KEEP_SEX_CHROM:-0}"
# Hard mask: Prefer explicit HARD_MASK_BED, then share path, then repo ref/ (_CNV_ROOT from lib_source_config).
# Do not use BASH_SOURCE here — config is often sourced via /dev/stdin (CR strip) and that breaks set -e.
if [[ -z "${HARD_MASK_BED:-}" || ! -s "${HARD_MASK_BED}" ]]; then
  HARD_MASK_BED=""
  for _mask_cand in \
    "${_CNV_ROOT:-}/ref/hard_mask.grch38.refseq.bed" \
    "${SHARE_ROOT}/cnv/ref/hard_mask.grch38.refseq.bed" \
    "${WORKDIR:-}/ref/hard_mask.grch38.refseq.bed"
  do
    if [[ -n "${_mask_cand}" && -s "${_mask_cand}" ]]; then
      HARD_MASK_BED="${_mask_cand}"
      break
    fi
  done
  export HARD_MASK_BED
fi
unset _mask_cand 2>/dev/null || true
