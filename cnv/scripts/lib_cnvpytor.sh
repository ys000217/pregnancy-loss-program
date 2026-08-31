#!/usr/bin/env bash
# CNVpytor needs a named genome whose contig IDs match the BAM (@SQ is NC_*).
# Bundled hg38 uses chr1 and will not call on this FASTA. Build GC+conf once.
cnvpytor_ref_dir="${WORKDIR}/ref"
cnvpytor_gc="${cnvpytor_ref_dir}/GRCh38.p14.gc.pytor"
cnvpytor_conf="${cnvpytor_ref_dir}/GRCh38.p14.cnvpytor_conf.py"
cnvpytor_genome_id="GRCh38p14_refseq"

cnvpytor_ensure_ref() {
  local fai="${REF_FASTA}.fai"
  if [[ ! -s "${REF_FASTA}" || ! -s "${fai}" ]]; then
    echo "ERROR REF_FASTA or .fai missing: ${REF_FASTA}" >&2
    exit 1
  fi
  mkdir -p "${cnvpytor_ref_dir}"
  # Directory lock (portable). Do not use RETURN trap + local var under set -u.
  local lock="${cnvpytor_ref_dir}/.gc.lockdir"
  while ! mkdir "${lock}" 2>/dev/null; do
    sleep 5
  done
  if [[ ! -s "${cnvpytor_gc}" ]]; then
    echo "Building CNVpytor GC file (once): ${cnvpytor_gc}"
    # Conda builds often omit bundled hg38 GC/mask. CNVpytor then logs ERROR and
    # returns 0 without running -gc. -download only satisfies that gate (chr1 files);
    # calling still uses -conf with NC_* names.
    if ! python3 -c "from cnvpytor.genome import Genome; import sys; sys.exit(0 if Genome.check_resources() else 1)"; then
      echo "Bundled CNVpytor GC/mask missing; writing placeholders (cnvpytor -download is broken on this conda build)"
      python3 "${_CNV_ROOT}/scripts/unstick_cnvpytor_resources.py"
    fi
    cnvpytor -root "${cnvpytor_gc}" -gc "${REF_FASTA}" -make_gc_file
    if [[ ! -s "${cnvpytor_gc}" ]]; then
      rmdir "${lock}" 2>/dev/null || true
      echo "ERROR GC pytor was not created: ${cnvpytor_gc}" >&2
      echo "Run: python3 scripts/unstick_cnvpytor_resources.py" >&2
      echo "Then: bash scripts/prepare_cnvpytor_ref.sh" >&2
      exit 1
    fi
  fi
  python3 "${_CNV_ROOT}/scripts/write_cnvpytor_conf.py" \
    --fai "${fai}" --gc "${cnvpytor_gc}" -o "${cnvpytor_conf}"
  rmdir "${lock}" 2>/dev/null || true
}

cnvpytor_primary_chroms() {
  awk -F '\t' '$1 ~ /^NC_0000/ {
    n = $1; sub(/^NC_0000/, "", n); sub(/\..*/, "", n);
    if (n+0 >= 1 && n+0 <= 24) print $1
  }' "${REF_FASTA}.fai"
}

cnvpytor_cmd() {
  cnvpytor -conf "${cnvpytor_conf}" "$@"
}

# BAM contig lengths match hg38 so CNVpytor auto-labels the sample 'hg38' and
# then opens bundled chr1 GC/mask (empty placeholders). Force our NC_* genome.
cnvpytor_rd_his_call() {
  local root="$1"
  local bam="$2"
  local outdir="$3"
  local sample="$4"
  shift 4
  local chroms=("$@")
  rm -f "${root}"
  echo "CNVpytor -rd ${sample} (${#chroms[@]} chroms)"
  cnvpytor_cmd -root "${root}" -rd "${bam}" -chrom "${chroms[@]}"
  # Must overwrite auto 'hg38' before -his; also turn mask off (placeholders are empty).
  echo "CNVpytor force genome ${cnvpytor_genome_id}"
  python3 "${_CNV_ROOT}/scripts/cnvpytor_force_genome.py" \
    "${cnvpytor_conf}" "${root}" "${cnvpytor_genome_id}"
  echo "CNVpytor -cgc from ${cnvpytor_gc}"
  cnvpytor_cmd -root "${root}" -cgc "${cnvpytor_gc}"
  local ls_out
  ls_out="$(cnvpytor_cmd -root "${root}" -ls 2>&1 || true)"
  echo "${ls_out}"
  if ! grep -q "${cnvpytor_genome_id}" <<<"${ls_out}"; then
    echo "ERROR pytor still not labeled ${cnvpytor_genome_id}" >&2
    exit 1
  fi
  local bin
  for bin in "${BIN_PRIMARY}" "${BIN_LARGE}"; do
    echo "CNVpytor his/partition/call bin=${bin}"
    cnvpytor_cmd -root "${root}" -his "${bin}" -chrom "${chroms[@]}"
    cnvpytor_cmd -root "${root}" -partition "${bin}" -chrom "${chroms[@]}"
    cnvpytor_cmd -root "${root}" -call "${bin}" -chrom "${chroms[@]}" > "${outdir}/${sample}.cnvpytor.${bin}.tsv"
  done
  if [[ ! -s "${outdir}/${sample}.cnvpytor.${BIN_PRIMARY}.tsv" ]]; then
    echo "ERROR CNVpytor wrote no calls: ${outdir}/${sample}.cnvpytor.${BIN_PRIMARY}.tsv" >&2
    echo "Run: cnvpytor -conf ${cnvpytor_conf} -root ${root} -ls" >&2
    exit 1
  fi
  echo "CNVpytor calls: $(wc -l < "${outdir}/${sample}.cnvpytor.${BIN_PRIMARY}.tsv") lines @ ${BIN_PRIMARY}"
}
