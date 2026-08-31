# SNV / SV mutation burden and genome-wide SV locus enrichment

Unbiased SNV burden (nonsynonymous mut/Mb), SV count burden, AnnotSV P/LP SVs
(including rare subsets), and genome-wide PASS SV per-locus Fisher enrichment.

Large files (VCF, AnnotSV, ANNOVAR) are not in git. Defaults read from
`D:/ONT/figure2` and `D:/ONT` (override with environment variables).

## Local inputs (not in git)

| File | Role |
|------|------|
| `figure2/sample_phenotype_648.tsv` | 648 samples: Condition (abnormal / normal / control) |
| `figure2/WGS_ONT_Intersection_648samples.vcf.gz` | Unbiased SNV callset |
| `figure2/S956.snp.annovar.hg38_multianno.txt.gz` | SNV ANNOVAR annotation |
| `clinical_649.GRCh38.correct.vcf` | Unbiased SV callset |
| `clinical_649.GRCh38.annotsv.tsv` | AnnotSV (ACMG, population AF) |

Environment variables:

- `BURDEN_ROOT`: phenotype and SNV data root (default `D:/ONT/figure2`)
- `BURDEN_OUT` / `BURDEN_PLOT`: output tables and plots (default `tables/`, `plots/`)

## Main results (summary)

- **SNV mut/Mb**: by Condition and case vs control (`tables/group_comparison_stats.tsv`)
- **SV total / P/LP / rare P/LP**: per-sample counts and group tests
- **Genome-wide PASS SV loci**: 55,388 sites; three pairwise Fishers with separate FDR
- **Abnormal-specific (strict)**: `fdr(ab vs ctrl)<0.05` AND `ab_rate>ctrl` AND `ab_rate>norm` AND `fdr(ab vs norm)<0.05` ? **0 loci**
- 63 sites with FDR(ab vs ctrl)<0.05; most are `case_vs_control` (normal also elevated), not abnormal-specific

The full per-locus table `sv_locus_enrichment_all_pass.tsv` (~14 MB) stays local
in `figure2/burden_analysis/tables/` and is not committed.

## Run

```bash
# Requires local access to the large input files
Rscript analyses/burden/scripts/run_all.R
# or stepwise
python analyses/burden/scripts/compute_burden.py
python analyses/burden/scripts/compute_sv_locus_enrichment.py
Rscript analyses/burden/scripts/analyze_and_plot.R
```

A working copy may also remain under `figure2/burden_analysis/` (`figure*/` is gitignored).
