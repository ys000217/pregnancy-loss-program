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

## Main results (committed under `tables/` / `plots/`)

### Catalog

| Metric | Value |
|--------|-------|
| Unique AnnotSV P/LP (ACMG 4+5, full) | **124** |
| Population-rare P/LP (no gnomAD AF or AF&lt;1%) | **45** |
| Strict rare (pop-rare + cohort AF&lt;5%) | **38** |
| Genome-wide PASS SV loci tested | **55,388** |

### Sample-level burden

- **SNV mut/Mb**: case vs control significant, but direction is case slightly *lower*; abnormal ≈ normal (not an abnormal-elevated burden story).
- **SV total (alt-carrying PASS SVs only)**: abnormal ≈ normal; vs control a modest difference in common polymorphic load (not P/LP).
- **P/LP SV counts (alt-carrying)**: no significant group difference.
- **Strict rare P/LP**: catalog size depends on cohort AF from true carriers.

### Per-locus enrichment

- **Strict abnormal-specific** (`fdr(ab vs ctrl)<0.05` AND `ab_rate>ctrl` AND `ab_rate>norm` AND `fdr(ab vs norm)<0.05`): **0 loci**.
- FDR(ab vs ctrl)&lt;0.05: 63 sites; most are shared case-vs-control patterns, not abnormal-specific.
- **All P/LP case vs control**: **2** loci with FDR&lt;0.05 — large chr1 DEL (ACMG5, high carrier rates) and **TMEM63C** small DEL (ACMG4). Treat as weak common-carrier signals, not rare pathogenic drivers.
- **Rare / strict-rare P/LP**: **0** FDR&lt;0.05 hits.

Full per-locus PASS table `sv_locus_enrichment_all_pass.tsv` (~14 MB) stays local under `figure2/burden_analysis/tables/` (not committed).

## Run

```bash
Rscript analyses/burden/scripts/run_all.R
# or stepwise
python analyses/burden/scripts/compute_burden.py
python analyses/burden/scripts/compute_sv_locus_enrichment.py
Rscript analyses/burden/scripts/analyze_and_plot.R
```

A working copy may also remain under `figure2/burden_analysis/` (`figure*/` is gitignored).

## EpiFactors catalog (this round)

Single script: `scripts/epifactors_catalog.py` (open in the IDE). Interpreter:

`analyses/burden/.venv/Scripts/python.exe`

```bash
analyses/burden/.venv/Scripts/python.exe analyses/burden/scripts/epifactors_catalog.py
# full rescan of ANNOVAR + VCFs:
analyses/burden/.venv/Scripts/python.exe analyses/burden/scripts/epifactors_catalog.py --catalog
```

Outputs: `tables/epifactors/` and `plots/epifactors/` (Nature-style PDF/PNG + source_data).
Group rates are descriptive only (no p-values). Panel is the full `figure2/EpiGenes_main.xlsx` library.

## Narrative figures (burden + SV loci + EpiFactors)

**8–10 weeks** (g8/g9/g10; n=47/356/110):

```bash
analyses/burden/.venv/Scripts/python.exe analyses/burden/scripts/epifactors_catalog.py --plot-only --gw8-10
analyses/burden/.venv/Scripts/python.exe analyses/burden/scripts/plot_narrative.py
```

Writes `plots/narrative/` (Fig. 1–5) and `tables/gw8_10/` (SV locus scan + derived EpiFactors).

**Full 648** with 30 CpG-cluster suspects relabeled abnormal (n=77/418/153):

```bash
analyses/burden/.venv/Scripts/python.exe analyses/burden/scripts/run_narrative_all_suspect.py
```

Writes `plots/narrative_all_suspect/` and `tables/suspect_abn/`.
