# Suggestive-87 × gsMap localization summary

- Date: 2026-08-25
- Trait name: `RPL_sug87`
- Background sumstats: `RPL_combine_gsmap.sumstats.gz` (204,641 SNPs)
- Boost: 70 / 87 suggestive hits mapped into sumstats, |Z| set to 8 (GWAS sign preserved)
- Mapping: 86/87 liftOver hg38→hg19; 85/87 rsID; 70/87 in HapMap3-overlapping sumstats

## Spatial LDSC (E16.5 mouse embryo)

| metric | value |
|---|---|
| n_spots | 121,767 |
| mean z² | 1.087 |
| median z² | 0.673 |
| max \|z\| | 4.65 |
| min spot p | 1.2e-3 |

## Cauchy (cell type)

Best `p_cauchy` = **0.116** (Sympathetic nerve). All cell types > 0.11 — **no significant localization**.

Consistent with prior RPL / OLA1-boost runs: S-LDSC needs polygenic fuel; boosting ~70 SNPs does not enable cell-type mapping.
