# Rebuild meqtl_hits.tsv from Delahaye 2018 S6 + miscarriage/RPL GWAS windows,
# then liftOver hg19 -> GRCh38.
# Citation: Delahaye F et al. PLoS Genet. 2018;14(11):e1007785. PMID:30452450
#
# Preferred path: use existing metadata/meqtl_hits.tsv (already GRCh38).
# This script regenerates from S6 if you need a full refresh; requires:
#   - readxl
#   - UCSC liftOver + hg19ToHg38.over.chain.gz (see docs/DATA_PATHS.md)

stop("Prefer the curated GRCh38 metadata/meqtl_hits.tsv. To fully rebuild, run the liftover pipeline documented in docs/REFERENCES.md (UCSC liftOver).")
