#!/bin/bash
# Reconstructed from successful historical record 16.
set -euo pipefail
cd "${MVA_PROJECT_ROOT:?Set controlled run directory}"
bcftools view -t 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y -Ou data/references/clinvar/clinvar_20260928.vcf.gz | bcftools norm -f data/references/ucsc/GRCh38.primary.no_chr.fa -c x -m -any -Oz --threads 4 -o data/references/clinvar/clinvar_20260928.primary.normalized.vcf.gz 2> logs/clinvar_normalization.log
bcftools index data/references/clinvar/clinvar_20260928.primary.normalized.vcf.gz
