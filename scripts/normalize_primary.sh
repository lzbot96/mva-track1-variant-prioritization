#!/bin/bash
# Reconstructed from successful historical record 6; only input path parameterized.
set -euo pipefail
cd "${MVA_PROJECT_ROOT:?Set controlled run directory}"
: "${MVA_INPUT_VCF:?Set relative input VCF path}"
bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M -Ou "$MVA_INPUT_VCF" 2> logs/primary_subset.log | bcftools norm -f data/references/ucsc/GRCh38.primary.no_chr.fa -c e -m -any --old-rec-tag ORIG -Oz --threads 4 -o work/primary.normalized.vcf.gz 2> logs/normalization.log
bcftools index work/primary.normalized.vcf.gz
