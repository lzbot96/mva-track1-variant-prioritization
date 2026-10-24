#!/bin/bash
# Reconstructed from successful early record 5. Explicit network stage.
set -euo pipefail
cd "${MVA_PROJECT_ROOT:?Set controlled run directory}"
mkdir -p data/references/clinvar
curl --fail --location --retry 3 --continue-at - --silent --show-error -o data/references/clinvar/clinvar_20260928.vcf.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20260928.vcf.gz
curl --fail --location --retry 3 --silent --show-error -o data/references/clinvar/clinvar_20260928.vcf.gz.tbi https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20260928.vcf.gz.tbi
curl --fail --location --retry 3 --silent --show-error -o data/references/clinvar/clinvar_20260928.vcf.gz.md5 https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20260928.vcf.gz.md5
