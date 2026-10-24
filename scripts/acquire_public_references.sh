#!/bin/bash
# Explicit network stage reconstructed from successful historical records 3,7,18.
# HPO/UCSC URLs are mutable: compare hashes against authorized archived receipts.
set -euo pipefail
cd "${MVA_PROJECT_ROOT:?Set controlled run directory}"
mkdir -p data/references/ucsc data/references/hpo data/references/ensembl
curl --fail --location --retry 3 --silent --show-error -o data/references/ucsc/hg38.refGene.gtf.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
curl --fail --location --retry 3 --silent --show-error -o data/references/ucsc/hg38.2bit https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
curl --fail --location --retry 3 --silent --show-error -o data/references/hpo/hp.obo https://purl.obolibrary.org/obo/hp.obo
curl --fail --location --retry 3 --silent --show-error -o data/references/hpo/genes_to_phenotype.txt https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt
curl --fail --location --retry 3 --silent --show-error -o data/references/ensembl/Homo_sapiens.GRCh38.116.gff3.gz https://ftp.ensembl.org/pub/release-116/gff3/homo_sapiens/Homo_sapiens.GRCh38.116.gff3.gz
