#!/bin/bash
set -euo pipefail
cd "${MVA_PROJECT_ROOT:?Set MVA_PROJECT_ROOT to a writable controlled run directory}"
"${PYTHON:-python3}" - <<'PY'
from pathlib import Path
import time,gzip
p=Path('data/references/ensembl/Homo_sapiens.GRCh38.116.gff3.gz')
deadline=time.monotonic()+3600
while not p.exists() or p.stat().st_size!=108378077:
 if time.monotonic()>deadline:raise SystemExit('Annotation resource not ready within one hour.')
 time.sleep(10)
with gzip.open(p,'rb') as f:
 while f.read(8*1024**2):pass
print('Ensembl GFF3 integrity checked.',flush=True)
PY
bcftools csq -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y -f data/references/ucsc/GRCh38.primary.no_chr.fa -g data/references/ensembl/Homo_sapiens.GRCh38.116.gff3.gz --local-csq --phase a --ncsq 64 --threads 4 -Oz -o work/primary.csq.vcf.gz work/primary.normalized.vcf.gz 2> logs/csq.log
bcftools index work/primary.csq.vcf.gz
bcftools annotate --pair-logic exact -a data/references/clinvar/clinvar_20260928.primary.normalized.vcf.gz -c INFO/ALLELEID,INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT,INFO/GENEINFO,INFO/MC,INFO/CLNHGVS -Oz --threads 4 -o work/primary.csq.clinvar.vcf.gz work/primary.csq.vcf.gz 2> logs/normalized_clinvar_annotate.log
bcftools index work/primary.csq.clinvar.vcf.gz
printf 'Primary chromosome consequence and ClinVar annotation complete.\n'
