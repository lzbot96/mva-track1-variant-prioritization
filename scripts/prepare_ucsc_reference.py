from pathlib import Path
import py2bit,hashlib,json
from project_paths import project_root, required_env
root = project_root()
source=root/'data/references/ucsc/hg38.2bit'
assert source.stat().st_size==835393456
ref=py2bit.open(str(source))
chroms=ref.chroms()
primary=[f'chr{i}' for i in range(1,23)]+['chrX','chrY','chrM']
dest=root/'data/references/ucsc/GRCh38.primary.no_chr.fa'
with dest.open('w') as f:
 for chrom in primary:
  f.write('>'+chrom[3:]+'\n')
  # bounded per-chromosome memory, uppercase to avoid soft-mask differences
  seq=ref.sequence(chrom).upper()
  for i in range(0,len(seq),60):f.write(seq[i:i+60]+'\n')
h=hashlib.sha256()
with source.open('rb') as f:
 for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
(root/'data/references/ucsc/reference_provenance.json').write_text(json.dumps({'source_url':'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit','sha256':h.hexdigest(),'source_bytes':source.stat().st_size,'primary_contigs':{c[3:]:chroms[c] for c in primary},'transformation':'Primary chromosomes only, strip chr prefix, uppercase; chrM becomes M. Alternate/decoy contigs preserved in original VCF but not in this reference.'},indent=2)+'\n')
print('Primary GRCh38 FASTA prepared from UCSC 2bit.')
