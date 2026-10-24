"""Reconstructed upstream MD5 verification; reads public reference only."""
import hashlib,json
from project_paths import project_root
R=project_root();p=R/'data/references/clinvar/clinvar_20260928.vcf.gz'
expected=p.with_name(p.name+'.md5').read_text().split()[0]
h=hashlib.md5()
with p.open('rb') as f:
    for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
if h.hexdigest()!=expected: raise SystemExit('ClinVar checksum mismatch')
(p.parent/'provenance.json').write_text(json.dumps({'url':'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20260928.vcf.gz','date':'2026-09-28','bytes':p.stat().st_size,'upstream_md5':expected,'md5_verified':True},indent=2)+'\n')
print('ClinVar matches upstream MD5.')
