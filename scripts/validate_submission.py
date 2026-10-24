"""Validate a review CSV against source genotypes and the GRCh38 reference."""
from pathlib import Path
import csv,json,pysam
from project_paths import project_root, required_env
R = project_root()
rows=list(csv.DictReader((R/'results/track1_submission_draft.csv').open()))
assert 0<len(rows)<=10
f=pysam.FastaFile(str(R/'data/references/ucsc/GRCh38.primary.no_chr.fa'))
v=pysam.VariantFile(str(R/'work/primary.normalized.vcf.gz'));alleles=[]
for row in rows:
 assert row['proband_id']=='PROBAND01' and 0<float(row['epcr'])<=1
 assert row['finding_type']=='primary'
 for n in ('1','2'):
  ch=row['chrom_'+n]
  if not ch:continue
  assert ch.startswith('chr');chrom=ch[3:];p=int(row['pos_'+n]);ref=row['ref_'+n];alt=row['alt_'+n]
  assert f.fetch(chrom,p-1,p-1+len(ref))==ref
  found=[r for r in v.fetch(chrom,p-1,p) if r.pos==p and r.ref==ref and alt in r.alts]
  assert len(found)==1
  gt=next(iter(found[0].samples.values())).get('GT');assert gt==(0,1)
  alleles.append([ch,p,ref,alt])
assert len(set(tuple(x) for x in alleles))==len(alleles)
assert [float(r['epcr']) for r in rows]==sorted([float(r['epcr']) for r in rows],reverse=True)
receipt={'valid':True,'rows':len(rows),'alleles':alleles,'checks':['schema limits','EPCR range and order','reference allele agreement','source VCF allele and heterozygous genotype','no duplicates'],'submitted':False}
(R/'results/submission_validation.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k!='alleles'}))
