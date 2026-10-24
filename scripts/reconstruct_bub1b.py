"""Independent canonical CDS reconstruction for BUB1B SNVs."""
from pathlib import Path
import json,pysam
from Bio.Seq import Seq
from project_paths import project_root, required_env
R = project_root()
g=json.loads((R/'data/references/bub1b_ensembl.json').read_text())
t=next(t for t in g['Transcript'] if t['id']==g['canonical_transcript'].split('.')[0])
assert t['strand']==1
tr=t['Translation'];coords=[]
for e in sorted(t['Exon'],key=lambda x:x['start']):
 coords.extend(range(max(e['start'],tr['start']),min(e['end'],tr['end'])+1))
f=pysam.FastaFile(str(R/'data/references/ucsc/GRCh38.primary.no_chr.fa'))
cds=''.join(f.fetch('15',p-1,p) for p in coords);index={p:i for i,p in enumerate(coords)}
assert len(cds)==3153 and str(Seq(cds).translate()).endswith('*')
v=pysam.VariantFile(str(R/'work/primary.normalized.vcf.gz'));rows=[]
for r in v.fetch('15',g['start']-1,g['end']):
 if r.pos not in index or len(r.ref)!=1 or len(r.alts[0])!=1:continue
 i=index[r.pos];assert cds[i]==r.ref
 codon=cds[i//3*3:i//3*3+3];a=list(codon);a[i%3]=r.alts[0];alt=''.join(a)
 rows.append({'position':r.pos,'ref':r.ref,'alt':r.alts[0],'coding_position':i+1,'protein_position':i//3+1,'ref_codon':codon,'alt_codon':alt,'ref_aa':str(Seq(codon).translate()),'alt_aa':str(Seq(alt).translate())})
(R/'results/candidates/bub1b_cds_reconstruction.json').write_text(json.dumps({'canonical_transcript':g['canonical_transcript'],'cds_length':len(cds),'variants':rows},indent=2)+'\n')
print('CDS reconstruction complete; controlled outputs saved locally.')
