from pathlib import Path
import py2bit,json
from project_paths import project_root, required_env
R = project_root()
pats=[l.split() for l in (R/'work/allele_kmers.tsv').read_text().splitlines()]
ref=py2bit.open(str(R/'data/references/ucsc/hg38.2bit'))
hits={n:[] for n,s in pats}
for chrom,size in ref.chroms().items():
 seq=ref.sequence(chrom).upper()
 for name,kmer in pats:
  start=0
  while True:
   pos=seq.find(kmer,start)
   if pos<0:break
   hits[name].append({'contig':chrom,'start_0based':pos});start=pos+1
(R/'results/read_validation/reference_kmer_specificity.json').write_text(json.dumps(hits,indent=2)+'\n')
for n,h in hits.items():print(n,len(h))
