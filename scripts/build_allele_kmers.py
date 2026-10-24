"""Reconstruct historical 31-base SNV anchor logic from reviewer-supplied TSV."""
import csv
from pathlib import Path
import pysam
from Bio.Seq import Seq
from project_paths import project_root, required_env
R = project_root()
f = pysam.FastaFile(str(R/'data/references/ucsc/GRCh38.primary.no_chr.fa'))
patterns=[]; names=set()
with (R/required_env('MVA_ANCHOR_CANDIDATES')).open() as stream:
    for row in csv.DictReader(stream,delimiter='\t'):
        name=row['name']; pos=int(row['pos']); ref=row['ref']; alt=row['alt']
        if not name or any(c.isspace() for c in name) or name in names:
            raise ValueError('Anchor names must be unique nonempty tokens')
        names.add(name)
        if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT' or ref==alt:
            raise ValueError('This reconstruction supports distinct single-base A/C/G/T alleles only')
        wt=f.fetch(row['chrom'],pos-16,pos+15)
        if len(wt)!=31 or wt[15]!=ref or any(c not in 'ACGT' for c in wt):
            raise ValueError('Reference window/allele mismatch or ambiguous reference')
        for allele,seq in [('REF',wt),('ALT',wt[:15]+alt+wt[16:])]:
            patterns += [(f'{name}_{allele}_F',seq),(f'{name}_{allele}_R',str(Seq(seq).reverse_complement()))]
if not patterns: raise ValueError('Empty candidate table')
(R/'work/allele_kmers.tsv').write_text(''.join(n+'\t'+s+'\n' for n,s in patterns))
print('Anchor table saved locally; no sequences emitted.')
