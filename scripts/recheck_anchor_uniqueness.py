#!/usr/bin/env python3
"""Recheck 31-base anchors against all contigs of a local 2bit reference.

Output coordinates and hashes may be case-specific; keep outputs private.
"""
from pathlib import Path
import argparse,datetime,hashlib,json,time

def main():
    import py2bit
    p=argparse.ArgumentParser();p.add_argument('--reference',type=Path,required=True);p.add_argument('--patterns',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise ValueError('Refusing to overwrite an earlier result')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    t0=time.monotonic()
    patterns={}
    for line in a.patterns.read_text().splitlines():
     fields=line.split()
     if len(fields)!=2:raise ValueError('Each anchor line must contain name and sequence')
     name,seq=fields
     if name in patterns or len(seq)!=31 or not set(seq)<=set('ACGT'):raise ValueError('Require unique labels and 31 uppercase A/C/G/T anchors')
     patterns[name]=seq
    if not patterns:raise ValueError('No anchors supplied')
    reference=py2bit.open(str(a.reference));hits={n:[] for n in patterns};total=0
    for contig,size in reference.chroms().items():
     seq=reference.sequence(contig).upper();total+=len(seq)
     for label,anchor in patterns.items():
      offset=0
      while True:
       found=seq.find(anchor,offset)
       if found<0:break
       hits[label].append({'contig':contig,'start_0based':found});offset=found+1
    reference.close()
    d={'performed_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-t0,'reference_bases_scanned':total,'reference_sha256':hashlib.sha256(a.reference.read_bytes()).hexdigest(),'patterns_sha256':hashlib.sha256(a.patterns.read_bytes()).hexdigest(),'exact_31mer_hits':hits,'scope':'All contigs in existing local hg38.2bit; exact anchors only. Not alignment, independent calling or phase.'};a.output.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'reference_bases_scanned':total,'elapsed_seconds':d['elapsed_seconds'],'pattern_hit_counts':{k:len(v) for k,v in hits.items()}}))

if __name__=="__main__":
    main()
