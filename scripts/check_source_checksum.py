#!/usr/bin/env python3
"""Compare current supplied-VCF bytes with the saved authenticated manifest.

Output hashes may identify controlled source data; keep outputs private.
"""
import argparse,hashlib,json,datetime
from pathlib import Path

def main():
    p=argparse.ArgumentParser();p.add_argument('--source-vcf',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise ValueError('Refusing to overwrite an earlier result')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    def entry(p):
     records=json.loads(p.read_text())['files'];hits=[x for x in records if Path(x['path']).name==a.source_vcf.name]
     if len(hits)!=1:raise ValueError('Source filename must identify exactly one manifest record')
     return hits[0]
    m=entry(a.manifest);r=entry(a.receipt);h=hashlib.sha256()
    with a.source_vcf.open('rb') as f:
     for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    actual=h.hexdigest();result={'checked_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_kind':'supplied VCF compressed bytes','current_sha256':actual,'authenticated_manifest_sha256':m['sha256'],'download_receipt_sha256':r['sha256'],'matches_authenticated_manifest':actual==m['sha256'],'matches_download_receipt':actual==r['sha256'],'current_size_bytes':a.source_vcf.stat().st_size,'matches_manifest_size':a.source_vcf.stat().st_size==m['size'],'historical_upstream_sha256_verified':r.get('upstream_sha256_verified'),'raw_FASTQ_hashes_recomputed':False};a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if not k.endswith('sha256')}))
    if not (result['matches_authenticated_manifest'] and result['matches_download_receipt'] and result['matches_manifest_size']):raise ValueError('Current source bytes differ from the saved reference records')

if __name__=="__main__":
    main()
