"""Public gene-level gnomAD acquisition; reconstructed stdlib transport, no patient input.
Gene-only requests; dated source/query/hash receipts are saved with each retrieval.
"""
import argparse,json,re,urllib.request,hashlib,datetime
from project_paths import project_root

def gene_query(gene, populations=False):
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_.-]*',gene) or re.match(r'chr(?:[0-9]+|[XYM])[-_:]',gene,re.I):
        raise ValueError('Expected a gene symbol, not a coordinate/variant query')
    counts = 'ac an af populations { id ac an }' if populations else 'ac an af'
    return 'query { gene(gene_symbol: '+json.dumps(gene)+', reference_genome: GRCh38) { gene_id variants(dataset: gnomad_r4) { variant_id exome { '+counts+' } genome { '+counts+' } } } }'

def fetch_gene(gene,dest,populations=False):
    body=json.dumps({'query':gene_query(gene,populations=populations)}).encode()
    request=urllib.request.Request('https://gnomad.broadinstitute.org/api',data=body,headers={'Content-Type':'application/json','User-Agent':'MVA-public-gene-reference/1.0'})
    with urllib.request.urlopen(request,timeout=90) as response: result=json.load(response)
    if result.get('errors') or not isinstance(result.get('data',{}).get('gene'),dict):
        raise RuntimeError('Public gene API returned errors or no gene; no cache saved')
    if not isinstance(result['data']['gene'].get('variants'),list):
        raise RuntimeError('Unexpected public gene response schema')
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2)+'\n')
    receipt={'ok':True,'source':'public gene-level gnomAD API','query_scope':'gene only','query':gene_query(gene,populations=populations),'retrieved_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'canonical_url':'https://gnomad.broadinstitute.org/gene/'+result['data']['gene']['gene_id']}
    dest.with_suffix('.receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--gene',required=True);p.add_argument('--output',required=True,help='relative public reference JSON path under run directory');p.add_argument('--populations',action='store_true',help='Include public group AC/AN; raw group maxima are not FAF95');a=p.parse_args()
    root=project_root();dest=(root/a.output).resolve()
    if root not in dest.parents:raise SystemExit('Output must be inside controlled run directory')
    try: fetch_gene(a.gene,dest,populations=a.populations)
    except Exception as exc:raise SystemExit('Public gene acquisition failed: '+type(exc).__name__)
    print('Public gene-level reference JSON saved.')
