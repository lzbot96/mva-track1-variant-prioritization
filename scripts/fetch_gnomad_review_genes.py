"""Fetch public gene-level frequency records, never patient genotypes."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess,json,time,sys,os
from gnomad_public import fetch_gene
from project_paths import project_root, required_env
R = project_root()
HELPER=os.environ.get('GNOMAD_HELPER')
out=R/'data/references/gnomad_review_genes';out.mkdir(exist_ok=True)
genes=json.loads((R/'results/candidates/review_gene_manifest.json').read_text())['genes']
genes=[g for g in genes if g not in ('BUB1B','CCDC84','FAM178A')]
def fetch(gene):
 dest=out/(gene+'.json')
 if dest.exists():return {'gene':gene,'cached':True}
 query='query { gene(gene_symbol: '+json.dumps(gene)+', reference_genome: GRCh38) { gene_id variants(dataset: gnomad_r4) { variant_id exome { ac an af } genome { ac an af } } } }'
 request={'query':query,'max_items':1,'max_depth':1,'timeout_sec':90,'save_raw':True,'raw_output_path':str(dest)}
 time.sleep(6)
 if HELPER:
  result=subprocess.run([sys.executable,HELPER],input=json.dumps(request),capture_output=True,text=True)
  try:summary=json.loads(result.stdout)
  except ValueError:summary={'ok':False,'error':'helper returned no JSON'}
 else:
  try:summary=fetch_gene(gene,dest)
  except Exception as exc:summary={'ok':False,'error_type':type(exc).__name__}
 (out/(gene+'.receipt.json')).write_text(json.dumps(summary,indent=2)+'\n')
 print('Public gene frequencies:',gene,'ok=',summary.get('ok'),flush=True)
 return {'gene':gene,'ok':summary.get('ok',False)}
with ThreadPoolExecutor(1) as pool:status=list(pool.map(fetch,genes))
(out/'download_status.json').write_text(json.dumps(status,indent=2)+'\n')
