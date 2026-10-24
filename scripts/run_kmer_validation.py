from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess,json,hashlib,tempfile,os
from project_paths import project_root
R = project_root()
out=R/'results/read_validation';out.mkdir(exist_ok=True)
marker=out/'complete.json'
marker.unlink(missing_ok=True)
files=sorted((R/'data/raw').glob('*.fastq.gz'))
if not files: raise SystemExit('No FASTQ inputs found; refusing empty success receipt.')
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def inventory():
 return [{'name':p.name,'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in files]
inputs=inventory();anchor_hash=digest(R/'work/allele_kmers.tsv')
with tempfile.TemporaryDirectory(prefix='.scan-',dir=out) as temporary:
 stage=Path(temporary)
 def run(path):
  stem=path.name.removesuffix('.fastq.gz')
  with (out/(stem+'.log')).open('w') as log:
   p=subprocess.run([str(R/'work/scan_allele_kmers'),str(R/'work/allele_kmers.tsv'),str(path),str(stage/(stem+'.tsv')),str(stage/(stem+'.matched.fastq'))],stderr=log)
  if p.returncode:raise RuntimeError('K-mer scan failed; inspect controlled logs')
  return {stem+suffix:digest(stage/(stem+suffix)) for suffix in ('.tsv','.matched.fastq')}
 with ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(run,files))
 if inventory()!=inputs or digest(R/'work/allele_kmers.tsv')!=anchor_hash:
  raise RuntimeError('Scan inputs changed during execution')
 artifacts={k:v for result in results for k,v in result.items()}
 for name in artifacts:os.replace(stage/name,out/name)
 receipt={'version':1,'files':len(files),'completed':True,'anchor_sha256':anchor_hash,'inputs':inputs,'output_sha256':artifacts}
 temp=stage/'complete.json';temp.write_text(json.dumps(receipt,indent=2)+'\n');os.replace(temp,marker)
print('Read validation complete; files:',len(files))
