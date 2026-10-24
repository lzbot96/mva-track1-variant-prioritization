"""Synthetic-only regressions for integrity receipts and scan evidence binding."""
import csv,gzip,hashlib,json,os,pathlib,subprocess,sys,tempfile
S=pathlib.Path(__file__).resolve().parents[1]/'scripts'
def invoke(name,r,success):
 p=subprocess.run([sys.executable,'-B',str(S/name)],env={**os.environ,'MVA_PROJECT_ROOT':str(r)},capture_output=True,text=True)
 assert (p.returncode==0)==success,(name,p.returncode,p.stderr)
 return p
with tempfile.TemporaryDirectory() as directory:
 r=pathlib.Path(directory)
 for d in ('data/raw','work','results/read_validation','results/qc'):(r/d).mkdir(parents=True)
 # Current content must be hashed, independently of saved receipts.
 digest=hashlib.sha256(b'original').hexdigest()
 (r/'data/raw/synthetic.bin').write_bytes(b'modified')
 (r/'data/authenticated_manifest.json').write_text(json.dumps({'repo':'synthetic/repository','revision':'synthetic','files':[{'path':'synthetic.bin','size':8,'sha256':digest}]}))
 receipt_data={'repo':'synthetic/repository','revision':'synthetic','files':[{'path':'synthetic.bin','sha256':digest}]}
 (r/'data/download_receipt.json').write_text(json.dumps(receipt_data))
 invoke('verify_downloads.py',r,False)
 (r/'data/raw/synthetic.bin').write_bytes(b'original');invoke('verify_downloads.py',r,True)
 receipt_data['revision']='different';(r/'data/download_receipt.json').write_text(json.dumps(receipt_data))
 invoke('verify_downloads.py',r,False)
 (r/'data/download_receipt.json').write_text('{"files":[]}');invoke('verify_downloads.py',r,False)
 assert json.loads((r/'data/integrity_report.json').read_text())['complete'] is False
 executable=r/'work/scan_allele_kmers'
 subprocess.run(['c++','-std=c++17',str(S/'scan_allele_kmers.cpp'),'-lz','-o',str(executable)],check=True)
 anchors=r/'work/allele_kmers.tsv';anchors.write_text('synthetic_REF_F\t'+'A'*31+'\n')
 reads=r/'data/raw/synthetic.fastq.gz'
 with gzip.open(reads,'wt') as f:f.write('@synthetic\n'+'A'*31+'\n+\n'+'I'*31+'\n')
 for count,matched in [(r/'missing/count.tsv',r/'matched.fastq'),(r/'counts.tsv',r/'missing/matched.fastq')]:
  p=subprocess.run([str(executable),str(anchors),str(reads),str(count),str(matched)],capture_output=True)
  assert p.returncode!=0,'Scanner accepted unwritable outputs'
 (r/'results/qc/fastq_stats.tsv').write_text('file\tnum_seqs\nsynthetic.fastq.gz\t1\n')
 invoke('run_kmer_validation.py',r,True);invoke('aggregate_read_validation.py',r,True)
 old_anchor=anchors.read_text();anchors.write_text(old_anchor.replace('REF','ALT'))
 invoke('aggregate_read_validation.py',r,False);anchors.write_text(old_anchor)
 counts=r/'results/read_validation/synthetic.tsv';original=counts.read_text()
 counts.write_text(original.replace('synthetic_REF_F','other_REF_F'));invoke('aggregate_read_validation.py',r,False)
 # Even a receipt with a matching digest must reject mismatched pattern identifiers.
 marker=r/'results/read_validation/complete.json';receipt=json.loads(marker.read_text())
 receipt['output_sha256']['synthetic.tsv']=hashlib.sha256(counts.read_bytes()).hexdigest();marker.write_text(json.dumps(receipt))
 invoke('aggregate_read_validation.py',r,False)
 counts.write_text(original)
 # A failed rerun cannot leave the prior completion receipt usable.
 with gzip.open(reads,'wt') as f:f.write('@truncated\n')
 invoke('run_kmer_validation.py',r,False);assert not marker.exists()
 invoke('aggregate_read_validation.py',r,False)
print('Synthetic integrity, scanner I/O, output/anchor binding, and stale-receipt regressions passed.')
