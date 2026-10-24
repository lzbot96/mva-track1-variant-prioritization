"""Aggregate exact anchors; read IDs are not deduplicated molecular fragments."""
from pathlib import Path
import csv,json,collections,hashlib
from project_paths import project_root, required_env
R = project_root()
D=R/'results/read_validation'
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
receipt=json.loads((D/'complete.json').read_text())
if receipt.get('version')!=1 or receipt.get('completed') is not True:
 raise ValueError('Missing current scan evidence; rerun scanner')
if digest(R/'work/allele_kmers.tsv')!=receipt['anchor_sha256']:
 raise ValueError('Anchor table changed since scan')
files=sorted((R/'data/raw').glob('*.fastq.gz'))
inputs=[{'name':p.name,'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in files]
if not files or inputs!=receipt['inputs'] or len(files)!=receipt['files']:
 raise ValueError('FASTQ inventory changed since scan')
artifacts={p.name.removesuffix('.fastq.gz')+suffix for p in files for suffix in ('.tsv','.matched.fastq')}
if set(receipt['output_sha256'])!=artifacts:
 raise ValueError('Incomplete scan artifact inventory')
for name in artifacts:
 if digest(D/name)!=receipt['output_sha256'][name]:raise ValueError('Scan output changed since completion')
expected={Path(r['file']).name.removesuffix('.fastq.gz'):int(r['num_seqs']) for r in csv.DictReader((R/'results/qc/fastq_stats.tsv').open(),delimiter='\t')}
if set(expected)!={p.name.removesuffix('.fastq.gz') for p in files}:raise ValueError('QC inventory differs from scan')
patterns=dict(line.split() for line in (R/'work/allele_kmers.tsv').read_text().splitlines())
totals=collections.defaultdict(collections.Counter); lanes={}; identifiers=collections.defaultdict(set); positions=collections.defaultdict(list)
for stem,n in expected.items():
 rows=list(csv.DictReader((D/(stem+'.tsv')).open(),delimiter='\t'))
 if len(rows)!=len(patterns) or {r['pattern'] for r in rows}!=set(patterns) or not all(int(r['file_records'])==n for r in rows):
  raise ValueError('Count table patterns or record totals differ from expected scan')
 lanes[stem]={r['pattern']:int(r['reads']) for r in rows}
 for r in rows:
  for k in ('reads','q20_reads','q30_reads'):totals[r['pattern']][k]+=int(r[k])
 with (D/(stem+'.matched.fastq')).open() as f:
  while True:
   h=f.readline()
   if not h:break
   seq=f.readline().strip();f.readline();qual=f.readline().strip()
   for name,p in patterns.items():
    ix=seq.find(p)
    if ix>=0:
     allele=name.rsplit('_',1)[0]
     identifiers[allele].add(h.split()[0])
     positions[name].append(ix+16)
summary={'method':'Exact 31-base reference/alternate anchors in all raw FASTQs, both orientations. No alignment or molecular deduplication; not coverage or clinical validation.','files':len(expected),'read_records_scanned':sum(expected.values()),'counts_by_pattern':totals,'counts_by_file':lanes,'unique_read_pair_identifiers_by_allele':{k:len(v) for k,v in identifiers.items()},'variant_position_in_matched_reads':{k:{'min':min(v),'max':max(v)} for k,v in positions.items()}}
(D/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('Read validation aggregation complete; files:',len(expected))
