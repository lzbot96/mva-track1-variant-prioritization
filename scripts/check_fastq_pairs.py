"""Exhaustive read-name pairing check using SeqKit's FASTQ parser."""
from concurrent.futures import ThreadPoolExecutor
from itertools import zip_longest
from pathlib import Path
import json
import subprocess
import time

from project_paths import project_root, required_env
ROOT = project_root()

def check(lane):
    paths=[ROOT/required_env('MVA_FASTQ_TEMPLATE').format(lane=lane,mate=mate) for mate in (1,2)]
    procs=[subprocess.Popen(['seqkit','seq','--name','--threads','2',str(p)],stdout=subprocess.PIPE,stderr=subprocess.PIPE) for p in paths]
    report={'lane':lane,'paired_records':0,'identifier_mismatches':0,'unequal_record_count':0,'mate_role_mismatches':0,'lane_label_mismatches':0}
    started=time.monotonic()
    for a,b in zip_longest(procs[0].stdout,procs[1].stdout):
        if a is None or b is None:
            report['unequal_record_count']+=1
            continue
        report['paired_records']+=1
        names=[a.strip().split(),b.strip().split()]
        ids=[v[0].removesuffix(b'/1').removesuffix(b'/2') for v in names]
        if ids[0]!=ids[1]:report['identifier_mismatches']+=1
        for mate,fields in enumerate(names,1):
            if len(fields)>1 and fields[1].split(b':')[0] in (b'1',b'2'):
                report['mate_role_mismatches']+=int(fields[1].split(b':')[0])!=mate
            parts=fields[0].split(b':')
            if len(parts)>=7 and parts[3].isdigit():report['lane_label_mismatches']+=int(parts[3])!=lane
    for proc in procs:
        err=proc.stderr.read()
        if proc.wait()!=0:raise RuntimeError(f'SeqKit failed in lane {lane}; no QC pass recorded.')
    report['seconds']=round(time.monotonic()-started,1)
    report['pass']=not any(report[k] for k in ('identifier_mismatches','unequal_record_count','mate_role_mismatches','lane_label_mismatches')) and report['paired_records']>0
    path=ROOT/f'results/qc/pairing_lane_{lane:03d}.json'
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(f'Lane {lane}: exhaustive pairing check complete; pass={report["pass"]}',flush=True)
    return report

with ThreadPoolExecutor(max_workers=2) as pool:
    reports=list(pool.map(check,range(1,5)))
(ROOT/'results/qc/pairing_summary.json').write_text(json.dumps(reports,indent=2)+'\n')
if not all(r['pass'] for r in reports):raise SystemExit('At least one lane failed pairing checks.')
