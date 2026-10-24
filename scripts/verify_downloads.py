"""Validate completed downloads without displaying patient data."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time
import zipfile

from project_paths import project_root, required_env
ROOT = project_root()
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--watch', action='store_true', help='Wait for pending downloads, up to 24 hours')
args = parser.parse_args()
manifest = json.loads((ROOT/'data/authenticated_manifest.json').read_text())
report_path = ROOT/'data/integrity_report.json'
report = {'revision':manifest['revision'], 'files':[], 'complete':False}
# Invalidate previous success before any verification can fail.
temp = report_path.with_suffix('.tmp')
temp.write_text(json.dumps(report,indent=2)+'\n')
temp.replace(report_path)
verified = set()
deadline = time.monotonic() + 24*3600
while True:
    try:
        receipt = json.loads((ROOT/'data/download_receipt.json').read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        receipt = {'files':[]}
    if receipt['files'] and any(receipt.get(key) != manifest.get(key) for key in ('repo', 'revision')):
        raise ValueError('Download receipt repository/revision differs from authenticated manifest')
    completed = {f['path']:f for f in receipt['files']}
    for meta in manifest['files']:
        name = meta['path']
        if name in verified or name not in completed:
            continue
        path = ROOT/'data/raw'/name
        if path.stat().st_size != meta['size']:raise ValueError('File size mismatch')
        entry = {'path':name, 'size_verified':True}
        if meta['sha256']:
            digest = hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda: stream.read(8*1024**2), b''):
                    digest.update(block)
            if digest.hexdigest() != meta['sha256']:
                raise ValueError('Upstream SHA-256 mismatch')
            entry['upstream_sha256_verified'] = True
        else:
            digest = hashlib.sha1()
            digest.update(f'blob {path.stat().st_size}\0'.encode())
            with path.open('rb') as stream:
                for block in iter(lambda: stream.read(8*1024**2), b''):
                    digest.update(block)
            if digest.hexdigest() != meta['blob_id']:raise ValueError('Upstream Git blob mismatch')
            entry['upstream_git_blob_verified'] = True
        if name.endswith('.gz'):
            with gzip.open(path,'rb') as stream:
                while stream.read(8*1024**2):
                    pass
            entry['gzip_crc_verified'] = True
        elif name.endswith('.docx'):
            with zipfile.ZipFile(path) as archive:
                if archive.testzip() is not None:raise ValueError('DOCX ZIP integrity error')
            entry['zip_crc_verified'] = True
        report['files'].append(entry)
        verified.add(name)
        report['complete'] = len(verified) == len(manifest['files'])
        temp = report_path.with_suffix('.tmp')
        temp.write_text(json.dumps(report,indent=2)+'\n')
        temp.replace(report_path)
        print(f'Integrity verified: {len(verified)}/{len(manifest["files"])} files',flush=True)
    if report['complete'] or not args.watch:
        break
    if time.monotonic() > deadline:
        raise SystemExit('Verification timed out while waiting for downloads; inspect receipt and rerun.')
    time.sleep(10)
print('All files verified.' if report['complete'] else 'Pending files remain.',flush=True)

raise SystemExit(0 if report['complete'] else 2)
