"""Download pinned challenge data; never print credentials or patient content."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
from datetime import datetime, timezone
from huggingface_hub import get_token, hf_hub_download

from project_paths import project_root, required_env
ROOT = project_root()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--core-only', action='store_true', help='Clinical document, VCF and index only')
    parser.add_argument('--inventory-only', action='store_true', help='Show file sizes without downloading')
    args = parser.parse_args()
    manifest = json.loads((ROOT / 'sources/dataset_manifest.json').read_text())
    files = [f for f in manifest['files'] if f['type'] == 'file' and f['path'] not in ('README.md', '.gitattributes')]
    if args.core_only:
        files = [f for f in files if not f['path'].endswith('.fastq.gz')]
    files.sort(key=lambda f: f['size'])
    dest = ROOT / 'data/raw'
    dest.mkdir(parents=True, exist_ok=True)
    print(f'Download inventory: {len(files)} files, {sum(f["size"] for f in files)} bytes', flush=True)
    if args.inventory_only:
        return 0
    if not get_token():
        print('Access required: accept the dataset conditions, then run .venv/bin/hf auth login locally. Do not paste tokens into chat.', file=sys.stderr)
        return 2
    remaining = sum(f['size'] for f in files if not (dest / f['path']).exists() or (dest / f['path']).stat().st_size != f['size'])
    if shutil.disk_usage(dest).free < remaining + 20 * 1024**3:
        print('Insufficient download headroom (remaining files plus 20 GiB).', file=sys.stderr)
        return 3
    receipt_path = ROOT / 'data/download_receipt.json'
    receipt = {'repo': manifest['repo'], 'revision': manifest['revision'], 'files': []}
    for f in files:
        try:
            path = Path(hf_hub_download(repo_id=manifest['repo'], repo_type='dataset', revision=manifest['revision'], filename=f['path'], local_dir=dest))
        except Exception as exc:
            # Avoid dumping request headers, signed URLs or credentials.
            print(f"Download failed: {type(exc).__name__}. Check access and network, then rerun to resume.", file=sys.stderr)
            return 4
        if path.stat().st_size != f['size']:
            print("Size mismatch in downloaded file", file=sys.stderr)
            return 5
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(8 * 1024**2), b''):
                digest.update(block)
        sha = digest.hexdigest()
        expected = f.get('lfs', {}).get('oid', '')
        upstream_available = len(expected) == 64 and all(c in '0123456789abcdef' for c in expected)
        if upstream_available and sha != expected:
            print("Upstream SHA-256 mismatch in downloaded file", file=sys.stderr)
            return 6
        receipt['files'].append({'path': f['path'], 'bytes': path.stat().st_size, 'sha256': sha, 'upstream_sha256_verified': upstream_available})
        receipt['updated_utc'] = datetime.now(timezone.utc).isoformat()
        receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
        print(f'Downloaded and size-checked: {len(receipt["files"])} files', flush=True)
    return 0

if __name__ == '__main__':
    sys.exit(main())
