"""Explicit network stage: authenticated file metadata for a pinned dataset revision."""
import json
from huggingface_hub import HfApi
from project_paths import project_root
R=project_root();m=json.loads((R/'sources/dataset_manifest.json').read_text())
info=HfApi().dataset_info(m['repo'],revision=m['revision'],files_metadata=True)
files=[]
for f in info.siblings:
    if f.rfilename in ('README.md','.gitattributes'):continue
    lfs=f.lfs
    files.append({'path':f.rfilename,'size':f.size,'blob_id':f.blob_id,'sha256':getattr(lfs,'sha256',None) if lfs else None})
(R/'data/authenticated_manifest.json').write_text(json.dumps({'repo':m['repo'],'revision':info.sha,'files':files},indent=2)+'\n')
print('Authenticated metadata saved locally; file count:',len(files))
