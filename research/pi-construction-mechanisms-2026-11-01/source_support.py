"""Optional exact source audit, with cache or pinned git submodule support."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent


def audit_source_files(paths=None):
    manifest=json.loads((HERE/'pinned_source_manifest.json').read_text())
    rows=manifest['source_files']
    if paths is not None:
        rows=[r for r in rows if r['path'] in paths]
        assert len(rows)==len(paths)
    candidates=[ROOT/'source_context/pi017',
                ROOT/'external/openai-math'/manifest['directory']]
    verified=[];missing=[]
    for r in rows:
        source=next((base/r['path'] for base in candidates if (base/r['path']).is_file()),None)
        if source is None:
            missing.append(r['path']);continue
        b=source.read_bytes()
        observed=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
        assert observed==r['sha'],('source version mismatch',r['path'])
        assert len(b)==r['size'] and len(b.splitlines())==r['total_lines']
        verified.append(dict(path=r['path'],git_blob_sha=observed,size=len(b),
                             lines=r['total_lines'],matches_pinned_closure=True))
    return dict(pin=manifest['upstream_commit'],
                status='passed' if not missing else ('partial' if verified else 'source_unavailable'),
                records=verified,missing=missing,
                scope='Byte identity only; no theorem or Lean verification.')
