"""Verify reading-receipt blobs against a pinned checkout or local source caches.

This checks source identity and line counts, not attention or proof correctness.
"""
from pathlib import Path
from collections import Counter
import argparse
import hashlib
import json

ROOT=Path(__file__).resolve().parent


def blob(data):return hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkout',type=Path)
    p.add_argument('--cache-root',action='append',type=Path,default=[])
    args=p.parse_args()
    if not args.checkout and not args.cache_root:p.error('pass --checkout or at least one --cache-root')
    index={}
    for root in args.cache_root:
        for path in root.rglob('*'):
            if path.is_file() and path.suffix in {'.tex','.bib','.md'}:
                data=path.read_bytes();index.setdefault(blob(data),[]).append(data)
    records=json.loads((ROOT/'ledger.json').read_text())
    counts=Counter();failures=[]
    for entry in records:
        record=entry['record']
        if not record:continue
        for source in record['source_files']:
            if args.checkout:
                path=args.checkout/source['path']
                candidates=[path.read_bytes()] if path.exists() else []
            else:candidates=index.get(source['sha'],[])
            good=any(blob(data)==source['sha'] and len(data.decode('utf-8').splitlines())==source['total_lines']
                     for data in candidates)
            counts['source_files']+=1
            counts['passed' if good else 'failed']+=1
            if not good:failures.append(dict(position=entry['position'],path=source['path'],sha=source['sha']))
    result=dict(status='passed' if not failures else 'failed',counts=dict(counts),failures=failures,
                scope='Exact source-byte and line-count audit; no independent theorem verification.')
    print(json.dumps(result,ensure_ascii=False,indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
