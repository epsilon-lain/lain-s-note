"""Integrity checks against independent snapshot counts and edge-case paths."""
from pathlib import Path
import contextlib
import hashlib
import io
import json
import re
import corpus


def main():
    root=Path(__file__).resolve().parent
    papers=json.loads((root/'catalog.json').read_text())
    families=json.loads((root/'families.json').read_text())
    manifest=json.loads((root/'manifest.json').read_text())
    inventory=json.loads((root/'snapshots/preprint_inventory.json').read_text())
    extras=json.loads((root/'uncatalogued.json').read_text())
    raw=(root/'snapshots/CONTENTS.md').read_text()
    assert len(papers)==raw.count('&emsp;[')==721
    assert len(families)==len(re.findall(r'^\*\*\d{3}\. ',raw,re.M))==373
    assert len(inventory['directories'])==751 and len(extras)==30
    assert len({p['pdf_path'] for p in papers})==721
    dirs={d['directory'] for d in inventory['directories']}
    mapped={p['directory'] for p in papers}
    assert mapped<=dirs and mapped.isdisjoint(d['directory'] for d in extras)
    assert mapped|{d['directory'] for d in extras}==dirs
    # The parser must handle parentheses inside a path and a suffix after link.
    assert any('CAT(0)' in p['pdf_path'] and p['pdf_path'].endswith('.pdf') for p in papers)
    secondary=[p for p in papers if p['upstream_note']=='— secondary writeup']
    assert len(secondary)==1 and 'Matrix Multiplication' in secondary[0]['title']
    assert all(p['abstract'] and p['category'] and corpus.PIN in p['pdf_url'] for p in papers)
    assert sum(manifest['category_counts'].values())==721
    for name,digest in manifest['snapshot_sha256'].items():
        assert hashlib.sha256((root/'snapshots'/name).read_bytes()).hexdigest()==digest
    out=io.StringIO()
    with contextlib.redirect_stdout(out):corpus.search('',family='17',as_json=True)
    query=json.loads(out.getvalue())
    assert query['total']==1 and query['results'][0]['family_id']=='017'
    out=io.StringIO()
    with contextlib.redirect_stdout(out):corpus.search('matroid',as_json=True,limit=1000)
    query=json.loads(out.getvalue())
    assert query['total']>0 and any(p['family_id']=='114' for p in query['results'])
    print(json.dumps(dict(status='passed',manuscripts=721,families=373,
                         preprint_directories=751,uncatalogued_directories=30,
                         snapshot_hashes_checked=len(manifest['snapshot_sha256']),
                         edge_case_pdf_paths_checked=2,offline_search_queries=2,
                         scope='Metadata and exact snapshot integrity. Full submodule download/Lean proofs were not run here.'),indent=2))


if __name__=='__main__':main()
