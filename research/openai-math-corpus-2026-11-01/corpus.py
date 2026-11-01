"""Pinned OpenAI math catalogue builder/search tool. Python standard library.

Metadata search is offline. Full source/PDF access is supplied by the pinned
external/openai-math Git submodule, initialized explicitly by the user.
"""
from pathlib import Path
import argparse
import collections
import hashlib
import html
import json
import re
import subprocess
import sys
import urllib.parse

ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parents[1]
PIN = '4eaee87548250b2682affda7d05d947c05554bae'
UPSTREAM = 'https://github.com/openai/math'
KEYWORDS = {
    'diophantine': ['diophantine', 'irrational', 'approximation', 'duffin', 'continued fraction', 'transcenden'],
    'determinant_rank': ['determinant', 'determinantal', 'rank', 'interpolation', 'multiplicity', 'zero test', 'identity test'],
    'matroid_optimization': ['matroid', 'common bases', 'entropy', 'log-concav', 'sampling', 'optimization'],
    'dynamics': ['dynamical', 'equidistribution', 'ergodic', 'shrinking target', 'recurrence'],
    'arithmetic_sum': ['exponential sum', 'gauss sum', 'lattice point', 'divisor', 'integer point'],
}


def plain(text):
    # Preserve the upstream raw abstract separately; this is a search rendering.
    text = re.sub(r'<(?:sub|sup)>(.*?)</(?:sub|sup)>', r'{\1}', text)
    return ' '.join(html.unescape(re.sub(r'<[^>]*>', '', text)).split())


def source_url(path):
    return f'{UPSTREAM}/blob/{PIN}/' + urllib.parse.quote(path, safe='/')


def parse_catalog(contents, overview, formalization, inventory):
    categories, category = {}, None
    for line in overview.splitlines():
        match = re.match(r'\\cataloguesection\{([^}]+)\}', line)
        if match:
            category = match[1]
        match = re.match(r'\\resultentry\{(\d+)\}', line)
        if match:
            categories[match[1].zfill(3)] = category
    formal_paths = set()
    for line in formalization.splitlines():
        match = re.match(r'\s+id:\s*(.+)', line)
        if match:
            value = match[1].strip().strip('\"\'')
            if value.startswith('../preprints/'):
                formal_paths.add(value[3:])
    headers = list(re.finditer(r'^\*\*(\d{3})\. (.*?)\*\*([^\n]*)', contents, re.M))
    papers, families = [], []
    directory_map = {d['directory']: d for d in inventory['directories']}
    for i, header in enumerate(headers):
        family_id, title, description = header.groups()
        block = contents[header.end():headers[i+1].start() if i+1<len(headers) else len(contents)]
        matches = list(re.finditer(r'^&emsp;\[(.*?)\]\((preprints/[^\n]+?\.pdf)\)([^\n]*)', block, re.M))
        family = dict(id=family_id, title=plain(title), description=plain(description),
                      category=categories.get(family_id), manuscripts=len(matches),
                      lean_scope_link=(source_url(f'lean/docs/{family_id}.md')
                                       if f'lean/docs/{family_id}.md' in description else None))
        families.append(family)
        for match in matches:
            paper_title, path, note = match.groups()
            remainder = block[match.end():]
            abstract = remainder.split('</td>', 1)[0].strip()
            directory = path.split('/')[1]
            if directory not in directory_map:
                raise ValueError(f'Catalog path directory missing from tree: {path}')
            haystack = plain(paper_title + ' ' + abstract + ' ' + description).casefold()
            tags = [tag for tag, terms in KEYWORDS.items() if any(t in haystack for t in terms)]
            papers.append(dict(family_id=family_id, family_title=family['title'],
                               category=family['category'], title=plain(paper_title),
                               abstract=plain(abstract), abstract_upstream_markup=abstract,
                               family_claim_description=family['description'],
                               pdf_path=path, pdf_url=source_url(path),
                               directory=directory, directory_tree_sha=directory_map[directory]['tree_sha'],
                               directory_url=f'{UPSTREAM}/tree/{PIN}/preprints/'+urllib.parse.quote(directory,safe=''),
                               formalization_source_listed=path in formal_paths,
                               family_lean_scope_url=family['lean_scope_link'],
                               upstream_note=plain(note), screening_tags=tags,
                               review_status='catalogue_metadata_only'))
    if len({p['pdf_path'] for p in papers}) != len(papers):
        raise ValueError('Duplicate manuscript PDF paths')
    declared = re.search(r'\*\*(\d+) manuscripts covering (\d+) result families', contents)
    if not declared or (len(papers),len(families)) != tuple(map(int,declared.groups())):
        raise ValueError('Parsed counts do not match upstream declaration')
    if any(f['category'] is None or f['manuscripts']==0 for f in families):
        raise ValueError('Missing category or constituent manuscripts')
    return papers, families


def build():
    snapshots = ROOT / 'snapshots'
    inventory = json.loads((snapshots/'preprint_inventory.json').read_text())
    papers, families = parse_catalog(
        (snapshots/'CONTENTS.md').read_text(), (snapshots/'overview.tex').read_text(),
        (snapshots/'formalization.yaml').read_text(), inventory)
    mapped = {p['directory'] for p in papers}
    extras = [d for d in inventory['directories'] if d['directory'] not in mapped]
    manifest = dict(schema_version=1, upstream_repository=UPSTREAM, upstream_commit=PIN,
                    upstream_license='Apache-2.0', user_count_reference=772,
                    catalogued_manuscripts=len(papers), result_families=len(families),
                    preprint_directories=len(inventory['directories']),
                    uncatalogued_directories=len(extras),
                    formalization_source_listed=sum(p['formalization_source_listed'] for p in papers),
                    families_with_lean_scope=sum(bool(f['lean_scope_link']) for f in families),
                    snapshot_date='2026-11-01',
                    verification_note='Upstream catalogue and formalization links are metadata, not independent proof verification.',
                    snapshot_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                     for p in sorted(snapshots.iterdir()) if p.is_file()},
                    category_counts=dict(sorted(collections.Counter(p['category'] for p in papers).items())),
                    screening_tag_counts={tag:sum(tag in p['screening_tags'] for p in papers) for tag in KEYWORDS})
    for filename, data in [('catalog.json',papers),('families.json',families),
                           ('uncatalogued.json',extras),('manifest.json',manifest)]:
        (ROOT/filename).write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
    lines = ['# 固定版本论文索引','',f'OpenAI math `{PIN}`：{len(papers)} 篇目录稿件、{len(families)} 个主题组。',
             '', '条目来自上游标题和摘要；列入此表不表示已独立验证证明。完整全文与源文件见仓库根目录的 `external/openai-math` 子模块。',
             '', '| 组 | 学科 | 稿件 | 上游形式化目录列出 |', '|---|---|---|---|']
    for p in papers:
        title=p['title'].replace('|','\\|').replace('[','\\[').replace(']','\\]')
        lines.append(f"| {p['family_id']} | {p['category']} | [{title}]({p['pdf_url']}) | {'是（需读范围）' if p['formalization_source_listed'] else '否'} |")
    (ROOT/'CATALOG.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(manifest,ensure_ascii=False,indent=2))


def search(query, family=None, tag=None, limit=20, as_json=False):
    papers=json.loads((ROOT/'catalog.json').read_text())
    terms=query.casefold().split()
    found=[]
    for p in papers:
        if family and p['family_id'] != family.zfill(3): continue
        if tag and tag not in p['screening_tags']: continue
        hay=' '.join([p['title'],p['abstract'],p['family_claim_description'],p['category']]).casefold()
        if all(term in hay for term in terms):found.append(p)
    if as_json:
        print(json.dumps({'total':len(found),'results':found[:limit]},ensure_ascii=False,indent=2))
    else:
        print(f'{len(found)} matches (showing at most {limit})')
        for p in found[:limit]:
            print(f"[{p['family_id']}] {p['title']}\n{p['pdf_url']}\n{p['abstract']}\n")


def fulltext(query):
    checkout=REPOSITORY/'external/openai-math'
    if not (checkout/'CONTENTS.md').exists():
        raise SystemExit('Initialize the pinned corpus first: git submodule update --init external/openai-math')
    pin=subprocess.check_output(['git','-C',str(checkout),'rev-parse','HEAD'],text=True).strip()
    if pin != PIN: raise SystemExit('Submodule revision differs from catalog pin; rebuild deliberately before using it.')
    try:
        result=subprocess.run(['rg','--line-number','--ignore-case','--fixed-strings','--max-count','3',
                               '--glob','*.tex','--glob','*.md','--glob','*.lean',
                               '--',query,str(checkout)],check=False)
    except FileNotFoundError:
        raise SystemExit('Full-text search requires ripgrep (rg); metadata search needs only Python.')
    if result.returncode not in (0,1):raise SystemExit(result.returncode)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    commands.add_parser('build')
    q=commands.add_parser('search');q.add_argument('query',nargs='?',default='')
    q.add_argument('--family');q.add_argument('--tag',choices=KEYWORDS);q.add_argument('--limit',type=int,default=20)
    q.add_argument('--json',action='store_true')
    q=commands.add_parser('fulltext');q.add_argument('query')
    args=parser.parse_args()
    if args.command=='build':build()
    elif args.command=='search':search(args.query,args.family,args.tag,args.limit,args.json)
    else:fulltext(args.query)
