"""Non-selective manuscript queue and evidence-based full-body reading ledger.

A coverage record reports source reading, not correctness of a theorem.
No relevance score is used when scheduling the queue.
"""
from pathlib import Path
from collections import Counter, deque
import argparse
import json
import re
from urllib.parse import quote

ROOT=Path(__file__).resolve().parent
CORPUS=ROOT.parent
PIN='4eaee87548250b2682affda7d05d947c05554bae'


def queue():
    papers=json.loads((CORPUS/'catalog.json').read_text())
    groups={}
    for paper in papers:
        groups.setdefault(paper['category'],deque()).append(paper)
    out=[]
    # All categories participate in every round until their queues are empty.
    # Preserve original catalogue order within each discipline.
    while any(groups.values()):
        for category, remaining in groups.items():
            if not remaining:continue
            p=remaining.popleft()
            out.append(dict(position=len(out)+1,family_id=p['family_id'],title=p['title'],
                            category=category,pdf_path=p['pdf_path'],directory=p['directory'],
                            directory_tree_sha=p['directory_tree_sha']))
    assert len(out)==len(papers)==721
    assert {p['pdf_path'] for p in out}=={p['pdf_path'] for p in papers}
    (ROOT/'queue.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    # Preserve uncatalogued versions as a separate supplementary queue.
    extras=json.loads((CORPUS/'uncatalogued.json').read_text())
    (ROOT/'supplementary_queue.json').write_text(json.dumps([
        dict(position=i+1,directory=p['directory'],directory_tree_sha=p['tree_sha'],
             status='pending_source_mapping') for i,p in enumerate(extras)],ensure_ascii=False,indent=2)+'\n')
    return out


def interval_coverage(ranges,total):
    """Validate inclusive 1-based intervals and merge overlap without double-count."""
    if not isinstance(total,int) or total<1:raise ValueError('positive total_lines required')
    checked=[]
    for pair in ranges:
        if len(pair)!=2 or any(not isinstance(n,int) for n in pair):
            raise ValueError('read ranges require two integer line numbers')
        a,b=pair
        if not 1<=a<=b<=total:raise ValueError(f'range {pair} outside1..{total}')
        checked.append((a,b))
    merged=[]
    for a,b in sorted(checked):
        if merged and a<=merged[-1][1]+1:merged[-1][1]=max(b,merged[-1][1])
        else:merged.append([a,b])
    covered=sum(b-a+1 for a,b in merged)
    return covered,merged


def validate_record(record,by_position):
    p=by_position[record['position']]
    if record['pdf_path']!=p['pdf_path']:raise ValueError('paper identity mismatch')
    if record['status'] not in {'full_body_read','partial_body_read','blocked'}:
        raise ValueError('unknown reading status')
    paths=set();coverage=[]
    for f in record.get('source_files',[]):
        path=f['path'] if f['path'].startswith('preprints/') else f"preprints/{p['directory']}/{f['path']}"
        if path in paths:raise ValueError('duplicate source path')
        paths.add(path)
        if not re.fullmatch(r'[0-9a-f]{40}',f['sha']):raise ValueError('missing source blob SHA')
        count,merged=interval_coverage(f.get('read_ranges',[]),f['total_lines'])
        coverage.append(dict(path=path,sha=f['sha'],total_lines=f['total_lines'],
                             source_url=f'https://github.com/openai/math/blob/{PIN}/'+quote(path,safe='/'),
                             read_lines=count,read_ranges=merged))
    total=sum(f['total_lines'] for f in coverage)
    read=sum(f['read_lines'] for f in coverage)
    full=(bool(coverage) and read==total and record.get('source_mapping_complete') is True
          and record.get('source_closure_complete') is True
          and not record.get('unread_files'))
    if record['status']=='full_body_read' and not full:
        raise ValueError(f"Unsupported full-body label at position{p['position']}")
    for field in ('sections_read','claim_as_stated','proof_architecture','construction_ideas','limits'):
        if not record.get(field):raise ValueError('missing '+field)
    normalized=dict(record,source_files=coverage,read_source_lines=read,listed_source_lines=total,
                    upstream_commit=PIN,independent_proof_verification=False)
    normalized.pop('source_cache',None)
    return normalized


def merge_continuation(record,extra,paper):
    """Union explicitly reported reading intervals; never promote mere downloads."""
    prefix=f"preprints/{paper['directory']}/"
    def canonical(path):return path if path.startswith('preprints/') else prefix+path
    sources={canonical(f['path']):f for f in record['source_files']}
    for f in extra['source_files']:
        key=canonical(f['path'])
        if key not in sources:raise ValueError('continuation outside declared source closure: '+key)
        original=sources[key]
        if original['sha']!=f['sha'] or original['total_lines']!=f['total_lines']:
            raise ValueError('continuation source version mismatch')
        _,ranges=interval_coverage(original.get('read_ranges',[])+f.get('read_ranges',[]),f['total_lines'])
        original['read_ranges']=ranges
    remaining=[]
    for key,f in sources.items():
        read,ranges=interval_coverage(f.get('read_ranges',[]),f['total_lines'])
        if read<f['total_lines']:
            start=1;missing=[]
            for a,b in ranges:
                if start<a:missing.append([start,a-1])
                start=b+1
            if start<=f['total_lines']:missing.append([start,f['total_lines']])
            remaining.append(dict(path=key,sha=f['sha'],unread_ranges=missing))
    # Unknown include paths in the primary receipt cannot disappear by merging.
    for f in record.get('unread_files',[]):
        path=f['path'] if isinstance(f,dict) else f
        if canonical(path) not in sources:remaining.append(f)
    record['unread_files']=remaining
    record.setdefault('continuation_receipts',[]).append(extra)
    if not remaining and record.get('source_mapping_complete') is True:
        record['source_closure_complete']=True
        record['status']='full_body_read'
        record['scope_note']=record.get('scope_note','')+' Body coverage assembled from explicitly identified readers; not independent theorem verification.'
    return record


def build():
    q=queue();by_position={p['position']:p for p in q}
    raw_records={}
    for path in sorted(ROOT.glob('batch_*/records.json')):
        raw=json.loads(path.read_text())
        if isinstance(raw,dict):raw=raw['records']
        for record in raw:
            if record['position'] in raw_records:raise ValueError('duplicate manuscript record')
            raw_records[record['position']]=record
    for path in sorted(ROOT.glob('continuation_*.json')):
        extra=json.loads(path.read_text())
        extra['position']=extra.get('position',extra.get('paper_position'))
        if extra['position'] not in raw_records:raise ValueError('continuation has no base receipt')
        extra['receipt_file']=path.name
        key=extra['position']
        raw_records[key]=merge_continuation(raw_records[key],extra,by_position[key])
    records={key:validate_record(record,by_position) for key,record in raw_records.items()}
    ledger=[]
    for p in q:
        record=records.get(p['position'])
        ledger.append(dict(p,status=record['status'] if record else 'not_started',
                           record=record))
    counts=Counter(p['status'] for p in ledger)
    summary=dict(upstream_commit=PIN,catalogued_manuscripts=721,
                 supplementary_uncatalogued_directories=30,
                 status_counts=dict(counts),
                 source_lines_read=sum(r['read_source_lines'] for r in records.values()),
                 source_lines_listed=sum(r['listed_source_lines'] for r in records.values()),
                 disciplines_in_queue=len({p['category'] for p in q}),
                 disciplines_with_body_reading=len({by_position[k]['category'] for k,r in records.items()
                                                    if r['read_source_lines']>0}),
                 ordering='round_robin_disciplines_then_catalogue_order_no_relevance_filter',
                 scope='Reported complete source reading does not certify theorem correctness or cited dependencies.')
    mode_file=ROOT/'reading_mode.json'
    if mode_file.exists():
        # Queue order and the user's current selection policy are distinct.
        summary['reading_mode']=json.loads(mode_file.read_text())
    (ROOT/'ledger.json').write_text(json.dumps(ledger,ensure_ascii=False,indent=2)+'\n')
    (ROOT/'progress.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    lines=['# 全库正文阅读台账','',
           f"目录稿件721；全文首读{counts['full_body_read']}；部分正文{counts['partial_body_read']}；"
           f"受阻{counts['blocked']}；未开始{counts['not_started']}。",'',
           '完整首读只表示已读相应正文与附录，不表示论文主张已独立证明。', '',
           '| 顺序 | 组 | 学科 | 稿件 | 状态 | 正文源码行数（已读／列入记录） |',
           '|---|---|---|---|---|---|']
    labels={'full_body_read':'正文首读完成','partial_body_read':'部分正文，待续',
            'blocked':'受阻，保留任务','not_started':'未开始'}
    for p in ledger:
        r=p['record'];n=f"{r['read_source_lines']}/{r['listed_source_lines']}" if r else '—'
        title=p['title'].replace('|','\\|')
        url=f"https://github.com/openai/math/blob/{PIN}/"+quote(p['pdf_path'],safe='/')
        lines.append(f"| {p['position']} | {p['family_id']} | {p['category']} | [{title}]({url}) | {labels[p['status']]} | {n} |")
    (ROOT/'LEDGER.md').write_text('\n'.join(lines)+'\n')
    ideas=[]
    note_lines=['# 跨领域构造线索：按阅读顺序保留','',
                '这里记录正文首读中的构造机制与可能迁移的接口。候选思路不等于已证明能用于 π；不因暂时没有接口而从全库队列排除论文。','']
    for p in ledger:
        r=p['record']
        if not r:continue
        all_ideas=list(r['construction_ideas'])
        for receipt in r.get('continuation_receipts',[]):
            all_ideas.extend(receipt.get('construction_ideas',[]))
        ideas.append(dict(position=p['position'],family_id=p['family_id'],title=p['title'],
                          reading_status=p['status'],construction_ideas=all_ideas))
        note_lines.extend([f"## {p['position']}. {p['title']}（{p['family_id']}）",'',
                           f"学科：{p['category']}；阅读状态：{labels[p['status']]}。",''])
        for idea in all_ideas:
            note_lines.append('- '+(idea if isinstance(idea,str) else json.dumps(idea,ensure_ascii=False)))
        note_lines.append('')
    (ROOT/'ideas.json').write_text(json.dumps(ideas,ensure_ascii=False,indent=2)+'\n')
    (ROOT/'IDEAS.md').write_text('\n'.join(note_lines)+'\n')
    print(json.dumps(summary,ensure_ascii=False,indent=2))


def pending(limit):
    ledger=json.loads((ROOT/'ledger.json').read_text())
    remaining=[p for p in ledger if p['status']!='full_body_read']
    print(json.dumps([dict(position=p['position'],family_id=p['family_id'],
                           title=p['title'],status=p['status'],
                           continuation=p['record'].get('unread_files',[]) if p['record'] else None)
                      for p in remaining[:limit]],ensure_ascii=False,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['build','pending'])
    parser.add_argument('--limit',type=int,default=20)
    args=parser.parse_args()
    if args.command=='build':build()
    else:pending(args.limit)
