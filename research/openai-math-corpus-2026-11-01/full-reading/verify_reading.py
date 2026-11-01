"""Check queue completeness and reject unsupported full-reading labels."""
import copy
import json
from pathlib import Path
import reading


def main():
    q=reading.queue()
    assert len(q)==721 and len({p['pdf_path'] for p in q})==721
    assert len({p['category'] for p in q[:17]})==17
    assert q[0]['family_id']=='001' and q[16]['family_id']=='362'
    assert len(json.loads((reading.ROOT/'supplementary_queue.json').read_text()))==30
    assert reading.interval_coverage([[2,4],[1,3],[6,6]],6)==(5,[[1,4],[6,6]])
    record=dict(position=1,pdf_path=q[0]['pdf_path'],status='full_body_read',
                source_files=[dict(path='test.tex',sha='a'*40,total_lines=6,read_ranges=[[1,6]])],
                source_mapping_complete=True,source_closure_complete=True,unread_files=[],sections_read=['test'],
                claim_as_stated='test',proof_architecture=['test'],construction_ideas=['test'],limits=['test'])
    result=reading.validate_record(record,{1:q[0]})
    assert result['read_source_lines']==6 and not result['independent_proof_verification']
    rejected=0
    for ranges,closed,unread in [([[1,5]],True,[]),([[1,6]],False,[]),([[1,6]],True,['appendix.tex'])]:
        bad=copy.deepcopy(record)
        bad['source_files'][0]['read_ranges']=ranges
        bad['source_closure_complete']=closed
        bad['unread_files']=unread
        try:reading.validate_record(bad,{1:q[0]})
        except ValueError:rejected+=1
    # Reading every declared line cannot justify "full" if include mapping
    # is still unfinished: an unlisted appendix would evade interval checks.
    bad=copy.deepcopy(record)
    bad.pop('source_mapping_complete')
    try:reading.validate_record(bad,{1:q[0]})
    except ValueError:rejected+=1
    assert rejected==4
    print(json.dumps(dict(status='passed',queue_size=721,first_round_disciplines=17,
                         supplementary_directories=30,unsupported_full_labels_rejected=4,
                         scope='Ledger consistency only; cannot mechanically prove that text was understood.'),indent=2))


if __name__=='__main__':main()
