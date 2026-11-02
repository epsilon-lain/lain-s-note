"""Replay all saved positive grid certificates with an independent lift formula.
Does not rerun or rely on the incomplete float proposal search.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
from verify_gain import true_lift,pi_interval,quality
ROOT=Path(__file__).resolve().parents[1]

def run():
    saved=json.loads((ROOT/'experiment'/'certificates.json').read_text())
    summary=json.loads((ROOT/'experiment'/'results.json').read_text())
    P=pi_interval();seen=set();steps=0;verified=[]
    for c in saved:
        om,kap=Q(c['omega']),Q(c['kappa']);q,p=c['iterate_q'],c['lift_p']
        j=kap*5000;assert j.denominator==1 and 1<=j<=600 and j not in seen;seen.add(int(j))
        x,r,A=map(Q,(c['periodic_point'],c['trapping_radius'],c['q_step_multiplier']))
        assert 0<A<1 and r>0 and (1-A)*r==Q(c['strict_margin'])
        outputs=[]
        for z in (x-r,x,x+r):
            for _ in range(q):z=true_lift(z,om,kap);steps+=1
            outputs.append(z-p)
        assert outputs==[x-A*r,x,x+A*r]
        assert x-r<outputs[0]<x<outputs[-1]<x+r
        rotation=Q(p,q);result=3+rotation
        assert str(rotation)==c['rotation'] and str(result)==c['output']
        assert result.denominator==c['actual_reduced_denominator']
        lower,upper=quality(result,P);stored_lo,stored_hi=map(Q,c['quality_interval'])
        assert stored_lo<=lower<=upper<=stored_hi
        verified.append((upper,result,j))
    unknown=set(summary['unresolved_grid_indices'])
    assert not seen&unknown and seen|unknown==set(range(1,601))
    assert len(saved)==summary['positive_certificates']==191
    best=min(verified)
    assert best[1]==Q(355,113) and int(best[2])==493
    out={'status':'PASS','exact_positive_certificates_replayed':len(saved),'absolute_lift_evaluations':steps,
         'attempted_grid_size':600,'unresolved_count':len(unknown),'best_among_certified':str(best[1]),
         'winning_grid_index':int(best[2]),'scope':'All saved positive certificates independently replayed; unresolved grid points are not exclusion results.'}
    (Path(__file__).with_name('saved_grid_checks.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':run()
