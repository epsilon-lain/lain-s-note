#!/usr/bin/env python3
"""Exact finite scale-window comparisons and rational pi intervals.

The window [Q,Q**Q] is explicitly finite; no extrapolation about pi follows.
"""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json
from denominator_envelope_lab import pi_interval, envelope

HERE=Path(__file__).resolve().parent


def nearest(x,n,coprime=False):
    z=n*x+F(1,2)
    p=z.numerator//z.denominator
    if not coprime or gcd(p,n)==1:return p
    step=1
    while True:
        choices=[k for k in (p-step,p+step) if gcd(k,n)==1]
        if choices:return min(choices,key=lambda k:(abs(x-F(k,n)),k))
        step+=1


def pi_candidates(N,coprime=False):
    lo,hi=pi_interval()
    result=[]
    for n in range(1,N+1):
        p=nearest((lo+hi)/2,n,coprime)
        assert nearest(lo,n,coprime)==nearest(hi,n,coprime)==p
        a=F(p,n)
        assert a<lo or a>hi
        elo,ehi=((lo-a,hi-a) if a<lo else (a-hi,a-lo))
        result.append(dict(n=n,p=p,reduced_p=a.numerator,reduced_q=a.denominator,
                           lower=n*n*elo,upper=n*n*ehi))
    return result


def winner(rows):
    best=min(rows,key=lambda r:r['upper'])
    assert all(best['upper']<r['lower'] for r in rows if r is not best)
    return best


def serialize(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):serialize(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [serialize(v) for v in x]
    return x


def main():
    rows=pi_candidates(3125)
    reduced=pi_candidates(3125,True)
    windows=[]
    for Q in (2,3,4,5):
        end=Q**Q
        a=winner(rows[Q-1:end]);b=winner(reduced[Q-1:end])
        windows.append(dict(Q=Q,end=end,nominal_denominator_winner=a,
                            reduced_denominator_winner=b))
    assert windows[0]['nominal_denominator_winner']['reduced_q']==1
    assert windows[0]['reduced_denominator_winner']['n']==2
    targets=[]
    for x in (F(3),F(22,7)):
        record=[dict(n=n,p=nearest(x,n),distance=abs(x-F(nearest(x,n),n)),
                     kappa=n*n*abs(x-F(nearest(x,n),n))) for n in (1,6,7,8,113)]
        zeros=[n for n in range(1,257) if F(nearest(x,n),n)==x]
        assert zeros==[n for n in range(1,257) if n%x.denominator==0]
        targets.append(dict(x=x,samples=record,exact_grid_hits_up_to_256=zeros))
    # Wavefront arrival time and denominator envelope have the same sublevel sets.
    wave_checks=0
    for D in (1,3,7):
        for j in range(101):
            x=F(j,100)
            for c in (F(1,20),F(1,5),F(1,2)):
                inside=False
                for d in range(1,D+1):
                    p=nearest(x,d)
                    inside |= abs(x-F(p,d))<=c/(d*d)
                assert (envelope(D,x)[0]<=c)==inside
                wave_checks+=1
    lo,hi=pi_interval()
    # All values below are plotting midpoints, separate from exact interval evidence.
    display=[]
    best_upper=None
    for row in rows[:256]:
        if best_upper is None or row['upper']<best_upper['lower']:
            best_upper=row
        else:
            assert best_upper['upper']<row['lower']
        display.append(dict(n=row['n'],pi_kappa=float((row['lower']+row['upper'])/2),
                            pi_running_envelope=float((best_upper['lower']+best_upper['upper'])/2)))
    result=dict(status='passed',pi_lower=lo,pi_upper=hi,
                pi_interval_method='Machin identity with alternating rational bounds',
                pi_window_winners=windows,rational_grid_targets=targets,
                pi_samples=[rows[n-1] for n in (1,7,8,112,113,114)],
                wavefront_sublevel_checks=wave_checks,
                plotting_midpoints_not_interval_endpoints=display,
                bounds='Finite n<=3125 only; no tail or badly-approximable conclusion.')
    (HERE/'model_results.json').write_text(json.dumps(serialize(result),indent=2)+'\n')
    print(json.dumps(dict(status='passed',pi_observation_scales=3125,
                          coprime_scales=3125,wavefront_checks=wave_checks)))


if __name__=='__main__':main()
