#!/usr/bin/env python3
"""Exact denominator-labelled envelopes, not an asymptotic claim about pi."""
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def farey(D):
    return sorted({F(p, d) for d in range(1, D+1) for p in range(d+1)})


def envelope(D, x):
    """Minimum d^2 distance to p/d; only nearest integer numerators needed."""
    out = []
    for d in range(1, D+1):
        y = d*x
        p = y.numerator//y.denominator
        for numerator in (p, p+1):
            out.append((d*d*abs(x-F(numerator,d)), F(numerator,d), d))
    return min(out)


def gap_segments(roots, u, v):
    # Every |x-r| is affine on this root-free gap. All pair crossings give
    # an independent finite partition on which the winning line is fixed.
    lines = [(F(r.denominator**2 if r <= u else -r.denominator**2),
              -F(r.denominator**2)*r if r <= u else F(r.denominator**2)*r,
              r) for r in roots]
    events = {u,v}
    for i,(s,b,_) in enumerate(lines):
        for ss,bb,_ in lines[i+1:]:
            if s != ss:
                cross = (bb-b)/(s-ss)
                if u < cross < v:
                    events.add(cross)
    result=[]
    cuts=sorted(events)
    for left,right in zip(cuts,cuts[1:]):
        mid=(left+right)/2
        slope,intercept,root=min(lines,key=lambda a:a[0]*mid+a[1])
        assert slope*mid+intercept == envelope(max(r.denominator for r in roots),mid)[0]
        if result and result[-1]['slope']==slope and result[-1]['intercept']==intercept:
            result[-1]['right']=right
        else:
            result.append(dict(left=left,right=right,slope=slope,intercept=intercept,root=root))
    return result


def atan_interval(inv, terms):
    partial=sum((F((-1)**k,(2*k+1)*inv**(2*k+1)) for k in range(terms)),F(0))
    next_term=F((-1)**terms,(2*terms+1)*inv**(2*terms+1))
    return min(partial,partial+next_term),max(partial,partial+next_term)


def pi_interval():
    # Machin identity; alternating-series remainder gives rational bounds.
    a,b=atan_interval(5,90)
    c,d=atan_interval(239,30)
    return 16*a-4*d,16*b-4*c


def best_pi(D):
    lo,hi=pi_interval()
    candidates=[]
    for d in range(1,D+1):
        mid=(lo+hi)*d/2
        p=(mid+F(1,2)).numerator//(mid+F(1,2)).denominator
        assert F(2*p-1,2) < d*lo < d*hi < F(2*p+1,2)
        r=F(p,d)
        assert r < lo or hi < r
        lower,upper=((d*d*(lo-r),d*d*(hi-r)) if r < lo
                     else (d*d*(r-hi),d*d*(r-lo)))
        candidates.append((lower,upper,p,d))
    winner=min(candidates,key=lambda a:a[1])
    assert all(winner[1]<other[0] for other in candidates if other != winner)
    return dict(D=D,p=winner[2],d=winner[3],lower=str(winner[0]),upper=str(winner[1]),
                display_midpoint=float((winner[0]+winner[1])/2),certified_unique_winner=True)


def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x


def main():
    gap_checks=0;segment_checks=0
    for D in range(1,9):
        roots=farey(D)
        for u,v in zip(roots,roots[1:]):
            a,b,c,d=u.numerator,u.denominator,v.numerator,v.denominator
            assert b*c-a*d==1 and b+d>D
            assert v-u==F(1,b*d)
            segments=gap_segments(roots,u,v)
            assert all(a['slope']>=b['slope'] for a,b in zip(segments,segments[1:]))
            gap_checks+=1;segment_checks+=len(segments)
    example=gap_segments(farey(3),F(1,3),F(1,2))
    assert [(r['left'],r['right'],r['slope'],r['intercept']) for r in example]==[
        (F(1,3),F(3,8),F(9),F(-3)),
        (F(3,8),F(2,5),F(1),F(0)),
        (F(2,5),F(1,2),F(-4),F(2))]
    endpoint_peak=F(6,13)
    actual_peak=max(envelope(3,x)[0] for x in [F(1,3),F(3,8),F(2,5),F(1,2)])
    assert actual_peak==F(2,5)<endpoint_peak
    rational_records=[]
    for x in (F(3),F(22,7)):
        for D in (1,6,7,113):
            value,r,d=envelope(D,x)
            rational_records.append(dict(x=x,D=D,value=value,winning_fraction=r,tagged_d=d))
    assert envelope(6,F(22,7))[0]>0 and envelope(7,F(22,7))[0]==0
    windows=[]
    for Q in (2,3,4):
        D=Q**Q
        points={F(p,d) for d in range(Q,D+1) for p in range(d+1)}
        tagged=sum(d+1 for d in range(Q,D+1))
        windows.append(dict(Q=Q,maximum_denominator=D,denominator_count=D-Q+1,
                            tagged_points_in_unit_interval=tagged,distinct_points=len(points)))
    alignment=[]
    for coordinates in ((F(1,2),F(2,3)),(F(2),F(-3))):
        B=lcm(*(x.denominator for x in coordinates))
        actual=[d for d in range(1,31) if all((d*x).denominator==1 for x in coordinates)]
        assert actual==[d for d in range(1,31) if d%B==0]
        alignment.append(dict(coordinates=coordinates,least_common_denominator=B,zeros_up_to_30=actual))
    payload=dict(status='passed',scope=__doc__,
                 farey_gap_checks=gap_checks,affine_segment_checks=segment_checks,
                 D3_nonlocal_integer_wins=example,
                 D3_actual_peak=actual_peak,D3_endpoint_only_peak=endpoint_peak,
                 rational_targets=rational_records,pi_targets=[best_pi(D) for D in (1,7,113)],
                 pi_interval_method='16 arctan(1/5)-4 arctan(1/239), alternating rational bounds',
                 denominator_windows=windows,shared_denominator_alignment=alignment,
                 asymptotic_pi_conclusion='none')
    destination=HERE/'denominator_envelope_results.json'
    destination.write_text(json.dumps(serial(payload),ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(status='passed',gap_checks=gap_checks,segments=segment_checks,
                          pi_winners=[(r['p'],r['d']) for r in payload['pi_targets']]),ensure_ascii=False))


if __name__=='__main__':main()
