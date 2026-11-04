#!/usr/bin/env python3
"""Exact small checks for the covering model, not a new pi record search."""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent


def floor(x):
    return x.numerator//x.denominator


def ceil(x):
    return -floor(-x)


def ball(p,q,c):
    return F(p,q)-c/(q*q),F(p,q)+c/(q*q)


def dangerous(j,c,q0,q1):
    """All reduced balls of q0<=q<q1 meeting the closed interval j."""
    lo,hi=j
    out=[]
    for q in range(q0,q1):
        for p in range(floor(lo*q-c/q)+1,ceil(hi*q+c/q)):
            if gcd(p,q)==1:
                left,right=ball(p,q,c)
                if left<hi and right>lo:
                    out.append((p,q,left,right))
    return out


def safe_children(j,c,Q):
    lo,hi=j
    ell=hi-lo
    assert ell==F(1,8*Q*Q)
    obstacles=dangerous(j,c,Q,2*Q)
    assert len(obstacles)<=1
    if obstacles:
        p,q,left,right=obstacles[0]
        parts=[]
        if lo<=left:
            parts.append((lo,min(left,hi)))
        if right<=hi:
            parts.append((max(right,lo),hi))
    else:
        parts=[j]
    unit=ell/4
    children=[]
    for left,right in parts:
        for k in range(floor((right-left)/unit)):
            child=(left+k*unit,left+(k+1)*unit)
            assert not dangerous(child,c,Q,2*Q)
            children.append(child)
            if len(children)==2:
                return children,len(obstacles)
    raise AssertionError('Two safe children must exist')


def atan_interval(inv,terms=32):
    partial=sum((F((-1)**k,(2*k+1)*inv**(2*k+1)) for k in range(terms)),F(0))
    following=partial+F((-1)**terms,(2*terms+1)*inv**(2*terms+1))
    return min(partial,following),max(partial,following)


def pi_interval():
    # tan(4 atan(1/5)-atan(1/239))=1; angle lies in (0,pi/2).
    t=F(1,5)
    t2=2*t/(1-t*t)
    t4=2*t2/(1-t2*t2)
    assert (t4-F(1,239))/(1+t4/F(239))==1
    a,b=atan_interval(5)
    d,e=atan_interval(239)
    return 16*a-4*e,16*b-4*d


def quality_bounds(p,q,pi_j):
    lo,hi=pi_j
    center=F(p,q)
    upper=q*q*max(abs(center-lo),abs(center-hi))
    lower=q*q*(lo-center if center<lo else center-hi if center>hi else 0)
    return lower,upper


def pair(x):
    return [x.numerator,x.denominator]


def run():
    result={}
    nodes=[(3,1,4,1)]
    checked=0
    prunes=0
    for depth in range(8):
        next_nodes=[]
        for a,b,r,s in nodes:
            assert b*r-a*s==1
            assert F(r,s)-F(a,b)==F(1,b*s)
            assert gcd(a+r,b+s)==1
            # Independently enumerate forbidden small denominators.
            for v in range(1,b+s):
                for u in range(floor(F(a*v,b))+1,ceil(F(r*v,s))):
                    assert not (F(a,b)<F(u,v)<F(r,s))
            c=F(1,10)
            if F(s,b)>=1/c:
                assert F(r,s)-F(a,b)<=c/(b*b)
                prunes+=1
            if F(b,s)>=1/c:
                assert F(r,s)-F(a,b)<=c/(s*s)
                prunes+=1
            next_nodes.extend([(a,b,a+r,b+s),(a+r,b+s,r,s)])
            checked+=1
        nodes=next_nodes
    # Include an explicit infinite-subtree pruning example.
    assert F(301,100)-3 <= F(1,100)
    result['farey_geometry']={'passed':True,'nodes_checked':checked,
                              'prunes_in_checked_nodes':prunes,
                              'explicit_pruned_interval':['3/1','301/100'],'c':'1/100'}

    c=F(1,64)
    leaves=[(F(25,8),F(13,4))]
    Q=1
    levels=[]
    for depth in range(7):
        new=[]
        hits=0
        for j in leaves:
            children,hit=safe_children(j,c,Q)
            new.extend(children)
            hits+=hit
        leaves=new
        Q*=2
        assert len(leaves)==2**(depth+1)
        assert all(j[1]-j[0]==F(1,8*Q*Q) for j in leaves)
        assert all(leaves[k][1]<=leaves[k+1][0] for k in range(len(leaves)-1))
        levels.append({'depth':depth+1,'leaves':len(leaves),'next_Q':Q,
                       'parents_with_dangerous_ball':hits})
    # Recheck final leaves from scratch against ALL prior denominator blocks.
    assert all(not dangerous(j,c,1,Q) for j in leaves)
    result['survival_tree']={'passed':True,'c':'1/64','levels':levels,
                             'all_q_below_final_Q_independently_checked':Q,
                             'is_a_subset_of_K_not_the_complete_complement':True}

    pj=pi_interval()
    assert F(31,10)<pj[0]<pj[1]<F(22,7)
    assert pj[1]-pj[0]<F(1,10**44)
    low,up=quality_bounds(355,113,pj)
    assert up<F(1,100)
    l,r=ball(355,113,F(1,100))
    assert l<pj[0]<pj[1]<r
    result['pi_cover_certificate']={'passed':True,'known_fraction':'355/113','c':'1/100',
                                    'pi_lower':pair(pj[0]),'pi_upper':pair(pj[1]),
                                    'quality_lower':pair(low),'quality_upper':pair(up),
                                    'not_a_new_discovery':True}

    # For c<1/2 every possible hit at denominator q has its center nearest q*pi.
    # We instead enumerate a slightly larger certified numerator range directly.
    threshold=F(1,1000)
    qcap=1000
    checked_q=0
    for q in range(1,qcap+1):
        p0=floor(pj[0]*q)-1
        p1=ceil(pj[1]*q)+1
        for p in range(p0,p1+1):
            if gcd(p,q)==1:
                left,right=ball(p,q,threshold)
                assert right<=pj[0] or left>=pj[1]
        checked_q+=1
    result['finite_nonhit']={'passed':True,'c':'1/1000','q_cap':qcap,
                             'denominators_checked':checked_q,
                             'no_infinite_tail_claim':True}
    return result


if __name__=='__main__':
    data=run()
    (ROOT/'verification.json').write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:{'passed':v['passed']} for k,v in data.items()},indent=2))
