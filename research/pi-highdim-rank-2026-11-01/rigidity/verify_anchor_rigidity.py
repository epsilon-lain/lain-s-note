"""Exact finite support/matching checks for ANCHOR_BLOCK_RIGIDITY.md."""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import json


def multiply(p,q,N):
    out=[F(0)]*(N+1)
    for i,x in enumerate(p):
        for j,y in enumerate(q[:N+1-i]):out[i+j]+=x*y
    return out


def tails_product(T,d,N):
    out=[F(1)]+[F(0)]*N
    for Ti,di in zip(T,d):
        tail=[F((-1)**k,k) if k>=Ti and k else F(0) for k in range(N+1)]
        for _ in range(di):out=multiply(out,tail,N)
    return out


def injection_count(edges):
    # Degree ordering makes the forced identity chain inexpensive to enumerate.
    order=sorted(range(len(edges)),key=lambda i:len(edges[i][1]))
    count=nonidentity=0
    def visit(n,used,moved):
        nonlocal count,nonidentity
        if n==len(order):
            count+=1;nonidentity+=moved;return
        i=order[n]
        for target in edges[i][1]:
            if target not in used:
                visit(n+1,used|{target},moved or target!=edges[i][0])
    visit(0,set(),False)
    return count,nonidentity


def check(H,w,T,theta=F(1,2),v0=F(1),expect_rigid=True):
    m=len(w)
    rows=[]
    ranges=[range(int(theta*H/wi)+1) for wi in w]
    for beta in product(*ranges):
        s=0
        while v0*s+sum(wi*bi for wi,bi in zip(w,beta))/theta<H:
            rows.append((s,beta));s+=1
    rowset=set(rows)
    edges=[];nonidentity_edges=outside=0
    strict=all(v0*Ti>wi/theta for wi,Ti in zip(w,T))
    for source in rows:
        s,beta=source
        target_set=set()
        for d in product(*(range(s//Ti+1) for Ti in T)):
            a=tuple(bi+di for bi,di in zip(beta,d))
            if sum(wi*ai for wi,ai in zip(w,a))>H:continue
            poly=tails_product(T,d,s)
            for k,coefficient in enumerate(poly):
                if not coefficient:continue
                ell=s-k
                target=(ell,a)
                scalar=coefficient
                for ai,bi in zip(a,beta):scalar*=comb(ai,bi)
                assert scalar
                target_set.add(target)
                if target!=source:
                    nonidentity_edges+=1
                    if strict:
                        source_weight=v0*s+sum(wi*bi for wi,bi in zip(w,beta))/theta
                        target_weight=v0*ell+sum(wi*ai for wi,ai in zip(w,a))/theta
                        assert target_weight<source_weight
                if target not in rowset:outside+=1
        assert source in target_set
        edges.append((source,target_set))
    count,moved=injection_count(edges)
    if expect_rigid:
        assert strict and outside==0 and count==1 and moved==0
    else:
        assert not strict and outside>0 and moved>0
    return {'H':H,'weights':list(w),'tail_orders':list(T),
            'theta':str(theta),'v0':str(v0),'source_rows':len(rows),
            'nonidentity_nonzero_edges':nonidentity_edges,
            'edges_outside_source_packet':outside,
            'injective_selections':count,'nonidentity_injective_selections':moved,
            'strict_filtration_condition':strict}


def main():
    cases=[check(5,(1,),(3,)),check(8,(1,),(3,)),
           check(6,(1,2),(3,5))]
    failure=check(2,(1,),(1,),expect_rigid=False)
    report={'status':'passed','arithmetic':'exact Fraction logarithmic-tail coefficients',
            'rigid_cases':cases,'failure_without_packet_preservation':failure,
            'scope':'Finite support checks; universal rigidity follows from strict weight descent and finite injectivity. No pi approximation conclusion.'}
    Path(__file__).with_name('anchor_rigidity_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
