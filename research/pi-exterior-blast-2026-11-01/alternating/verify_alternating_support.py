"""Exact checks of alternating boundary support; standard library only."""
from fractions import Fraction as F
from itertools import product, combinations
from math import comb, ceil, prod
from pathlib import Path
import json


def total(c,k):
    return sum(a*b for a,b in zip(c,k))


def alternating(c,g,d):
    return sum((-1)**sum(k)*prod(comb(n,j) for n,j in zip(g,k))
               for k in product(*(range(n+1) for n in g)) if total(c,k)<d)


def slab(c,g,d,i):
    reduced=list(g);reduced[i]-=1
    return sum((-1)**sum(k)*prod(comb(n,j) for n,j in zip(reduced,k))
               for k in product(*(range(n+1) for n in reduced))
               if total(c,k)<d<=total(c,k)+c[i])


def polynomial(c,g):
    p=[1]
    for cost,mult in zip(c,g):
        q=[0]*(cost*mult+1)
        for k in range(mult+1):q[cost*k]=(-1)**k*comb(mult,k)
        out=[0]*(len(p)+len(q)-1)
        for i,x in enumerate(p):
            for j,y in enumerate(q):out[i+j]+=x*y
        p=out
    assert sum(p)==0
    acc=0;out=[]
    for value in p[:-1]:
        acc+=value;out.append(acc)
    return out



def determinant(a):
    a=[[F(x) for x in row] for row in a];value=F(1)
    for k in range(len(a)):
        pivot=next((i for i in range(k,len(a)) if a[i][k]),None)
        if pivot is None:return F(0)
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];value=-value
        x=a[k][k];value*=x
        for i in range(k+1,len(a)):
            ratio=a[i][k]/x
            for j in range(k+1,len(a)):a[i][j]-=ratio*a[k][j]
    return value


def vandermonde(xs):
    return prod(xs[j]-xs[i] for i in range(len(xs)) for j in range(i+1,len(xs)))


def boundary_entry(r,b,a,zeta):
    # Independently form U^{-1} T_{L,O}, using binomial transport.
    return sum(comb(k,b)*zeta**(k-b)*comb(a,k)*(-zeta)**(a-k)
               for k in range(b,r+1))


def cauchy_checks():
    entries=0;minors=0;full=0;units=[]
    for n in range(1,9):
        for r in range(n):
            for zeta in (F(1,5),F(-2,7)):
                for b in range(r+1):
                    for a in range(r+1,n+1):
                        direct=boundary_entry(r,b,a,zeta)
                        formula=(-1)**(a+r)*zeta**(a-b)*F((r+1)*comb(r,b)*comb(a,r+1),a-b)
                        assert direct==formula
                        entries+=1
                for k in range(1,min(r+1,n-r,3)+1):
                    for bs in combinations(range(r+1),k):
                        for aa in combinations(range(r+1,n+1),k):
                            direct=determinant([[boundary_entry(r,b,a,zeta) for a in aa] for b in bs])
                            formula=(-1)**(sum(aa)+k*r+k*(k-1)//2)*zeta**(sum(aa)-sum(bs))*(r+1)**k
                            formula*=prod(comb(r,b) for b in bs)*prod(comb(a,r+1) for a in aa)
                            formula*=F(vandermonde(bs)*vandermonde(aa),prod(a-b for b in bs for a in aa))
                            assert direct==formula and direct
                            minors+=1
                if n-r>=r+1:
                    for aa in combinations(range(r+1,n+1),r+1):
                        direct=determinant([[boundary_entry(r,b,a,zeta) for a in aa] for b in range(r+1)])
                        E=sum(aa)-r*(r+1)//2
                        denominator=prod(prod(range(1,k+1)) for k in range(r+1))
                        assert direct==(-zeta)**E*F(vandermonde(aa),denominator)
                        full+=1
                    aa=tuple(range(r+1,2*r+2))
                    direct=determinant([[boundary_entry(r,b,a,zeta) for a in aa] for b in range(r+1)])
                    assert direct==(-zeta)**((r+1)**2)
                    units.append({'n':n,'r':r,'zeta':str(zeta),'exponent':(r+1)**2})
    return {'entry_factorization_checks':entries,'cauchy_minor_checks':minors,
            'full_layer_vandermonde_checks':full,'unit_coefficient_examples':units}


def run():
    slab_checks=0;antichain_checks=0;dominance_checks=0
    for c in [(F(1),F(3)),(F(2,3),F(7,4)),(F(2),F(5),F(9))]:
        for g in product(range(4),repeat=len(c)):
            n=sum(g)
            if not n:continue
            L=total(c,g)
            i=min((i for i in range(len(c)) if g[i]),key=lambda i:c[i])
            points={F(0),L,L+1}
            for k in product(*(range(x+1) for x in g)):
                x=total(c,k)
                points.update([x,x+F(1,24),x-F(1,24)])
            for d in points:
                val=alternating(c,g,d)
                assert val==slab(c,g,d,i)
                assert abs(val)<=comb(n-1,(n-1)//2)
                slab_checks+=1;antichain_checks+=1
                active=[i for i in range(len(c)) if g[i]]
                h=active[-1]
                if len(active)>1 and sum(c[i]*g[i] for i in range(h))<c[h]:
                    k=ceil(d/c[h])-1;e=d-k*c[h]
                    expected=0
                    if 0<=k<=g[h]:expected=(-1)**k*comb(g[h],k)*alternating(c[:h],g[:h],e)
                    assert val==expected
                    dominance_checks+=1
    poly_checks=0;reflection_checks=0;odd_cases=[];feasible=[]
    for c in [(1,),(1,3),(3,101),(1,9,99),(2,5),(2,4,7)]:
        for g in product(range(1,4),repeat=len(c)):
            p=polynomial(c,g);N=sum(g);L=total(c,g)
            assert len(p)==L
            for d in range(1,L+1):
                assert p[d-1]==alternating(c,g,d)
                assert p[L-d]==(-1)**(N-1)*p[d-1]
                poly_checks+=1;reflection_checks+=1
            if all(x%2 for x in c):
                assert sum((-1)**i*x for i,x in enumerate(p))==2**(N-1)
                best=max(abs(x) for x in p)
                assert best*L>=2**(N-1)
                d=next(i+1 for i,x in enumerate(p) if abs(x)==best and 2*(i+1)>=L+1)
                w=tuple((x+1)//2 for x in c);H=N+d;s=N
                assert s<H and 2*total(w,g)>=H and total(w,g)<=H
                odd_cases.append({'costs':list(c),'gamma':list(g),'N':N,'L':L,
                                  'upper_half_d':d,'max_abs_coefficient':best,
                                  'lower_bound':str(F(2**(N-1),L))})
                feasible.append({'weights':list(w),'H':H,'s':s,'gamma':list(g)})
    result={'status':'passed','slab_identity_checks':slab_checks,
            'antichain_bound_checks':antichain_checks,
            'dominance_reduction_checks':dominance_checks,
            'generating_polynomial_coefficient_checks':poly_checks,
            'reflection_checks':reflection_checks,'odd_cost_cases':odd_cases,
            'feasible_source_witnesses':feasible,'one_dimensional_boundary_minors':cauchy_checks(),
            'scope':'Exact scalar boundary sums and feasible row/type budgets. No compatible mixed-determinant estimate, rank threshold verification, or pi approximation conclusion.'}
    Path(__file__).with_name('alternating_support_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:len(v) if isinstance(v,list) else v for k,v in result.items()},indent=2))


if __name__=='__main__':run()
