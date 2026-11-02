"""Small exact multivariate determinant checks; standard library only."""
from fractions import Fraction as F
from itertools import permutations, product
from functools import lru_cache
from math import comb, factorial
from pathlib import Path
import json


def add(a,b):
    out=dict(a)
    for k,v in b.items():
        out[k]=out.get(k,F(0))+v
        if not out[k]: del out[k]
    return out


def mul(a,b):
    out={}
    for ka,va in a.items():
        for kb,vb in b.items():
            k=tuple(x+y for x,y in zip(ka,kb))
            out[k]=out.get(k,F(0))+va*vb
    return {k:v for k,v in out.items() if v}


def determinant(a,nvars):
    out={}
    for p in permutations(range(len(a))):
        inv=sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
        term={(0,)*nvars:F((-1)**inv)}
        for i,j in enumerate(p):
            term=mul(term,a[i][j])
            if not term:break
        out=add(out,term)
    return out


@lru_cache(None)
def stirling(n,k):
    if n==k==0:return 1
    if not n or not k or k>n:return 0
    return k*stirling(n-1,k)+stirling(n-1,k-1)


def time_coefficient(h,n,shifted):
    return F(factorial(h)*stirling(n,h) if shifted else h**n,factorial(n))


def center_entry(row,col,shifted):
    j,s,beta=row
    h,alpha=col
    if any(b>a for a,b in zip(alpha,beta)):return {}
    out={}
    for d in product(*(range(a-b+1) for a,b in zip(alpha,beta))):
        if sum(d)>s:continue
        k=tuple(a-b-dd for a,b,dd in zip(alpha,beta,d))
        val=time_coefficient(h,s-sum(d),shifted)*j**sum(k)
        for a,b,dd in zip(alpha,beta,d):val*=comb(a,b)*comb(a-b,dd)
        if val:out[k]=out.get(k,F(0))+val
    return out


def slope_entry(row,col,shifted):
    # Variables are (W,zeta_1,...,zeta_m); W is a formal common center.
    j,s,beta=row
    h,alpha=col
    if any(b>a for a,b in zip(alpha,beta)):return {}
    out={}
    for d in product(*(range(a-b+1) for a,b in zip(alpha,beta))):
        if sum(d)>s:continue
        wdegree=sum(a-b-dd for a,b,dd in zip(alpha,beta,d))
        val=time_coefficient(h,s-sum(d),shifted)*j**wdegree
        for a,b,dd in zip(alpha,beta,d):val*=comb(a,b)*comb(a-b,dd)
        for gamma in product(*(range(dd+1) for dd in d)):
            coefficient=val*(-1)**sum(gamma)
            for dd,g in zip(d,gamma):coefficient*=comb(dd,g)
            key=(wdegree,)+gamma
            if coefficient:out=add(out,{key:coefficient})
    return out


def run_case(rows,cols,shifted):
    m=len(rows[0][2])
    d=tuple(sum(a[i] for h,a in cols)-sum(b[i] for j,s,b in rows) for i in range(m))
    p=determinant([[center_entry(r,c,shifted) for c in cols] for r in rows],m)
    actual=determinant([[slope_entry(r,c,shifted) for c in cols] for r in rows],m+1)
    predicted={}
    delta_bound=sum(s for j,s,b in rows)-(sum(h for h,a in cols) if shifted else 0)
    q={}
    for k,val in p.items():
        delta=tuple(dd-kk for dd,kk in zip(d,k))
        assert all(x>=0 for x in delta) and sum(delta)<=delta_bound
        q[delta]=val
        for gamma in product(*(range(dd+1) for dd in delta)):
            coefficient=val*(-1)**sum(gamma)
            for dd,g in zip(delta,gamma):coefficient*=comb(dd,g)
            predicted=add(predicted,{(sum(k),)+gamma:coefficient})
    assert predicted==actual
    radial={}
    for (wdegree,*gamma),val in actual.items():
        radial=add(radial,{(wdegree,sum(gamma)):val})
    grouped={}
    for delta,val in q.items():grouped[sum(delta)]=grouped.get(sum(delta),F(0))+val
    grouped={k:v for k,v in grouped.items() if v}
    if grouped:
        top=max(grouped)
        assert {degree for wdegree,degree in radial}==set(range(top+1))
    else:assert not radial
    # Every general zeta coefficient is exactly the sum of its homogeneous
    # slice moments, retained separately by distinct formal powers of W.
    for gamma in product(*(range(max(dd,0)+1) for dd in d)):
        moments={}
        for delta,val in q.items():
            if any(g>dd for g,dd in zip(gamma,delta)):continue
            v=val
            for dd,g in zip(delta,gamma):v*=comb(dd,g)
            moments[sum(delta)]=moments.get(sum(delta),F(0))+v
        coeff={k[0]:v for k,v in actual.items() if k[1:]==gamma}
        expected={sum(d)-k:(-1)**sum(gamma)*v for k,v in moments.items() if v}
        assert coeff==expected
    return dict(m=m,rows=len(rows),basis='shifted' if shifted else 'raw',
                center_polynomial_terms=len(p),aggregate_terms=len(actual),
                time_bound=delta_bound,radial_degree=max(grouped) if grouped else None)


def main():
    examples=[]
    rows1=[(j,s,(0,)) for j in range(2) for s in range(3)]
    cols1=[(h,(a,)) for h,a in [(0,0),(1,0),(2,0),(0,1),(1,1),(2,2)]]
    rows2=[(j,s,(0,0)) for j in range(2) for s in range(2)]
    cols2=[(0,(0,0)),(0,(1,0)),(1,(1,0)),(0,(0,1))]
    candidates=[cols2,
      [(0,(0,0)),(1,(0,0)),(0,(1,0)),(1,(1,0))],
      [(0,(0,0)),(0,(1,0)),(0,(0,1)),(0,(1,1))],
      [(0,(0,0)),(0,(1,0)),(0,(2,0)),(1,(0,1))],
      [(0,(1,0)),(1,(1,0)),(0,(0,1)),(1,(0,1))]]
    for shifted in (False,True):
        examples.append(run_case(rows1,cols1,shifted))
        for cols in candidates:examples.append(run_case(rows2,cols,shifted))
    # A complete packet with genuine transverse derivative rows.
    rows3=[(0,s,(0,0)) for s in range(3)]+[(0,0,(1,0)),(0,0,(0,1))]
    cols3=[(0,(0,0)),(1,(0,0)),(2,(0,0)),(0,(1,0)),(0,(0,1))]
    for shifted in (False,True):examples.append(run_case(rows3,cols3,shifted))
    result=dict(status='passed',arithmetic='exact rational multivariate polynomials',
                cases=examples,number_of_cases=len(examples),
                scope='Exact formal identities, support and radial binomial transform. Transcendence supplies the written nonvanishing theorem; no numerical pi claim.')
    Path(__file__).with_name('reciprocal_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
