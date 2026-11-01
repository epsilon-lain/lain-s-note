"""Exact checks of the m=1 infinite-H construction.

Proof: support_rank_constraints.md, Section 1. Rational r values below test
the polynomial determinant identity, not an approximation to pi.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json


def multiply(a,b,N):
    c=[F(0)]*(N+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b[:N+1-i]):
            c[i+j]+=x*y
    return c


def power(a,k,N):
    out=[F(1)]+[F(0)]*N
    for _ in range(k): out=multiply(out,a,N)
    return out


def determinant(a):
    a=[r[:] for r in a]
    answer=F(1)
    for i in range(len(a)):
        pivot=next((j for j in range(i,len(a)) if a[j][i]),None)
        if pivot is None: return F(0)
        if pivot!=i:
            a[i],a[pivot]=a[pivot],a[i]
            answer=-answer
        v=a[i][i]
        answer*=v
        for j in range(i+1,len(a)):
            ratio=a[j][i]/v
            for k in range(i+1,len(a)): a[j][k]-=ratio*a[i][k]
    return answer


def check(H,K,theta,w0,v0,w1,r,curve):
    rows=[]; columns=[]; bs=[]
    s=0
    while v0*s<H:
        q=theta*(H-v0*s)/w1
        B=(q.numerator+q.denominator-1)//q.denominator-1
        bs.append(B)
        rows.extend((j,s,beta) for j in range(K) for beta in range(B+1))
        columns.extend((s,k) for k in range(K*(B+1)))
        s+=1
    assert len(rows)==len(columns)
    assert all(w0*h+w1*k<=H for h,k in columns)
    N=max(s for _,s,_ in rows)
    if curve=='log':
        G=[F(0)]+[F((-1)**(k+1),k) for k in range(1,N+1)]
    else:
        G=[F(0)]*(N+1)
        if N>=1:G[1]=F(3,2)
        if N>=2:G[2]=F(-2,3)
    cache={}
    matrix=[]
    for j,s,beta in rows:
        line=[]
        for h,k in columns:
            if k<beta:
                line.append(F(0));continue
            key=(j,h,k-beta)
            if key not in cache:
                shifted=G[:];shifted[0]=j*r
                xpow=power(shifted,k-beta,N)
                ypow=[F(comb(h,d)) if d<=h else F(0) for d in range(N+1)]
                cache[key]=multiply(ypow,xpow,N)
            line.append(comb(k,beta)*cache[key][s])
        matrix.append(line)
    direct=determinant(matrix)
    predicted=F(1)
    for B in bs:
        for i in range(K):
            for j in range(i+1,K): predicted*=((j-i)*r)**((B+1)**2)
    assert direct==predicted and direct!=0
    return {'H':H,'K':K,'theta':str(theta),'w0':str(w0),'v0':str(v0),
            'w1':str(w1),'r_formal_placeholder':str(r),'curve':curve,
            'matrix_size':len(rows),'B_s':bs,
            'determinant_matches_exact_product':True,
            'determinant_numerator_digits':len(str(abs(direct.numerator))),
            'determinant_denominator_digits':len(str(direct.denominator))}


def main():
    cases=[]
    for H,K,theta in ((8,2,F(1,3)),(12,2,F(1,3)),(16,3,F(1,4))):
        for curve in ('log','arbitrary_zero_constant'):
            cases.append(check(H,K,theta,F(1),F(2),F(2),F(3,2),curve))
    for curve in ('log','arbitrary_zero_constant'):
        cases.append(check(8,2,F(1,4),F(1,4),F(1),F(1),F(3,2),curve))
    report={'status':'passed','arithmetic':'exact Fraction determinant elimination',
            'cases':cases,
            'scope':'Checks raw monomial matrices and exact product formulas. The infinite-H claim is proved in the accompanying note; no inference about pi badly approximable.'}
    Path(__file__).with_name('explicit_period_minor_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','cases':len(cases),'sizes':[r['matrix_size'] for r in cases]}))


if __name__=='__main__':main()
