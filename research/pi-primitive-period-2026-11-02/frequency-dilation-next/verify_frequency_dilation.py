"""Exact integer-frequency scaling and primitive-content checks.

Standard library only. Period polynomials are reconstructed from the proved
Laurent-Hankel identity, then checked against independent raw matrices.
"""
from fractions import Fraction as F
from math import comb, factorial, gcd, lcm, prod
from pathlib import Path
import json


def determinant(matrix):
    a=[list(map(F,row)) for row in matrix]
    result=F(1)
    for i in range(len(a)):
        p=next((j for j in range(i,len(a)) if a[j][i]),None)
        if p is None:return F(0)
        if p!=i:a[i],a[p]=a[p],a[i];result=-result
        pivot=a[i][i];result*=pivot
        for j in range(i+1,len(a)):
            ratio=a[j][i]/pivot
            for k in range(i+1,len(a)):a[j][k]-=ratio*a[i][k]
    return result


def trim(poly):
    while len(poly)>1 and not poly[-1]:poly.pop()
    return poly


def evaluate(poly,z):
    result=F(0)
    for c in reversed(poly):result=result*z+c
    return result


def interpolate_from_one(values):
    out=[F(0)]*len(values)
    basis=[F(1)]
    differences=list(map(F,values))
    for k in range(len(values)):
        for i,c in enumerate(basis):out[i]+=c*differences[0]
        differences=[b-a for a,b in zip(differences,differences[1:])]
        nxt=[F(0)]*(len(basis)+1)
        for i,c in enumerate(basis):
            nxt[i]-=c
            nxt[i+1]+=c/F(k+1)
        basis=nxt
    return trim(out)


def hankel_value(n,z):
    z=F(z)
    moments=[F(0)]
    for k in range(1,2*n-1):moments.append(1/z-k*moments[-1]/z)
    factor=prod(factorial(k) for k in range(n))
    return z**(n*n)*determinant([[moments[i+j] for j in range(n)] for i in range(n)])/factor**2


def polynomial(n):
    bound=n*(n-1)
    poly=interpolate_from_one([hankel_value(n,z) for z in range(1,bound+2)])
    for z in (-1,F(1,2),bound+3):assert evaluate(poly,z)==hankel_value(n,z)
    assert all(c==0 for k,c in enumerate(poly) if k<n)
    assert poly[n]==(-1)**n*(1-n)
    assert all(c==0 for k,c in enumerate(poly) if (k-n)%2)
    return poly


def raw_value(n,b,z):
    rows=[(j,s) for j in range(2) for s in range(n)]
    cols=[(h,a) for h in (0,b) for a in range(n)]
    H=max(2*n,2*n+b-2,5)
    assert F(H,n)>1 and F(4,H)<1 and F(4*n,H*H)<1
    assert all(h+2*a<=H for h,a in cols)
    assert [(s,beta) for s in range(n+1) for beta in range(2)
            if F(H,n)*s+H*beta<H]==[(s,0) for s in range(n)]
    matrix=[]
    for j,s in rows:
        row=[]
        for h,a in cols:
            row.append(sum(F(comb(a,k))* (j*z)**(a-k)*F(h**(s-k),factorial(s-k))
                           for k in range(min(s,a)+1)))
        matrix.append(row)
    return determinant(matrix)


def primes_of(b):
    primes=[];p=2
    while p*p<=b:
        if b%p==0:
            primes.append(p)
            while b%p==0:b//=p
        p+=1
    if b>1:primes.append(b)
    return primes


def valuation(a,p):
    a=abs(a);out=0
    assert a
    while a%p==0:a//=p;out+=1
    return out


def interval_mul(a,b):
    products=[x*y for x in a for y in b]
    return min(products),max(products)


def interval_eval(poly,x):
    value=(F(0),F(0))
    for c in reversed(poly):
        lo,hi=interval_mul(value,x);value=(lo+c,hi+c)
    return value


def atan_interval(q,N=80):
    total=sum(F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(N))
    nxt=F((-1)**N,(2*N+1)*q**(2*N+1))
    return min(total,total+nxt),max(total,total+nxt)


def two_pi_interval():
    a,b=atan_interval(5),atan_interval(239)
    interval=(32*a[0]-8*b[1],32*a[1]-8*b[0])
    assert interval[0]>6
    return interval


def absolute_at_period(poly,n,period):
    # Every nonzero degree has the same parity as n.
    real=[c*(-1 if ((k-n)//2)%2 else 1) if (k-n)%2==0 else F(0) for k,c in enumerate(poly)]
    lo,hi=interval_eval(real,period)
    assert lo*hi>0
    return (lo,hi) if lo>0 else (-hi,-lo)


def outward_integer(interval):
    return [str(interval[0].numerator//interval[0].denominator),
            str(-(-interval[1].numerator//interval[1].denominator))]


def main():
    period=two_pi_interval();records=[];raw_checks=content_checks=0
    for n in range(2,11):
        poly=polynomial(n)
        delta=len(poly)-1
        d=lcm(*(c.denominator for c in poly))
        integer=[int(d*c) for c in poly]
        h=gcd(*integer)
        primitive=[a//h for a in integer]
        assert gcd(*primitive)==1
        multiplier=F(d,h)
        if n<=8:
            for b in (1,2,3,5):
                for z in (F(1),F(-3,2)):
                    assert raw_value(n,b,z)==evaluate(poly,b*z)
                    raw_checks+=1
        low_sum=sum(abs(a) for a in primitive[:-1])
        lead=abs(primitive[-1])
        threshold=max(1,-(-low_sum//(3*lead)))
        assert 6*threshold>=2*F(low_sum,lead)
        assert F(1,2)*lead*6**delta>1
        cases=[]
        for b in (1,2,3,4,5,6,8,10,12):
            dilated=[a*b**k for k,a in enumerate(primitive)]
            g=gcd(*dilated)
            reduced=[a//g for a in dilated]
            assert gcd(*reduced)==1
            assert g%b**n==0 and b**delta%g==0
            assert F(multiplier,b**delta)<=multiplier/g<=F(multiplier,b**n)
            for p in primes_of(b):
                assert valuation(g,p)==min(valuation(a,p)+k*valuation(b,p)
                                          for k,a in enumerate(primitive) if a)
            assert all(p in primes_of(b) for p in primes_of(g))
            d_b=lcm(*((c*b**k).denominator for k,c in enumerate(poly)))
            if h==1:assert d_b==d//gcd(d,g)
            actual=absolute_at_period(reduced,n,period)
            cases.append({'frequency':b,'primitive_content':str(g),
                          'primitive_multiplier':str(multiplier/g),
                          'least_integer_denominator':str(d_b),
                          'primitive_actual_absolute_integer_enclosure':outward_integer(actual)})
            content_checks+=1
        records.append({'n':n,'degree':delta,'lowest_coefficient':str(poly[n]),
                        'least_integer_denominator':str(d),'integer_clearing_content':str(h),
                        'primitive_multiplier':str(multiplier),'fixed_dimension_threshold':str(threshold),
                        'polynomial':{str(k):str(c) for k,c in enumerate(poly) if c},
                        'frequency_cases':cases})
    result={'status':'passed','polynomial_reconstructions':len(records),
            'independent_raw_scaling_checks':raw_checks,'primitive_content_cases':content_checks,
            'records':records,
            'scope':'Exact scaling, primitive coefficient content, valid budgets and fixed-n frequency threshold. '
                    'No optimal joint n,b choice, infinite root-convergence, or rational approximation claim.'}
    Path(__file__).with_name('frequency_dilation_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))


if __name__=='__main__':main()
