"""Exact arithmetic checks for normalized reflection factors.

Python standard library only. General identities and comparison estimates
are proved in PRIMITIVE_REFLECTION_FACTORS.md; no finite extrapolation.
"""
from fractions import Fraction as F
from functools import cache
from math import comb, factorial, gcd, prod
from pathlib import Path
import json


@cache
def A(n):
    return prod(factorial(j) for j in range(n))


def det(matrix):
    a = [list(map(F, row)) for row in matrix]
    sign = 1
    previous = F(1)
    for k in range(len(a)-1):
        p = next((i for i in range(k,len(a)) if a[i][k]),None)
        if p is None:
            return F(0)
        if p != k:
            a[k],a[p] = a[p],a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):
                a[i][j] = (a[i][j]*pivot-a[i][k]*a[k][j])/previous
            a[i][k] = F(0)
        previous = pivot
    return sign*a[-1][-1]


def evaluate(poly,z):
    value = F(0)
    for c in reversed(poly):
        value = value*z+c
    return value


def interpolate(values):
    out = [F(0)]*len(values)
    basis = [F(1)]
    differences = list(map(F,values))
    for k in range(len(values)):
        for j,c in enumerate(basis):
            out[j] += c*differences[0]
        differences = [b-a for a,b in zip(differences,differences[1:])]
        nxt = [F(0)]*(len(basis)+1)
        for j,c in enumerate(basis):
            nxt[j] -= c*F(k,k+1)
            nxt[j+1] += c/F(k+1)
        basis = nxt
    while len(out)>1 and out[-1] == 0:
        out.pop()
    return out


def reflection(n,z,sign):
    return det([[sum(F(comb(a,k))*F(z,2)**(a-k)
                     *(F(1,2)**(s-k)+sign*(-1)**a*F(-1,2)**(s-k))
                     /factorial(s-k) for k in range(min(a,s)+1))
                 for a in range(n)] for s in range(n)])


def raw(n,z):
    return det([[sum(F(comb(a,k))*F(j*z)**(a-k)*F(h**(s-k),factorial(s-k))
                     for k in range(min(a,s)+1))
                 for h in (0,1) for a in range(n)]
                for j in (0,1) for s in range(n)])


def multiply(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            c[i+j]+=x*y
    return c


def at_imaginary(poly,x):
    assert all(c==0 for k,c in enumerate(poly) if k%2)
    return sum(c*(-1)**(k//2)*x**k for k,c in enumerate(poly) if c)


def sufficient_height(b):
    n=24*b
    while True:
        r=n//2
        B=A(n)//A(r)**2
        Delta=r*(r-1)
        ratio=4*b**(2*Delta)*F(44,7)**(n*n-n)*F(B**2*A(n)**2,A(2*n))
        if ratio<=1:
            break
        n+=2
    return {'b':b,'even_n':n,'factor_degree':Delta,
            'squared_sufficient_condition_certified':True,
            'ratio_numerator_bits':ratio.numerator.bit_length(),
            'ratio_denominator_bits':ratio.denominator.bit_length()}


def main():
    records=[]
    denominator_checks=0
    raw_checks=0
    for n in (4,6,8,10,12):
        r=n//2
        B=A(n)//A(r)**2
        Delta=r*(r-1)
        S=n*(n-1)//2
        factors=[]
        for sign in (1,-1):
            poly=interpolate([reflection(n,z,sign) for z in range(r*r+1)])
            assert len(poly)-1==r*r
            assert next(k for k,c in enumerate(poly) if c)==r
            leading_sign=(-1)**(r*(r-1)//2 if sign==1 else r*(r+1)//2)
            assert poly[-1]==F(leading_sign,B)
            assert all(c==0 for k,c in enumerate(poly) if (k-r)%2)
            for z in (F(-3,2),F(2,3),F(r*r+2)):
                assert evaluate(poly,z)==reflection(n,z,sign)
            monic=[leading_sign*B*c for c in poly[r:]]
            assert len(monic)-1==Delta and monic[-1]==1
            assert all(c.denominator==1 for c in monic)
            assert gcd(*(int(c) for c in monic))==1
            M=B*factorial(n)*2**n*n*n*5**S
            coefficient_majorant=sum(k*abs(c)*8**(k-1) for k,c in enumerate(monic) if k)
            assert coefficient_majorant<=M
            for b in range(2,16):
                for a in range(6*b,7*b+1):
                    if gcd(a,b)!=1:
                        continue
                    value=at_imaginary(monic,F(a,b))
                    assert value!=0 and value.denominator==b**Delta
                    assert abs(value)>=F(1,b**Delta)
                    denominator_checks+=1
            factors.append(monic)
        product=[F(0)]*n+[(-1)**r*c for c in multiply(*factors)]
        for z in (F(-1),F(1,2),F(3)):
            assert evaluate(product,z)==B**2*raw(n,z)
            raw_checks+=1
        records.append({'n':n,'reflection_degree':r*r,'reflection_order':r,
                        'factor_degree':Delta,'normalization_B':str(B),
                        'factors':[{'sign':sign,'coefficients':{str(k):str(c) for k,c in enumerate(p) if c}}
                                   for sign,p in zip(('+','-'),factors)],
                        'integer_monic_even':True,'derivative_majorant_check':True})
    result={'status':'passed','even_dimensions':len(records),
            'independent_raw_product_checks':raw_checks,
            'exact_rational_value_denominator_checks':denominator_checks,
            'records':records,
            'finite_sufficient_heights':[sufficient_height(b) for b in (2,3,7,10)],
            'scope':'Finite exact checks supplement the general proof. No claim of general irreducibility, '
                    'subquadratic factor degrees, q^-2 lower bounds, or pi bad approximability.'}
    Path(__file__).with_name('primitive_reflection_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))


if __name__=='__main__':
    main()
