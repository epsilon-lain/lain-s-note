"""Rational majorants of the true squared integration kernel.

No pi values enter construction of the majorant or Gram moments.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial


@lru_cache(None)
def beta(a, b):
    return F(factorial(a-1)*factorial(b-1), factorial(a+b-1))


def moments(m, degree, K=None):
    c=4*m+1
    s=2*c
    K=2*s if K is None else K
    r=F(s+K+1, 2*(K+2))
    assert r < 1
    coefficients=[F(comb(s+k-1,k),2**(s+k)) for k in range(K+1)]
    coefficients.append(F(comb(s+K,K+1),2**(s+K+1))/(1-r))
    out=[]
    for ell in range(2*degree+1):
        # (1-x)^B(1-x^2)^k=(1-x)^(B+k)(1+x)^k.
        # This beta sum is positive term by term.
        value=sum((a*sum((comb(k,j)*beta(4*m+ell+j+1,4*m+k+1)
                           for j in range(k+1)),F(0))
                   for k,a in enumerate(coefficients)),F(0))
        out.append(value)
    return out


def inner(a,b,mu):
    return sum((F(x)*y*mu[i+j] for i,x in enumerate(a) if x
                for j,y in enumerate(b) if y),F(0))


def majorant(m,x,K=None):
    c=4*m+1;s=2*c;K=2*s if K is None else K
    y=(1-x*x)/2;r=F(s+K+1,2*(K+2))
    return (sum((comb(s+k-1,k)*y**k for k in range(K+1)),F(0))
            +comb(s+K,K+1)*y**(K+1)/(1-r))/2**s


if __name__=='__main__':
    for m in range(1,4):
        for j in range(101):
            x=F(j,100)
            assert majorant(m,x)>=(1+x*x)**(-2*(4*m+1))
        print('rational majorant checks passed',m,flush=True)
