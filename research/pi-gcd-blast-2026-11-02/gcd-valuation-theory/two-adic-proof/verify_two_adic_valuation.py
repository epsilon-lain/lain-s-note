"""Exact finite audits for TWO_ADIC_NUMERATOR_VALUATION.md.

Standard library only. The all-index theorem is proved in the note; these
checks independently reconstruct its finite polynomial certificate.
"""
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import json


def valuation(x):
    x = F(x)
    if not x:
        raise ValueError('valuation of zero is not needed')
    return ((x.numerator & -x.numerator).bit_length()
            -(x.denominator & -x.denominator).bit_length())


def multiply(a,b):
    c = [F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            c[i+j] += x*y
    return c


def polynomial_certificate():
    polynomial = [F(1)]
    h = [F(0)]*10
    for j in range(1,6):
        polynomial = multiply(polynomial,[F(1),F(-1,j)])
        square = multiply(polynomial,polynomial)
        square[0] -= 1
        assert square[0] == 0
        t = F(2**j,(2*j+1)*comb(2*j,j))
        for k,c in enumerate(square[1:]):
            h[k] += t*c
    transformed = [c*2**k for k,c in enumerate(h)]
    transformed[0] -= 2
    w = [c*F(155925,4) for c in transformed]
    expected = [-133890,106091,-112150,95560,-57235,
                25853,-8200,1800,-240,16]
    assert w == list(map(F,expected))
    assert 155925%2 == 1
    assert all(valuation(c)>=2 for c in transformed if c)
    return expected


def generating_checks(order=80):
    c = [F(factorial(j)**2,factorial(2*j+1)) for j in range(order+1)]
    for k in range(1,order+1):
        assert (4*k+2)*c[k] == k*c[k-1]
    d = []
    for k in range(order+1):
        coefficient = sum((c[j]*comb(j,k-j)*4**(2*j-k)*(-1)**(k-j)
                           for j in range((k+1)//2,k+1)),F(0))
        d.append(coefficient)
        assert 2*d[k]-(d[k-1] if k else 0)-2*c[k] == 0
    return {'differential_coefficients':order,
            'duplication_coefficients':order+1}


def sum_checks(limit=32):
    count = 0
    for m in range(limit+1):
        even = sum((F(2**j*factorial(m)**2,
                      factorial(m-j)**2*factorial(2*j))
                    for j in range(m+1)),F(0))
        odd = sum((F(2**j*factorial(m)**2,
                     factorial(m-j)**2*factorial(2*j+1))
                   for j in range(m+1)),F(0))
        assert valuation(even) == m%2
        assert valuation(odd) == (0 if m%2==0 else 1+valuation(m+1))
        if m>=3 and m%2:
            r2 = F((m*(m-1))**2,6)
            r3 = F((m*(m-1)*(m-2))**2,90)
            assert r2+r3 == F((m*(m-1))**2*((m-2)**2+15),90)
            assert valuation(r2+r3)>=4
        count += 2
    return count


def main():
    w = polynomial_certificate()
    generating = generating_checks()
    sums = sum_checks()
    a = [1,1]
    limit = 2000
    for n in range(1,limit):
        a.append((2*n+1)*a[n]+n*n*a[n-1])
    for n in range(limit+1):
        m,r = divmod(n,4)
        predicted = [2*m,2*m,2*m+2,2*m+3+valuation(m+1)][r]
        assert valuation(a[n]) == predicted,(n,valuation(a[n]),predicted)
    partial = F(0)
    for j in range(101):
        t = F(2**j,(2*j+1)*comb(2*j,j))
        assert valuation(t) == j-j.bit_count()
        partial += t
        assert valuation(partial) >= (j+1)-(j+1).bit_count()
    for j in range(6,10001):
        assert j-j.bit_count()-(j.bit_length()-1)>=2
    result = {'status':'passed',
              'arithmetic':'Python standard-library exact integers and rationals',
              'integer_recurrence_valuation_checks':limit+1,
              'finite_polynomial_W_ascending':w,
              **generating,
              'even_odd_reversed_sum_checks':sums,
              'zero_series_partial_sum_residue_checks':101,
              'tail_valuation_inequality_checks':9995,
              'scope':'Finite exact audits. The all-index valuation formula follows from the proof, including the convergent-series duplication and five-term congruence; no odd-prime compensation theorem or pi bad-approximability claim.'}
    Path(__file__).with_name('two_adic_valuation_verification.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
