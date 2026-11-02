"""Independent exact spot checks for the spatial-coupling certificate."""
from fractions import Fraction as F
from itertools import product
from math import comb, factorial, gcd, prod
from pathlib import Path
import json


def determinant(matrix):
    a = [list(map(F, row)) for row in matrix]
    value = F(1)
    for i in range(len(a)):
        p = next((j for j in range(i, len(a)) if a[j][i]), None)
        if p is None:
            return F(0)
        if p != i:
            a[i], a[p] = a[p], a[i]
            value = -value
        pivot = a[i][i]
        value *= pivot
        for j in range(i+1, len(a)):
            ratio = a[j][i]/pivot
            for k in range(i+1, len(a)):
                a[j][k] -= ratio*a[i][k]
    return value


def entry(row, col, z):
    j, s, beta = row
    h, alpha = col
    if any(a < b for a, b in zip(alpha, beta)):
        return F(0)
    degree = sum(alpha)-sum(beta)
    return prod(comb(a, b) for a, b in zip(alpha, beta))*sum(
        F(comb(degree, k))*F(j*z)**(degree-k)*
        F(h**(s-k), factorial(s-k)) for k in range(min(s, degree)+1))


def rows(H, v):
    beta_order = [(0,0), (1,0), (0,1)]
    expected = {(s, b) for b in product(range(H+1), repeat=2)
                for s in range(H) if s+sum(x*y for x,y in zip(v,b)) < H}
    result = [(j,s,b) for j in (0,1) for b in beta_order
              for s in range(H) if s+sum(x*y for x,y in zip(v,b)) < H]
    assert {(s,b) for j,s,b in result if j == 0} == expected
    return result


def S_imaginary(x):
    return -13*x**6+720*x**4-6480*x*x-63936


def main():
    alpha = [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(1,3),(1,4),(1,6)]
    r18 = rows(6,(4,5))
    c18 = [(h,a) for h in (0,1) for a in alpha]
    assert len(r18) == len(c18) == 18
    assert max(F(h)+F(2,5)*a[0]+F(1,2)*a[1] for h,a in c18) == F(22,5)
    r14 = rows(4,(2,3))
    c14 = [(h,a) for h in (0,1) for a in
           [(k,0) for k in range(5)]+[(k,1) for k in range(2)]]
    assert len(r14) == len(c14) == 14
    comparisons = 0
    for z in (F(1), F(-2), F(1,2), F(7,3)):
        raw = determinant([[entry(r,c,z) for c in c18] for r in r18])
        assert raw == z**21*(13*z**6+720*z**4+6480*z*z-63936)/864
        raw = determinant([[entry(r,c,z) for c in c14] for r in r14])
        assert raw == -2*z**13
        comparisons += 2
    lo,hi = F(6302063,10**6), F(6302064,10**6)
    assert S_imaginary(lo)>0>S_imaginary(hi)
    assert 39 < lo*lo < hi*hi < 40
    assert 13*39**2-480*39+2160>0 and 26*39-480>0
    clo,chi = F(62954,10**4),F(62956,10**4)
    assert clo**4-36*clo**2-144<0<chi**4-36*chi**2-144
    assert F(44,7) < clo < chi < lo
    assert all((13*u**3+720*u*u+6480*u-63936)%11 for u in range(11))
    assert 63936%13 != 0
    denominator_checks = 0
    for b in range(2,81):
        for a in range(6*b,6*b+11):
            if gcd(a,b) != 1:
                continue
            x=F(a,b)
            residual=S_imaginary(x)
            assert residual and residual.denominator == b**6//gcd(b,13)
            full=x**21*residual  # i^21/i = 1.
            assert full.denominator == b**27//gcd(b,13)
            denominator_checks += 2
    result={'status':'passed','independent_raw_determinant_checks':comparisons,
            'exact_rational_value_denominator_checks':denominator_checks,
            'root_bracket_and_comparator':'passed',
            'cubic_irreducibility_and_nonsquare_norm_checks':'passed',
            'scope':'Exact supplementary checks; the main full polynomial identity certificate and general profile proof are separate.'}
    Path(__file__).with_name('coupled_independent_review_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
