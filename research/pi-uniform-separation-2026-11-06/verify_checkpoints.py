"""Independent exact checks for the inherited recurrence and pi^2 lower bound.

Python standard library only. Does not import or write any earlier research file.
Finite checks support the formulas; the all-index proofs are in CHECKPOINTS.md
and the cited inherited proof files.
"""
from fractions import Fraction as F
from math import gcd, isqrt, factorial
from pathlib import Path
import json


def atan_interval(t, count):
    s = sum((F((-1)**k, (2*k+1)*t**(2*k+1)) for k in range(count)), F())
    next_term = F((-1)**count, (2*count+1)*t**(2*count+1))
    return min(s, s+next_term), max(s, s+next_term)


def valuation(n, p):
    assert n
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


def verify():
    al, au = atan_interval(5, 320)
    bl, bu = atan_interval(239, 100)
    pl, pu = 16*al-4*bu, 16*au-4*bl
    assert 3 < pl < pu < 4

    u, v = [1, 3], [F(0), F(5)]
    for n in range(1, 25):
        un = F((11*n*n+11*n+3)*u[n]+n*n*u[n-1], (n+1)**2)
        assert un.denominator == 1
        u.append(un.numerator)
        v.append(((11*n*n+11*n+3)*v[n]+n*n*v[n-1])/(n+1)**2)
    rows, identity_checks = [], 0
    for n in range(1, 25):
        center, neighbor = 6*v[n]/u[n], 6*v[n+1]/u[n+1]
        delta = F(30, (n+1)**2*u[n]*u[n+1])
        assert abs(center-neighbor) == delta
        assert min(center, neighbor) < pl*pl < pu*pu < max(center, neighbor)
        A, B = center.numerator, center.denominator
        assert gcd(A, B) == 1
        assert not (isqrt(A)**2 == A and isqrt(B)**2 == B)
        budget = F(1, 2*B*delta)
        qmax = isqrt(budget.numerator//budget.denominator)
        assert B*qmax*qmax*delta <= F(1, 2)
        assert B*(qmax+1)**2*delta > F(1, 2)
        for q in {1, max(1, qmax//2), qmax}:
            assert q > 0
            floorlo = (q*pl).numerator//(q*pl).denominator
            floorhi = (q*pu).numerator//(q*pu).denominator
            assert floorlo == floorhi
            for p in (floorlo, floorlo+1):
                assert 3*q <= p <= 4*q
                h = B*p*p-A*q*q
                assert h != 0
                # At either rational Machin endpoint the signed identity is exact.
                for endpoint in (pl, pu):
                    rhs = F(h, B)-q*q*(endpoint*endpoint-center)
                    assert rhs == p*p-q*q*endpoint*endpoint
                if F(p, q) < pl:
                    quality_lower = q*(q*pl-p)
                else:
                    assert F(p, q) > pu
                    quality_lower = q*(p-q*pu)
                assert quality_lower > F(1, 16*B)
                identity_checks += 1
        if n in (1, 2, 4, 8, 16, 24):
            rows.append({'n': n, 'B': str(B), 'q_max': str(qmax),
                         'uniform_on_this_finite_range': str(F(1, 16*B))})

    # The h=0 exception is tested rather than silently discarded.
    square_cases = 0
    for p, q in ((3, 1), (22, 7), (355, 113)):
        A, B = p*p, q*q
        assert gcd(A, B) == 1 and B*p*p-A*q*q == 0
        assert isqrt(A)**2 == A and isqrt(B)**2 == B
        square_cases += 1

    aa, bb = [1, 1], [0, 1]
    for n in range(1, 512):
        aa.append((2*n+1)*aa[n]+n*n*aa[n-1])
        bb.append((2*n+1)*bb[n]+n*n*bb[n-1])
    for n, a in enumerate(aa):
        k, r = divmod(n, 4)
        expected = (2*k, 2*k, 2*k+2, 2*k+3+valuation(k+1, 2))[r]
        assert valuation(a, 2) == expected
    surplus = []
    for n, p in ((72, 17), (107, 73), (166, 131)):
        gn = gcd(aa[n], 4*bb[n])
        vg, vf = valuation(gn, p), valuation(factorial(n), p)
        assert vg > vf
        surplus.append({'n': n, 'p': p, 'gcd_valuation': vg,
                        'factorial_valuation': vf})
    return {'status': 'PASS', 'source_centers': 24,
            'norm_identity_and_finite_lower_bound_checks': identity_checks,
            'rational_square_zero_norm_cases': square_cases,
            'two_adic_indices': len(aa),
            'odd_surplus_counterexamples': surplus,
            'source_height_cost_examples': rows,
            'scope': 'Finite exact checks; no all-denominator pi lower bound or contrary conclusion.'}


if __name__ == '__main__':
    result = verify()
    output = Path(__file__).with_name('checkpoint_checks.json')
    output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))
