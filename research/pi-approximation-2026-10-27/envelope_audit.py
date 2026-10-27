"""Exact checks and finite certificates for the root-envelope obstruction.

Python standard library only. The infinite theorem is proved in
ROOT_ENVELOPE_OBSTRUCTION.md, not by this finite computation.
"""
from fractions import Fraction as F
from math import gcd, isqrt
from pathlib import Path
import json

from verify import pi_bounds, rounded_enclosure

ROOT = Path(__file__).parent


def sqrt_interval(a, b, places=80):
    """Enclose sqrt(b/a) using integer arithmetic."""
    scale = 10 ** places
    k = isqrt(b * scale * scale // a)
    lo, hi = F(k, scale), F(k + 1, scale)
    assert lo * lo <= F(b, a) < hi * hi
    return lo, hi


def ceil_sqrt(r):
    k = isqrt(r.numerator // r.denominator)
    return k if k * k == r else k + 1


def finite_height_certificate(cutoff, pi_lo, pi_hi):
    """Certify q^2 G_infinity(q) <= 1/[7(K+1)] for all q>=q0.

    Enumerates a superset of the m-family: all primitive irrational roots
    sqrt(B/A) in (3,4), A<=K, on either side of pi. Higher A are bounded
    analytically, not sampled. All separation and q0 decisions are exact.
    """
    worst_ratio, worst_pair = F(0), None
    count = 0
    for a in range(1, cutoff + 1):
        for b in range(9 * a + 1, 16 * a):
            if gcd(a, b) != 1 or isqrt(a * b) ** 2 == a * b:
                continue
            lo, hi = sqrt_interval(a, b)
            assert 3 < lo <= hi < 4
            separation = max(pi_lo - hi, lo - pi_hi)
            if separation <= 0:
                raise ArithmeticError('Increase enclosure precision')
            c_upper = 1 / (a * (4 + lo))
            ratio = c_upper / separation
            if ratio > worst_ratio:
                worst_ratio, worst_pair = ratio, [a, b]
            count += 1
    q0 = ceil_sqrt(worst_ratio)
    assert q0 * q0 >= worst_ratio
    if q0:
        assert (q0 - 1) ** 2 < worst_ratio
    return {'height_cutoff': cutoff, 'enumerated_primitive_roots': count,
            'q_from_which_small_heights_contribute_zero': q0,
            'uniform_upper_bound_on_q_squared_infinite_envelope':
                str(F(1, 7 * (cutoff + 1))),
            'worst_root_A_B': worst_pair,
            'max_c_upper_over_eta_lower_exact': str(worst_ratio)}


def main():
    # Independently regenerate the supplied pi enclosure.
    lo, hi = pi_bounds()
    l, u, scale = rounded_enclosure(lo, hi)
    supplied = json.loads((ROOT / 'pi_interval.json').read_text())
    assert (l, u, scale) == tuple(map(int,
        (supplied['lower'], supplied['upper'], supplied['scale'])))
    pi_lo, pi_hi = F(l, scale), F(u, scale)

    # Verify the signed residual identity and simple worked certificate.
    a, b = 14891952, 146977675
    assert gcd(a, b) == 1
    rl, ru = sqrt_interval(a, b)
    eta_hi = pi_hi - rl
    c_lo = 1 / (a * (4 + ru))
    assert 0 < pi_lo - ru <= eta_hi < F(13, 10**14)
    assert c_lo > F(9, 10**9)
    # The following rational surrogate identity is exact, for either tail.
    checks = 0
    for p in (F(22, 7), F(355, 113), pi_lo, pi_hi):
        for r in (rl, ru):
            delta = a * (p * p - r * r)
            eta = p - r
            c = 1 / (a * (4 + r))
            assert c / eta == (p + r) / ((4 + r) * delta)
            checks += 1
    # Exercise the square-adjustment proof on synthetic square floors.
    for s in range(2, 1001):
        assert (s - 1) ** 2 < s * s - 1 < s * s
        assert isqrt(s * s - 1) ** 2 != s * s - 1

    certificates = [finite_height_certificate(k, pi_lo, pi_hi)
                    for k in (1, 5, 10, 25, 50, 100, 200)]
    out = {'pi_enclosure_regenerated': True,
           'residual_identity_checks': checks,
           'synthetic_square_adjustment_checks': 999,
           'worked_example_verified': True,
           'certificates': certificates,
           'interpretation': 'Upper bounds on this certificate envelope, '
               'not upper bounds on the actual rational approximation error.'}
    (ROOT / 'envelope_audit_results.json').write_text(
        json.dumps(out, indent=2) + '\n', encoding='utf-8')
    for row in certificates:
        print('K=', row['height_cutoff'], 'q>=',
              row['q_from_which_small_heights_contribute_zero'],
              'q^2 G_infinity <=',
              row['uniform_upper_bound_on_q_squared_infinite_envelope'])


if __name__ == '__main__':
    main()
