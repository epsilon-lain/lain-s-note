"""Exact finite checks for the reciprocal flow and observation-clock note.

Run beside verify_transport.py. No floating point is used for assertions.
General theorems are proved in OBSERVATION_CLOCK_RECIPROCAL_FLOW.md.
"""

from fractions import Fraction as F
from math import factorial, gcd
from pathlib import Path
import json

from verify_transport import certify_digits, interval_record


def flow(z, elapsed):
    assert 1 - elapsed * z > 0
    return z / (1 - elapsed * z)


def value(n, elapsed):
    assert elapsed < n
    return 1 / (F(n) - elapsed)


def positive_absolute_interval(low, high):
    assert low <= high and (low > 0 or high < 0)
    return (low, high) if low > 0 else (-high, -low)


def main():
    counts = {'flow_composition': 0, 'affine_ratio': 0,
              'sample_reduction': 0, 'target_divisibility': 0,
              'factorial_coverage': 0, 'normalized_identity': 0}
    for n in range(1, 9):
        z = F(1, n)
        for a in range(9):
            for b in range(9):
                first, second = F(a * n, 20), F(b * n, 20)
                assert flow(flow(z, first), second) == flow(z, first + second)
                counts['flow_composition'] += 1
    for n0 in range(-3, 4):
        for p0 in range(1, 5):
            for a in range(-2, 3):
                for b in range(1, 5):
                    for s in (F(0), F(1, 3), F(7, 2)):
                        N, P = n0 + a*s, p0 + b*s
                        assert N/P == F(a, b) + F(b*n0-a*p0, b*P)
                        assert (N+a)/(P+b)-N/P == (a*P-b*N)/(P*(P+b))
                        counts['affine_ratio'] += 1
    for n in (1, 4):
        for M in range(1, 49):
            for k in range(n*M):
                g = gcd(M, k)
                out = value(n, F(k, M))
                assert (out.numerator, out.denominator) == (M//g, (n*M-k)//g)
                counts['sample_reduction'] += 1
        for p in range(1, 41):
            for q in range(1, 41):
                if gcd(p, q) != 1 or n*p < q:
                    continue
                s = F(n*p-q, p)
                assert s.denominator == p and value(n, s) == F(p, q)
                for M in range(1, 41):
                    assert ((M*s).denominator == 1) == (M % p == 0)
                    counts['target_divisibility'] += 1
                assert (factorial(p)*s).denominator == 1
                counts['factorial_coverage'] += 1
    assert value(4, F(2, 3)) == F(3, 10)
    for B in range(1, 61):
        for A in range(B):
            if gcd(A, B) != 1:
                continue
            s = F(A, B)
            Q = B-A
            x = value(1, s)
            assert x.denominator == Q
            for alpha in (F(7, 3), F(22, 7), F(355, 113)):
                target_time = 1-1/alpha
                delta = s-target_time
                assert Q**2*abs(x-alpha) == (1-alpha*delta)*B**2*abs(delta)
                counts['normalized_identity'] += 1
    low, high, _ = certify_digits(16)
    sigma_low, sigma_high = 1-1/low, 1-1/high
    examples = []
    for p, q in ((22, 7), (333, 106), (355, 113), (103993, 33102)):
        r, s = F(p, q), F(p-q, p)
        assert value(1, s) == r and s.denominator == p
        spatial = positive_absolute_interval(r-high, r-low)
        temporal = positive_absolute_interval(s-sigma_high, s-sigma_low)
        factor_low, factor_high = low*q/p, high*q/p
        # Enclose both sides of the identity independently; interval products
        # may overestimate because the factors are correlated through pi.
        right_low = factor_low*p*p*temporal[0]
        right_high = factor_high*p*p*temporal[1]
        assert right_low <= q*q*spatial[1]
        assert right_high >= q*q*spatial[0]
        examples.append({'rational_value': str(r), 'event_time': str(s),
                         'minimum_uniform_clock_M': p,
                         'spatial_normalized_error': interval_record(q*q*spatial[0], q*q*spatial[1]),
                         'temporal_normalized_error': interval_record(p*p*temporal[0], p*p*temporal[1]),
                         'exact_conversion_factor': interval_record(factor_low, factor_high)})
    result = {'status': 'exact finite checks passed; P2848 unresolved',
              'checks': counts,
              'unit_interval_counterexample': {'n': 4, 'time': '2/3', 'value': '3/10'},
              'examples': examples,
              'limits': 'No finite check establishes an infinite liminf; no Lean run.'}
    destination = Path(__file__).with_name('observation_clock_results.json')
    destination.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
