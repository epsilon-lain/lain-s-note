"""Certified finite pi sweep-area records and exact geometry checks.

Run with python3 needle_area_checks.py. Uses the neighboring research folder's
exact Machin bounds; no external dependencies or model calls.
"""
from fractions import Fraction as F
from pathlib import Path
from random import Random
import sys

sys.path.insert(0, str(Path(__file__).parent / 'pi-approximation-2026-10-27'))
from verify import pi_bounds, rounded_enclosure


def area(x, a, q):
    return F(q * q, 2) * abs(x - F(a, q))


def main():
    lo, hi = pi_bounds()
    lower, upper, scale = rounded_enclosure(lo, hi)
    best = None
    records = []
    for q in range(1, 100001):
        a = (2*q*lower + scale) // (2*scale)
        assert a == (2*q*upper + scale) // (2*scale)
        left, right = q*(q*lower-a*scale), q*(q*upper-a*scale)
        assert left*right > 0
        low = F(min(abs(left), abs(right)), 2*scale)
        high = F(max(abs(left), abs(right)), 2*scale)
        if best is None or high < best[0]:
            best = (low, high)
            records.append((q, a, low, high))
        elif low <= best[1] and high >= best[0]:
            raise ArithmeticError('insufficient precision to compare records')
    print('Certified minimum records for q <= 100000:')
    for q, a, low, high in records:
        print(f'q={q}, a={a}, area ~ {float((low+high)/2):.15g}')
    assert [(q, a) for q, a, _, _ in records] == [(1, 3), (7, 22), (113, 355)]

    for a, q, b, s in [(333, 106, 355, 113), (333, 106, 22, 7)]:
        assert b*q-a*s == 1 and F(a, q) < lo < hi < F(b, s)
        for x in (lo, hi):
            assert 2*(F(s, q)*area(x, a, q)+F(q, s)*area(x, b, s)) == 1
        print(f'bracket {a}/{q} < pi < {b}/{s}: exact identity passed')

    rng = Random(2848)
    for _ in range(30000):
        Q = rng.randrange(2, 2000)
        q, s = rng.randrange(Q, 2*Q), rng.randrange(Q, 2*Q)
        a, b = rng.randrange(-2*q, 2*q), rng.randrange(-2*s, 2*s)
        if a*s == b*q:
            continue
        x = F(rng.randrange(-500, 501), 113)
        assert area(x, a, q)+area(x, b, s) >= F(1, 4)
    print('30,000 seeded rational trials of the band inequality passed')


if __name__ == '__main__':
    main()
