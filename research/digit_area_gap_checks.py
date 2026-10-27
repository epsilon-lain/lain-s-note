"""Exact finite checks for the nonnegative-area gap lemma."""
from fractions import Fraction as F
from pathlib import Path
from random import Random
import sys

sys.path.insert(0, str(Path(__file__).parent / 'pi-approximation-2026-10-27'))
from verify import pi_bounds


def partial_area(x, q, terms=30):
    p = q+1
    r = x
    digits = []
    for _ in range(terms):
        product = p*r
        d = product.numerator//product.denominator
        digits.append(d)
        r = product-d
    area = sum((F(abs(d-digits[0]), p**j)
                for j, d in enumerate(digits, 1)), F())
    return digits[0], q*q*(area+F(1, p**terms))


def main():
    lo, hi = pi_bounds()
    assert lo < hi < F(355, 113)
    area_upper = F(113**2, 2)*(F(355, 113)-lo)
    ratio = F(33173, 113)
    assert 2*area_upper*(ratio+1/ratio) < 1
    assert 355*7-22*113 == -1
    print('Certified: no new pi area record for 114 <= q <= 33173')
    print(f'355/113 area < {float(area_upper):.15g}')

    rng = Random(284800)
    checked = 0
    for _ in range(5000):
        q = rng.randrange(2, 250)
        s = rng.randrange(q+1, 2*q+1)
        x = F(rng.randrange(1, 10000), 10000)
        n, bound_q = partial_area(x, q)
        m, bound_s = partial_area(x, s)
        if n*s == m*q:
            continue
        assert 1 <= F(s, q)*bound_q+F(q, s)*bound_s
        checked += 1
    print(f'Area gap: {checked} exact rational finite-prefix comparisons passed')

    for k in range(2, 7):
        # q_k^2 times the standard upper bound for the Liouville tail, divided by 2.
        from math import factorial
        q = 10**factorial(k)
        tail = F(2, 10**factorial(k+1))
        assert F(q*q, 2)*tail == F(1, 10**((k-1)*factorial(k)))
    print('Liouville tail bounds checked for k=2,...,6')


if __name__ == '__main__':
    main()
