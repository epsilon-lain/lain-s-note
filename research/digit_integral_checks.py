"""Finite exact checks for the integral discrepancy model, Python stdlib only."""
from fractions import Fraction as F
from pathlib import Path
from random import Random
import sys

sys.path.insert(0, str(Path(__file__).parent / 'pi-approximation-2026-10-27'))
from verify import pi_bounds, rounded_enclosure, certified_digits


def pi_examples():
    low, high = pi_bounds()
    lo, hi, scale = rounded_enclosure(low, high)
    for p in (8, 114):
        q = p - 1
        digits = certified_digits(lo-3*scale, hi-3*scale, scale, p, 70)
        e_low = sum((F(abs(d-digits[0]), p**j)
                     for j, d in enumerate(digits, 1)), F())
        e_high = e_low + F(1, p**70)
        print(f'base={p} first_digits={digits[:8]} '
              f'q²E≈{float(q*q*(e_low+e_high)/2):.14g} '
              f'tail_bound={float(q*q*(e_high-e_low)):.3g}')


def exact_checks():
    rng = Random(10101)
    for _ in range(5000):
        p = rng.randrange(3, 80)
        q = p - 1
        x = F(rng.randrange(1, 10000), 10000)
        n = rng.randrange(1, q)
        c = F(n, q)
        tail = x
        digits = []
        for __ in range(60):
            product = p*tail
            d = product.numerator // product.denominator
            digits.append(d)
            tail = product-d
        error = abs(x-c)
        discrepancy = sum((F(abs(d-n), p**j)
                           for j, d in enumerate(digits, 1)), F())
        remainder = F(1, p**60)
        bound = max(c, 1-c)
        coefficient = (1-bound)/(1+bound)
        assert coefficient*discrepancy <= error <= discrepancy+remainder
        first = digits[0]
        first_discrepancy = sum((F(abs(d-first), p**j)
                                 for j, d in enumerate(digits, 1)), F())
        assert first_discrepancy <= discrepancy+remainder
    print('5000 exact rational finite-prefix checks passed')


if __name__ == '__main__':
    pi_examples()
    exact_checks()
