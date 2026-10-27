"""Check adjacent-base wave bounds with exact rational arithmetic."""
from fractions import Fraction as F
from pathlib import Path
from random import Random
import sys

sys.path.insert(0, str(Path(__file__).parent / 'pi-approximation-2026-10-27'))
from verify import pi_bounds, rounded_enclosure, certified_digits


def wave_prefix(x, q, length=45):
    p = q+1
    r = x
    digits = []
    for _ in range(length):
        v = p*r
        d = v.numerator//v.denominator
        r = v-d
        digits.append(d)
    low = sum((F(abs(d-digits[0]), p**j)
               for j, d in enumerate(digits, 1)), F())
    return low, low+F(1, p**length)


def pi_wave(q):
    low, high = pi_bounds()
    lo, hi, scale = rounded_enclosure(low, high)
    p = q+1
    digits = certified_digits(lo-3*scale, hi-3*scale, scale, p, 60)
    low = sum((F(abs(d-digits[0]), p**j)
               for j, d in enumerate(digits, 1)), F())
    return q*q*low, q*q*(low+F(1, p**60))


def main():
    pi_lo, _ = pi_bounds()
    b113_lo, b113_hi = pi_wave(113)
    b114_lo, b114_hi = pi_wave(114)
    lower = 114*(pi_lo-3-b113_hi/113)
    assert lower > 16 and b114_lo > lower
    print(f'B(113)≈{float((b113_lo+b113_hi)/2):.12g}')
    print(f'B(114)≈{float((b114_lo+b114_hi)/2):.12g}')
    print(f'adjacent inequality certifies B(114)>{float(lower):.10g}')

    rng = Random(2948)
    for _ in range(5000):
        q = rng.randrange(2, 300)
        x = F(rng.randrange(1, 10000), 10000)
        _, e1 = wave_prefix(x, q)
        _, e2 = wave_prefix(x, q+1)
        assert q*e1+(q+1)*e2 >= min(x, 1-x)
    print('5000 rational finite-prefix checks passed')

    p, q = 3, 2
    for _ in range(20):
        assert p*p-2*q*q == 1
        p, q = 3*p+4*q, 2*p+3*q
    print('20 positive Pell identities checked')


if __name__ == '__main__':
    main()
