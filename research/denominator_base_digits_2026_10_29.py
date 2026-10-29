"""Certified variable-base digit experiments for rational approximations to pi.

All pi bounds come from the adjacent dyadic rectangle script and exact Machin
alternating-series bounds. This is a finite check, not a proof of P2848.
"""

from fractions import Fraction as F
from dyadic_rectangles_2026_10_29 import pi_bounds


def certified_fractional_digits(lo, hi, base, count):
    assert base >= 2 and 3 < lo < hi < 4
    left, right = lo - 3, hi - 3
    digits = []
    for _ in range(count):
        left *= base
        right *= base
        digit = left.numerator // left.denominator
        assert digit == right.numerator // right.denominator
        digits.append(digit)
        left -= digit
        right -= digit
    return digits


def inspect(lo, hi, q):
    p_lo = (q * lo + F(1, 2)).numerator // (q * lo + F(1, 2)).denominator
    p_hi = (q * hi + F(1, 2)).numerator // (q * hi + F(1, 2)).denominator
    assert p_lo == p_hi
    p = p_lo
    if F(p, q) < lo:
        error_lo, error_hi = q * q * (lo - F(p, q)), q * q * (hi - F(p, q))
    else:
        assert F(p, q) > hi
        error_lo, error_hi = q * q * (F(p, q) - hi), q * q * (F(p, q) - lo)
    digits = certified_fractional_digits(lo, hi, q, 4)
    second, third = digits[1:3]
    if second == 0:
        index = F(third, q)
    elif second == q - 1:
        index = F(q - 1 - third, q)
    else:
        index = F(1)
    assert index <= min(F(1), error_lo) <= min(F(1), error_hi) < index + F(1, q)
    return p, digits, index, error_lo, error_hi


def main():
    lo, hi = pi_bounds()
    samples = {}
    best_hi = None
    records = []
    for q in range(2, 10001):
        p, digits, index, error_lo, error_hi = inspect(lo, hi, q)
        if best_hi is None or error_hi < best_hi:
            best_hi = error_hi
            records.append((q, p, digits, float((error_lo + error_hi) / 2)))
        if q in (6, 7, 106, 113):
            samples[q] = (p, digits, index, error_lo, error_hi)

    assert samples[7][0] == 22 and samples[7][1] == [0, 6, 6, 3]
    assert samples[113][0] == 355 and samples[113][1] == [15, 112, 112, 69]
    assert [row[0] for row in records] == [2, 7, 113]
    for q, (p, digits, index, low, high) in samples.items():
        print(f'q={q}, p={p}, first_digits={digits}, index={index}, score≈{float((low+high)/2):.12g}')
    print('Certified records for q=2..10000:', records)
    print('All clipped-score bounds verified for q=2..10000.')


if __name__ == '__main__':
    main()
