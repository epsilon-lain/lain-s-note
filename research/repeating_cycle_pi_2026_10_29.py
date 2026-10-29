"""Exact checks for repeating-decimal approximations to pi.

Uses rational Machin bounds from the adjacent dyadic-rectangle verifier.
This is a finite calculation; no assertion on pi's bad approximability follows.
"""

from fractions import Fraction as F
from math import gcd

from dyadic_rectangles_2026_10_29 import pi_bounds


def repeating_block(p, q):
    assert gcd(p, q) == 1 and gcd(q, 10) == 1 and q > 1
    initial = p % q
    remainder = initial
    block = []
    while True:
        digit, remainder = divmod(10 * remainder, q)
        block.append(str(digit))
        if remainder == initial:
            return ''.join(block)


def error_interval(pi_lo, pi_hi, approximation):
    if approximation < pi_lo:
        return pi_lo - approximation, pi_hi - approximation
    assert approximation > pi_hi
    return approximation - pi_hi, approximation - pi_lo


def prefix_rational(pi_lo, pi_hi, length):
    scale = 10 ** length
    left = ((pi_lo - 3) * scale).__floor__()
    right = ((pi_hi - 3) * scale).__floor__()
    assert left == right
    result = 3 + F(left, scale - 1)
    error_lo, error_hi = error_interval(pi_lo, pi_hi, result)
    assert 0 < error_lo < error_hi < F(1, scale)
    factor = gcd(left, scale - 1)
    assert result.denominator == (scale - 1) // factor
    return result, factor, error_lo, error_hi


def main():
    lo, hi = pi_bounds()
    assert repeating_block(22, 7) == '142857'
    assert repeating_block(23, 7) == '285714'
    block_113 = repeating_block(355, 113)
    assert len(block_113) == 112 and block_113.startswith('1415929203')
    assert pow(10, 112, 113) == 1 and pow(10, 56, 113) != 1 and pow(10, 16, 113) != 1
    assert gcd(142857, 999999) == 142857
    result, factor, error_lo, error_hi = prefix_rational(lo, hi, 6)
    assert result == F(285599, 90909) and factor == 11
    for p, q in [(22, 7), (23, 7), (355, 113), (285599, 90909)]:
        low, high = error_interval(lo, hi, F(p, q))
        print(f'{p}/{q}: abs_error≈{float((low+high)/2):.12g}, q²error≈{float(q*q*(low+high)/2):.12g}')
    for n in range(1, 41):
        result, factor, low, high = prefix_rational(lo, hi, n)
        xlo, xhi = lo - 3, hi - 3
        block = (xlo * 10**n).__floor__()
        tail_lo = 10**n*xlo - block
        tail_hi = 10**n*xhi - block
        # The exact identity, verified by matching rational interval bounds.
        signed_lo = F(10**n - 1, factor**2) * (tail_lo - xhi)
        signed_hi = F(10**n - 1, factor**2) * (tail_hi - xlo)
        diff_lo = result.denominator**2 * (lo - result)
        diff_hi = result.denominator**2 * (hi - result)
        assert signed_lo <= diff_lo <= diff_hi <= signed_hi
    print('All exact checks passed: periods and repeating pi prefixes n=1..40.')


if __name__ == '__main__':
    main()
