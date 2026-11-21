"""Exact test of the affine construction from selected 1/pi seeds.

Run next to verify_primitive_window_extrema.py. All comparisons use integers
or rational intervals derived from Machin's formula.
"""

from fractions import Fraction
from math import gcd, isqrt

from verify_primitive_window_extrema import convergents, reciprocal_pi


def verify_public_gap():
    # The giant digit and the length of the reported prefix are external inputs.
    giant_index = 383_059_989_913  # pi, starting at a_0
    term_count = 653_520_000_000
    steps = term_count - 1 - giant_index
    blocks = (steps - 1) // 10
    assert 1618**10 > 100 * 1000**10
    assert blocks > 0
    print("Conditional denominator gap: factor >", f"10^{2 * blocks}")

    max_order = 10**17
    half_log_q_lower = (giant_index - 1) // 4
    gcd_bound_per_order = 2 * (isqrt(2 * max_order) + 1) * (2 * max_order).bit_length()
    assert gcd_bound_per_order == 51_876_777_136 < half_log_q_lower
    assert 12 * (845_188_380_288 + 2)**2 < 2**100
    print("Padé denominator inequality certified through order", max_order)


def scan(points, index, following_digit, low, high, scale):
    a, r = points[index]
    ap, rp = points[index - 1]
    assert abs(a * rp - ap * r) == 1
    b, s = points[index + 1]
    assert (b, s) == (following_digit * a + ap, following_digit * r + rp)
    assert gcd(s + 1, b * r + a) == 1

    cutoff = 4 * (following_digit + 2) * r**2
    count = 0
    best_gain = Fraction(0)
    best_index = best_gcd = 0
    for j in range(index + 1, len(points)):
        p, q = points[j]
        if q < cutoff:
            continue
        count += 1
        g = gcd(p * r + a, q + 1)
        direction = a * q - p * r
        assert g == gcd(direction, q + 1)
        assert gcd(g, r) == 1
        n, d = (p * r + a) // g, r * (q + 1) // g
        assert gcd(n, d) == 1
        # |normalized linear form| > 3/4 of the seed's exact error.
        seed_lower = min(abs(r * low - a * scale), abs(r * high - a * scale))
        fit_lower = min(abs(d * low - n * scale), abs(d * high - n * scale))
        assert (r * low - a * scale) * (r * high - a * scale) > 0
        assert (d * low - n * scale) * (d * high - n * scale) > 0
        assert 4 * d * fit_lower * g**2 > 3 * (q + 1) * r * seed_lower
        gain = Fraction(g**2, q + 1)
        if gain > best_gain:
            best_gain, best_index, best_gcd = gain, j, g
    assert best_gain > 0
    assert best_gain < Fraction(3, 4)
    ceiling = next(k for k in range(10000) if best_gain >= Fraction(1, 10**k))
    assert best_gain < Fraction(1, 10 ** (ceiling - 1))
    print(f"A={following_digit}: {count} candidates, best j={best_index}, gcd={best_gcd}")
    print(f"    10^-{ceiling} <= gain < 10^-{ceiling - 1}")


if __name__ == "__main__":
    scale = 1 << 10300
    low, high = reciprocal_pi(scale)
    points = convergents(low, high, scale, 1 << 4993)
    # Read input convergents only through index 2000, with a_0(1/pi)=0.
    points = points[:2001]
    for index, digit in ((4, 292), (431, 20776)):
        scan(points, index, digit, low, high, scale)
    verify_public_gap()
