#!/usr/bin/env python3
"""Exact arithmetic checks for the Machin follow-up."""

from decimal import Decimal, localcontext
from fractions import Fraction
from functools import cache
from math import comb, gcd, lcm


def arctangent_interval(m, terms=150):
    total = sum(Fraction((-1)**j, (2 * j + 1) * m**(2 * j + 1)) for j in range(terms))
    endpoint = total + Fraction((-1)**terms, (2 * terms + 1) * m**(2 * terms + 1))
    return min(total, endpoint), max(total, endpoint)


def pi_interval():
    lo5, hi5 = arctangent_interval(5)
    lo239, hi239 = arctangent_interval(239)
    return 16 * lo5 - 4 * hi239, 16 * hi5 - 4 * lo239


def error_bounds(height, numerator, interval):
    lo, hi = (height * endpoint - numerator for endpoint in interval)
    return (0 if lo <= 0 <= hi else min(abs(lo), abs(hi))), max(abs(lo), abs(hi))


def as_decimal(x):
    x = Fraction(x)
    with localcontext() as context:
        context.prec = 12
        return Decimal(x.numerator) / Decimal(x.denominator)


@cache
def pade_coefficients(n):
    k = 2 * n + 1
    leading = comb(2 * k, k)
    q = [Fraction((-1)**j * comb(k, j) * comb(2 * k - 2 * j, k), leading)
         for j in range(n + 1)]
    p = [sum(q[j] / (2 * (i - j) + 1) for j in range(i + 1))
         for i in range(n + 1)]
    q_den = lcm(*(coefficient.denominator for coefficient in q))
    assert q_den == leading // 2**k.bit_count()
    norm = Fraction(2**(2 * k), (2 * k + 1) * leading**2)
    return p, q, lcm(q_den, *(coefficient.denominator for coefficient in p)), norm


@cache
def arm(m, n):
    p, q, multiplier, _ = pade_coefficients(n)
    height = multiplier * sum(c * (-1)**j * m**(2 * n + 1 - 2 * j) for j, c in enumerate(q))
    numerator = multiplier * sum(c * (-1)**j * m**(2 * n - 2 * j) for j, c in enumerate(p))
    assert height.denominator == numerator.denominator == 1
    return height.numerator, numerator.numerator


def row(n, r):
    a, b = arm(5, n)
    c, d = arm(239, r)
    return a * c, 16 * c * b - 4 * a * d


def determinant(v, w):
    return v[0] * w[1] - w[0] * v[1]


def valuation_two(value):
    value = Fraction(value)
    assert value != 0
    numerator, denominator = abs(value.numerator), value.denominator
    return (numerator & -numerator).bit_length() - (denominator & -denominator).bit_length()


def check_even_valuations():
    checks = 0
    even_orders = range(2, 33, 2)
    for m in (5, 239):
        rationals = {}
        for n in range(2, 33):
            height, numerator = arm(m, n)
            rationals[n] = Fraction(numerator, height)
            if n % 2:
                even = n - 1
                multiplier = pade_coefficients(n)[2]
                q_value = Fraction(height, multiplier * m**(2 * n + 1))
                assert valuation_two(rationals[n] - rationals[even]) == (
                    3 * even - 2 * even.bit_count() - valuation_two(q_value)
                )
        for n in even_orders:
            for next_n in range(n + 2, 33, 2):
                assert valuation_two(rationals[next_n] - rationals[n]) == (
                    2 * (n - n.bit_count())
                )
                checks += 1
    print("Exact even-order differences across both arms:", checks)


def check_odd_spacing():
    count = 0
    for m in (5, 239):
        values = [Fraction(1), Fraction(m)]
        for j in range(1, 4 * 40 + 3):
            values.append(m * values[-1] + Fraction(j * j, 4 * j * j - 1) * values[-2])
        ratios = [values[4 * k + 3] / (4 * values[4 * k]) for k in range(41)]
        for k, ratio in enumerate(ratios):
            assert valuation_two(ratio) >= 1
            for j in range(k + 1, len(ratios)):
                assert valuation_two(ratios[j] - ratio) == 1 + valuation_two(j - k)
                count += 1
    print("Exact odd-order residue-spacing identities:", count)


def common_convergents(interval, count):
    lo, hi = interval
    p, q, old_p, old_q = 1, 0, 0, 1
    conv = []
    for _ in range(count):
        digit = lo.numerator // lo.denominator
        assert digit == hi.numerator // hi.denominator and lo != digit
        p, old_p = digit * p + old_p, p
        q, old_q = digit * q + old_q, q
        conv.append((p, q))
        lo, hi = 1 / (hi - digit), 1 / (lo - digit)
    return conv


def check_convergents(interval):
    conv = common_convergents(interval, 42)
    passing = []
    for j in range(1, 36):
        p, q = conv[j]
        next_p, next_q = conv[j + 1]
        current = error_bounds(q, p, interval)
        following = error_bounds(next_q, next_p, interval)
        assert abs(p * next_q - next_p * q) == 1
        assert 2 * current[0] * next_q > 1 and current[1] * next_q < 1
        if 2 * current[1] * q <= 1:
            assert 2 * following[1] * q <= 1
            passing.append((j, Fraction(next_q, q)))
        else:
            assert 2 * current[0] * q > 1
    assert len(passing) == 26
    print(f"Common continued fraction; small-scale pairs {len(passing)}/35")
    print("Passing indices:", [n for n, _ in passing])
    print("Required C range at those small scales:", min(c for _, c in passing),
          max(c for _, c in passing))

    p, q = conv[3]
    next_p, next_q = conv[4]
    e = error_bounds(q, p, interval)
    u = error_bounds(next_q, next_p, interval)
    k = 146
    second_height = next_q - k * q
    lower = 2 * (u[0] + (k + 1) * e[0]) * second_height
    upper = 2 * (u[1] + (k + 1) * e[1]) * second_height
    assert Fraction(1477, 10) < lower < upper < Fraction(1478, 10)
    print("Uniform C necessary from j=3:", as_decimal(lower))


def main():
    interval = pi_interval()
    max_order = 12

    for m in (5, 239):
        for n in range(1, max_order):
            _, _, multiplier, norm = pade_coefficients(n)
            _, _, next_multiplier, _ = pade_coefficients(n + 1)
            assert abs(determinant(arm(m, n), arm(m, n + 1))) == (
                m * multiplier * next_multiplier * norm
            )
    print("Exact adjacent Padé determinant identities:", 2 * (max_order - 1))
    check_even_valuations()
    check_odd_spacing()

    assert len(str(abs(determinant(row(4, 2), row(5, 2))))) == 42
    best = None
    for n in range(1, max_order + 1):
        for r in range(1, max_order + 1):
            height, numerator = row(n, r)
            divisor = gcd(height, numerator)
            height //= divisor
            numerator //= divisor
            lo, hi = error_bounds(height, numerator, interval)
            lo, hi = height * lo, height * hi
            assert lo > 1
            if best is None or hi < best[0]:
                best = lo, hi, n, r
            else:
                assert lo > best[1]
    lo, hi, n, r = best
    assert hi - lo < lo / 10**10
    print(f"Primitive Padé grid: {max_order**2} rows all have q*abs(q*pi-p) > 1")
    print(f"Smallest primitive defect: indices ({n},{r}), approx {as_decimal(lo)}")
    check_convergents(interval)


if __name__ == "__main__":
    main()
