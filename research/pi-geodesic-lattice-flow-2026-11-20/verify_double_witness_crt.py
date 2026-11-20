#!/usr/bin/env python3
"""Exact checks of cross-degree divisibility and the CRT mixed slopes."""

from decimal import Decimal, localcontext
from fractions import Fraction
from math import comb, factorial, gcd, isqrt, prod


def slopes(max_n, max_j):
    remainders = [Fraction(0)]
    for m in range(1, 2 * max_n + max_j + 1):
        remainders.append(Fraction(1, 2 * m - 1) - remainders[-1])

    def slope(n, j):
        rational_part = sum(
            (-1)**r * comb(n, r) * remainders[n + j + r]
            for r in range(n + 1)
        )
        return (-1)**(n + j + 1) * Fraction(4, 2**n) * rational_part

    return slope


def valuation_two(value):
    value = abs(value)
    assert value
    return (value & -value).bit_length() - 1


def residue(value, prime):
    return value.numerator * pow(value.denominator, -1, prime) % prime


def primes_through(limit):
    return [
        candidate for candidate in range(2, limit + 1)
        if all(candidate % divisor for divisor in range(2, isqrt(candidate) + 1))
    ]


def determinant(a, b):
    return a.numerator * b.denominator - b.numerator * a.denominator


def crt_pair(n, a, b, prime_pool, center):
    exponent = 2 + n - n.bit_count()
    delta = b - a
    product = prod(prime_pool)
    signed_sum = (1 << exponent) * product
    weight = 0
    modulus = 1
    for prime in prime_pool:
        assert a.denominator % prime == 0 and a.denominator % prime**2 != 0
        assert residue(prime * a, prime) == residue(prime * b, prime) != 0
        target = -(signed_sum // prime) * residue(prime * a, prime)
        target *= pow(residue(delta, prime), -1, prime)
        weight += modulus * ((target - weight) * pow(modulus, -1, prime) % prime)
        modulus *= prime
    assert modulus == product
    shift = round((signed_sum * center - weight) / product)
    weight += shift * product
    weights = (weight, weight + product)
    u, v = (a + Fraction(t, signed_sum) * delta for t in weights)
    assert all(0 < t < signed_sum for t in weights)
    assert v - u == delta / (1 << exponent)
    assert gcd(u.denominator * v.denominator, 2 * product) == 1
    cross = abs(determinant(u, v))
    assert gcd(cross, 2) == 1
    return u, v, cross


def beta_center(n, exponent):
    """Lower bound for I_(n,0)/(I_(n,0)+I_(n,1)), with error < 2^-exponent."""
    moment = Fraction(1)
    center = Fraction(0)
    for j in range(exponent):
        center += moment / 2**(j + 1)
        moment *= Fraction(2 * (n + j + 1), 4 * n + 2 * j + 3)
    return center


def five_interval_checks(slope):
    checked = 0
    for n in (20, 32, 48, 64, 80, 96):
        for prime in primes_through(2 * n):
            a, h = divmod(n, prime)
            if a > 2 or prime <= 5 or prime**2 <= 4 * n + 15:
                continue
            if 4 * h + 16 < prime and a > 0:
                block = 0
            elif 2 * h > prime and 4 * h + 16 < 3 * prime:
                block = 1
            else:
                continue
            contribution = sum(
                comb(a, v) * sum(
                    Fraction((-1)**ell, 2 * ell + 1)
                    for ell in range(a + v + block)
                )
                for v in range(a + 1)
            )
            expected = Fraction(-4 * (-1)**((prime + 1) // 2), 2**a)
            expected *= contribution
            for j in range(9):
                assert residue(prime * slope(n, j), prime) == residue(expected, prime)
                checked += 1
    return checked


def arctan_bounds(base, terms=400):
    total = sum(
        Fraction((-1)**i, (2 * i + 1) * base**(2 * i + 1))
        for i in range(terms)
    )
    end = total + Fraction((-1)**terms, (2 * terms + 1) * base**(2 * terms + 1))
    return min(total, end), max(total, end)


def height_error_bounds(a, b, pi_bounds):
    pi_lo, pi_hi = pi_bounds
    lower, upper = [], []
    for value in (a, b):
        lo = value.denominator * pi_lo - value.numerator
        hi = value.denominator * pi_hi - value.numerator
        lower.append(0 if lo <= 0 <= hi else min(abs(lo), abs(hi)))
        upper.append(max(abs(lo), abs(hi)))
    height = max(a.denominator, b.denominator)
    return height * max(lower), height * max(upper)


def scientific(value):
    with localcontext() as context:
        context.prec = 12
        return f"{Decimal(value.numerator) / Decimal(value.denominator):.6E}"


def main():
    max_n, max_j = 96, 8
    slope = slopes(max_n, max_j)
    fixed = (slope(2, 0), slope(2, 1))
    assert abs(determinant(*fixed)) == 840
    mixes = tuple(
        (a * fixed[0] - b * fixed[1]) / (a - b)
        for a, b in [(-20, 1), (-61, 3)]
    )
    assert mixes == (Fraction(20672, 6615), Fraction(25, 8))
    assert abs(determinant(*mixes)) == 1

    for n in range(2, max_n + 1):
        exponent = 2 + n - n.bit_count()
        for j in range(max_j):
            a, b = slope(n, j), slope(n, j + 1)
            assert a.denominator % 2 and b.denominator % 2
            change = Fraction((-1)**(n + j) * 4 * factorial(n), prod(
                2 * (n + j + r) + 1 for r in range(n + 1)
            ))
            assert b - a == change
            order = valuation_two(a.numerator)
            if n % 2:
                assert order == exponent - 1
            else:
                assert order >= exponent
        if n < max_n:
            assert slope(n + 1, 0) == (slope(n, 1) + slope(n, 2)) / 2

    a_lo, a_hi = arctan_bounds(5)
    b_lo, b_hi = arctan_bounds(239)
    pi_bounds = 16 * a_lo - 4 * b_hi, 16 * a_hi - 4 * b_lo
    for n in (2, 3, 4, 6, 8, 12, 24, 48, 72, 96):
        primes = primes_through(4 * n + 1)
        pool = [p for p in primes if 3 * p > 4 * n + 1 and p <= 2 * n - 1]
        second_pool = [p for p in primes if p > 5 and 5 * p > 4 * n + 3 and p <= n]
        a, b = slope(n, 0), slope(n, 1)
        exponent = 2 + n - n.bit_count()
        center = beta_center(n, exponent + 2)
        delta = b - a
        simple_poles = [
            p for p in primes if p > 2
            and a.denominator % p == 0 and a.denominator % p**2 != 0
            and delta.denominator % p != 0 and delta.numerator % p != 0
        ]
        assert set(pool + second_pool) <= set(simple_poles)
        _, _, first_cross = crt_pair(n, a, b, pool, center)
        assert first_cross % prod(second_pool) == 0
        u, v, cross = crt_pair(n, a, b, simple_poles, center)
        assert gcd(cross, 2 * prod(simple_poles)) == 1
        assert min(a, b) < min(u, v) < max(u, v) < max(a, b)
        assert max(abs(endpoint - z) for endpoint in pi_bounds for z in (u, v)) <= (
            2 * abs(delta) / 2**exponent
        )
        lower, upper = height_error_bounds(u, v, pi_bounds)
        assert 2 * upper >= cross and lower > 0
        assert upper - lower < lower / 10**10
        print(
            f"n={n:2d} primes_removed={len(simple_poles):2d} "
            f"det={scientific(Fraction(cross))} HE≈{scientific(lower)}"
        )
    print(f"Five-interval residue identities: {five_interval_checks(slope)}")
    print("PASS exact adjacent differences, cross-degree recurrence, and CRT reductions")


if __name__ == "__main__":
    main()
