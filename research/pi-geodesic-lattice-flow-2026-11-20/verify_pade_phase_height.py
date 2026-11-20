#!/usr/bin/env python3
"""Exact checks for the normalized diagonal Padé gcd and phase height."""
from math import comb, gcd, isqrt, lcm, prod


def coefficients(n):
    return [
        comb(2*n-k, k) * prod(range(1, 2*(n-k), 2))
        for k in range(n+1)
    ]


def direct_value(n, p, q):
    real = imag = 0
    for k, coefficient in enumerate(coefficients(n)):
        term = coefficient * q**(n-k) * p**k
        if k % 4 == 0:
            real += term
        elif k % 4 == 1:
            imag += term
        elif k % 4 == 2:
            real -= term
        else:
            imag -= term
    return real, imag


def prime_factors(n):
    divisors = []
    for candidate in range(2, isqrt(n) + 1):
        if n % candidate:
            continue
        divisors.append(candidate)
        while n % candidate == 0:
            n //= candidate
    return divisors + ([n] if n > 1 else [])


def valuation(n, prime):
    exponent = 0
    while n % prime == 0:
        n //= prime
        exponent += 1
    return exponent


def check_pair(p, q, max_order):
    a_old, b_old = 1, 0
    a, b = q, p
    odd_factorial = 1
    factors = prime_factors(p)
    for n in range(1, max_order + 1):
        if n > 1:
            a, a_old = (2*n-1)*q*a - p*p*a_old, a
            b, b_old = (2*n-1)*q*b - p*p*b_old, b
        odd_factorial *= 2*n-1
        if n <= 6:
            assert (a, b) == direct_value(n, p, q)
        assert b*a_old - a*b_old == p**(2*n-1)

        content = gcd(a, b)
        assert p**(2*n-1) % content == 0
        a0, b0 = a//content, b//content
        extra = int(a0 % 2 == 1 and b0 % 2 == 1)
        assert extra == int(p % 2 == q % 2 == 1 and n % 3 == 1)

        # After the primitive Gaussian correction the coordinate denominator
        # equals the denominator norm. Compute both using rational integers.
        norm = a0*a0 + b0*b0
        coordinate_denominator = lcm(
            norm//gcd(norm, 2*b0*b0),
            norm//gcd(norm, 2*a0*b0),
        )
        denominator_norm = norm // (2**extra)
        assert coordinate_denominator == denominator_norm

        for prime in factors:
            if prime*prime > 2*n-1 and prime != 2:
                assert valuation(content, prime) == valuation(odd_factorial, prime)
            if prime != 2 and prime**valuation(p, prime) > 2*n-1:
                assert valuation(content, prime) == valuation(odd_factorial, prime)


def main():
    pairs = [(p, q) for p in range(1, 81) for q in range(1, 51) if gcd(p, q) == 1]
    pairs.extend([(355, 113), (104348, 33215)])
    for p, q in pairs:
        check_pair(p, q, 30)
    print(f"{len(pairs) * 30} exact recurrence, gcd and denominator cases passed")

    real, imag = direct_value(5, 3, 1)
    assert (real, imag) == (-1620, 243)
    assert gcd(real, imag) == 81
    assert prod(range(1, 10, 2)) == 945
    print("Constant-coefficient upper bound fails at (n,p,q)=(5,3,1)")


if __name__ == "__main__":
    main()
