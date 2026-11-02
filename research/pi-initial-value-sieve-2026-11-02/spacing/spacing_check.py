#!/usr/bin/env python3
"""Exact checks for a common-large-prime quantitative spacing certificate."""
from math import prod, isqrt


def sequence(limit):
    u = [1, 2]
    for n in range(1, limit):
        numerator = 2 * (2 * n + 1) * u[n] + 4 * n * u[n - 1]
        assert numerator % (n + 1) == 0
        u.append(numerator // (n + 1))
    return u


def primes(limit):
    flags = bytearray(b"\x01") * (limit + 1)
    flags[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p * p::p] = b"\x00" * ((limit - p * p) // p + 1)
    return [p for p in range(2, limit + 1) if flags[p]]


def transfer(n, h):
    if h == 0:
        return 1, 0
    f0, g0 = 1, 0
    f1, g1 = 2 * (2 * n + 1), 4 * n
    for k in range(1, h):
        alpha, beta = 2 * (2 * n + 2 * k + 1), 4 * (n + k) ** 2
        f0, f1 = f1, alpha * f1 + beta * f0
        g0, g1 = g1, alpha * g1 + beta * g0
    return f1, g1


def main():
    u = sequence(1100)
    for n in range(1, 101):
        for h in range(1, 41):
            f, g = transfer(n, h)
            assert prod(range(n + 1, n + h + 1)) * u[n + h] == f * u[n] + g * u[n - 1]
            assert 0 < g <= (5 * (n + h)) ** h

    ps = primes(1100)
    for n in range(2, 1001):
        for p in ps:
            if p == 2 or p > n:
                continue
            if p * p <= n:
                continue
            if u[n] % p == 0 and u[n - 1] % p == 0:
                assert u[n // p] % p == 0

    # Exact product certificate; use p>sqrt(n), a wider band than required.
    for n in range(2, 250):
        for h in range(1, 16):
            m = n + h
            common = [p for p in ps if p != 2 and p <= n and p * p > m
                      and u[n] % p == 0 and u[m] % p == 0]
            good = [p for p in common if u[n - 1] % p]
            bad = [p for p in common if u[n - 1] % p == 0]
            _, g = transfer(n, h)
            assert g % prod(good) == 0
            k = prod(u[b] for b in range(1, isqrt(n) + 1))
            assert k % prod(bad) == 0
    print("PASS: transfer identities, height, adjacent-prefix exceptions, common-prime certificates")


if __name__ == "__main__":
    main()
