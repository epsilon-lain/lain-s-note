#!/usr/bin/env python3
"""Independent exact tests of the two-index prime-overlap certificate."""

import json
from math import gcd, prod
from pathlib import Path


def primes_to(n):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = b"\x00" * (((n - p * p) // p) + 1)
    return [p for p in range(3, n + 1, 2) if sieve[p]]


def cube_floor(n):
    b = 0
    while (b + 1) ** 3 <= n:
        b += 1
    return b


def run(n_max=400, h_max=30):
    limit = n_max + h_max
    u = [1, 2]
    for n in range(1, limit):
        numerator = 2 * (2 * n + 1) * u[n] + 4 * n * u[n - 1]
        assert numerator % (n + 1) == 0
        u.append(numerator // (n + 1))
    primes = primes_to(limit)
    pp = [1] * (limit + 1)
    theta_product = [1] * (limit + 1)
    selected = [set() for _ in range(limit + 1)]
    for n in range(1, limit + 1):
        b = cube_floor(n)
        selected[n] = {
            p for p in primes
            if p <= n and p * p > n and p * (b + 1) > n and u[n] % p == 0
        }
        pp[n] = prod(selected[n])
        theta_product[n] = prod(p for p in primes if p <= n)
        assert u[n] <= 6**n

    pair_checks = adjacent_exception_checks = 0
    for n in range(1, n_max + 1):
        k_product = prod(u[1:cube_floor(n) + 1])
        for p in selected[n]:
            if u[n - 1] % p == 0:
                assert u[n // p] % p == 0
                assert k_product % p == 0
                adjacent_exception_checks += 1
        f_previous, f_current = 1, 2 * (2 * n + 1)
        g_previous, g_current = 0, 4 * n
        d = 1
        for h in range(1, h_max + 1):
            if h > 1:
                k = h - 1
                f_previous, f_current = (
                    f_current,
                    2 * (2 * n + 2 * k + 1) * f_current
                    + 4 * (n + k)**2 * f_previous,
                )
                g_previous, g_current = (
                    g_current,
                    2 * (2 * n + 2 * k + 1) * g_current
                    + 4 * (n + k)**2 * g_previous,
                )
            m = n + h
            d *= m
            assert d * u[m] == f_current * u[n] + g_current * u[n - 1]
            assert 0 < g_current <= (5 * m)**h
            common = gcd(pp[n], pp[m])
            assert (g_current * k_product) % common == 0
            assert (theta_product[m] * g_current * k_product) % (pp[n] * pp[m]) == 0
            pair_checks += 1
    return {
        "n_max": n_max,
        "h_max": h_max,
        "pair_checks": pair_checks,
        "adjacent_exception_checks": adjacent_exception_checks,
        "status": "PASS",
        "scope": "Finite exact identity checks; the asymptotic theorem requires the written proof.",
    }


if __name__ == "__main__":
    result = run()
    Path(__file__).with_name("spacing_certificate_checks.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))
