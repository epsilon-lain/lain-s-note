#!/usr/bin/env python3
"""Finite exact checks accompanying EDGE_CAPACITY_AND_MIDPOINT.md."""

from collections import defaultdict
from fractions import Fraction
from math import factorial
from pathlib import Path
import json


LIMIT = 2000


def primes_through(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            sieve[p*p::p] = b"\x00" * ((limit - p*p)//p + 1)
    return [p for p in range(2, limit + 1) if sieve[p]]


def trace(p):
    residues = {y*y % p for y in range(p)}
    order = 1
    for x in range(p):
        value = x*(x*x + 2*x + 2) % p
        order += 1 if value == 0 else 2 if value in residues else 0
    return p + 1 - order


def add(P, Q, p):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    z, w = Q
    if x == z:
        if (y + w) % p == 0:
            return None
        slope = (3*x*x + 4*x + 2) * pow(2*y, -1, p) % p
    else:
        slope = (w-y) * pow(z-x, -1, p) % p
    xx = (slope*slope - 2 - x - z) % p
    yy = (slope*(x-xx) - y) % p
    return xx, yy


def main():
    primes = primes_through(LIMIT)
    A = [1, 1]
    for r in range(1, LIMIT):
        A.append((2*r + 1)*A[-1] + r*r*A[-2])

    edge_pairs = 0
    for n in range(3, LIMIT + 1):
        fibers = defaultdict(list)
        for p in primes:
            if p == 2 or not n < 2*p or not p < n:
                continue
            k = min(n-p, 2*p-n-1)
            assert k >= 0
            assert 3*k <= n-1
            assert p == n-k or 2*p == n+1+k
            fibers[k].append(p)
            assert len(fibers[k]) <= 2
            # For any K >= k, one of these windows contains p. It suffices
            # to check their smallest such K=k, since the windows expand.
            assert n-k <= p < n or n < 2*p <= n+1+k
            edge_pairs += 1
        # Checking the count at each actual cutoff implies it at all cutoffs.
        count = 0
        for k in sorted(fibers):
            count += len(fibers[k])
            assert count <= 2*(k+1)

    midpoints = []
    four_torsion_checks = 0
    for p in primes:
        if p == 2:
            continue
        h = (p-1)//2
        U = Fraction((1 << h)*A[h], factorial(h))
        assert U.denominator == 1
        ap = trace(p)
        assert ap % p == U.numerator % p
        assert ap*ap <= 4*p
        if p >= 5:
            assert (A[h] % p == 0) == (ap == 0)
        if ap == 0:
            assert p % 8 == 7
            midpoints.append({"p": p, "n": (3*p-1)//2, "residue_mod32": p % 32})
        # This examines only x^2=2 candidates, then checks their actual doubles.
        square_roots = defaultdict(list)
        for y in range(p):
            square_roots[y*y % p].append(y)
        for x in square_roots[2 % p]:
            rhs = x*(x*x + 2*x + 2) % p
            for y in square_roots[rhs]:
                if y:
                    assert add((x, y), (x, y), p) == (0, 0)
                    formula_x = (x*x-2)**2 * pow(4*y*y, -1, p) % p
                    assert formula_x == 0
                    four_torsion_checks += 1

    report = {
        "limit": LIMIT,
        "status": "all exact assertions passed",
        "edge_prime_pairs_checked": edge_pairs,
        "odd_primes_with_independent_point_counts": len(primes)-1,
        "four_torsion_candidates_checked": four_torsion_checks,
        "midpoint_trace_zero_primes": midpoints,
        "scope": "Finite verification only. No global S(n) bound or mod32 theorem asserted.",
    }
    path = Path(__file__).with_name("edge_midpoint_verification.json")
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
