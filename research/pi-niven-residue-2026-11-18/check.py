"""Exact arithmetic for the rational-endpoint test. No floating point or pi input."""
from fractions import Fraction
from math import factorial, gcd
from pathlib import Path
import json


def degree(q):
    n = 2
    fac = factorial(4)
    while fac <= 4 * 3**n * q**3:
        fac *= (2*n+1)*(2*n+2)*(2*n+3)*(2*n+4)
        n += 2
    return n


def endpoint(p, q, n):
    older = (1, 0)
    previous = (2*q, -p)
    for j in range(2, n+1):
        factor = 2*(2*j-1)*q
        current = tuple(factor*a-p*p*b for a, b in zip(previous, older))
        if current[1]*previous[0]-current[0]*previous[1] != -p**(2*j-1):
            raise ArithmeticError("Wronskian identity failed")
        older, previous = previous, current
    return previous


def split(p, q, n):
    k = max(1, n//4)
    left_mod, right_mod = p**(2*k), q**(n-2*k+2)
    low = high = 0
    for j in range(n//2+1):
        coeff = factorial(2*n-2*j)//(factorial(2*j)*factorial(n-2*j))
        term = (-1)**j * coeff * p**(2*j) * q**(n-2*j)
        if j < k:
            low = (low+term) % left_mod
        else:
            high = (high+term) % right_mod
    combined = low + left_mod*((high-low)*pow(left_mod,-1,right_mod) % right_mod)
    return combined, left_mod*right_mod


def test_pair(p, q):
    if q < 2 or gcd(p,q) != 1 or not 3*q <= p <= 16*q//5:
        raise ValueError("expected a reduced fraction in [3,16/5] with q >= 2")
    n = degree(q)
    u, v = endpoint(p, q, n)
    residue, modulus = split(p, q, n)
    t = q**n * (factorial(2*n)//factorial(n))
    if not (residue == u % modulus and gcd(u,q) == 1):
        raise ArithmeticError("residue test failed")
    centered = min(residue, modulus-residue)
    if 2*abs(u)<modulus and centered != abs(u):
        raise ArithmeticError("center lift failed")
    return n, Fraction(q*q*abs(u),t), Fraction(abs(v),t), (2*abs(u)<modulus)


if __name__ == "__main__":
    path = Path(__file__).resolve().parent.parent / "pi-first-hit-window-2026-11-18/certificate.json"
    previous = json.loads(path.read_text())
    count = 0
    for q in range(2,161):
        for p in range(3*q,16*q//5+1):
            if gcd(p,q)==1:
                test_pair(p,q)
                count += 1
    print(f"Grid verified: {count} reduced rational endpoints")
    for row in previous["records"][1:]:
        p, q = row["p"], row["q"]
        n, score, imag, centered = test_pair(p,q)
        grid = 10**12
        old_lo = Fraction(row["cost_floor_1e15"],10**15)
        old_hi = Fraction(row["cost_ceil_1e15"],10**15)
        if not (old_hi/5-Fraction(1,q) <= score <= 2*old_lo+Fraction(1,q)
                and Fraction(4,5)<imag<2):
            raise ArithmeticError("comparison with independent pi certificate failed")
        print(f"n_CF={row['n']}, order={n}, floor(10^12 S)={score*grid//1}, centered={centered}")
