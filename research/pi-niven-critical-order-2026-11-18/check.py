"""Exact checks for the beta balancing and rational midpoint correction."""

from fractions import Fraction
from math import factorial, gcd
from pathlib import Path
import json
import runpy


def choose_degree(p, q, eta):
    """Find the nearest degree on the logarithmic scale by integer comparisons."""
    r = Fraction(p, q)
    target_squared = q*q / eta
    n, a = 1, 2/r
    while True:
        below = (2*n+1)*a*a
        next_a = (4*n+2)*a/r
        above = (2*n+3)*next_a*next_a
        if above >= target_squared:
            return n if target_squared**2 <= below*above else n+1
        n, a = n+1, next_a


def midpoint_terms(p, q, n, m):
    r = Fraction(p, q)
    term = q**n * r**(2*n+1) * factorial(n) / factorial(2*n+1)
    for j in range(m+2):
        yield term
        term *= r*r / (8*(j+1)*(2*n+2*j+3))


def endpoints(p, q, degree):
    """Evaluate H_n both from the integer recurrence and its coefficient expansion."""
    older, current = (1, 0), (2*q, -p)
    for n in range(2, degree+1):
        older, current = current, tuple(
            2*(2*n-1)*q*a-p*p*b for a, b in zip(current, older)
        )
    real = sum(
        (-1)**(k//2) * (factorial(2*degree-k)//(factorial(k)*factorial(degree-k)))
        * p**k * q**(degree-k)
        for k in range(0, degree+1, 2)
    )
    if current[0] != real:
        raise ArithmeticError("recurrence and coefficient expansion disagree")
    return current


def integer_grid_checks():
    count = 0
    for q in range(2, 51):
        for p in range(3*q, 16*q//5+1):
            if gcd(p, q) != 1:
                continue
            for n in (1, 3, 6):
                for m in range(4):
                    terms = list(midpoint_terms(p, q, n, m))
                    correction = sum((-1)**j*t for j, t in enumerate(terms[:-1])) / 2
                    u, _ = endpoints(p, q, n)
                    d = (u-correction).denominator
                    if d % q**(n+2*m+1) or terms[-1]*2 > terms[-2]:
                        raise ArithmeticError("incorrect denominator or tail bound")
                    count += 1
    return count


def record_checks():
    previous = Path(__file__).resolve().parent.parent / "pi-first-hit-window-2026-11-18"
    data = json.loads((previous / "certificate.json").read_text())
    verifier = runpy.run_path(str(previous / "verify.py"))
    verifier["verify"](data)
    for row in data["records"][1:]:
        p, q = row["p"], row["q"]
        # Take the certified error upper bound as a concrete test threshold.
        eta = Fraction(row["cost_ceil_1e15"], 10**15)
        n = choose_degree(p, q, eta)
        u, v = endpoints(p, q, n)
        r = Fraction(p, q)
        height = q**n * factorial(2*n) // factorial(n)
        scale = p**(n-1)
        abs_ratio = Fraction(abs(u), scale)
        if not abs_ratio**2 < 20**2*eta or not abs(v) < 2*height:
            raise ArithmeticError("critical-scale inequality failed")
        a = Fraction(factorial(2*n), factorial(n)) / r**n
        beta_bound = r / ((2*n+1)*a*a)
        if not abs(u) <= height*(beta_bound+2*eta/(q*q)):
            raise ArithmeticError("beta inequality failed")
        print(f"CF {row['n']}: degree {n}, floor(10^9 |U|/p^(n-1)) = {abs_ratio*10**9//1}")


if __name__ == "__main__":
    print(f"Denominator and moment ratio checked in {integer_grid_checks()} cases.")
    record_checks()
