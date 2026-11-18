"""Verify exact Padé coefficients, residue bounds, and signed tails."""

import json
from fractions import Fraction
from math import factorial, gcd, isqrt, prod
from pathlib import Path


def denominators(n, count):
    values = [1]
    for k in range(1, count + 1):
        values.append(values[-1] * 8*k*(2*n + 2*k + 1))
    return values


def determinant(rows):
    """Compute an integer determinant with fraction-free elimination."""
    matrix = [list(row) for row in rows]
    size = len(matrix)
    sign = previous = 1
    for column in range(size - 1):
        pivot_row = next((i for i in range(column, size) if matrix[i][column]), None)
        if pivot_row is None:
            return 0
        if pivot_row != column:
            matrix[pivot_row], matrix[column] = matrix[column], matrix[pivot_row]
            sign = -sign
        pivot = matrix[column][column]
        for i in range(column + 1, size):
            for j in range(column + 1, size):
                numerator = pivot*matrix[i][j] - matrix[i][column]*matrix[column][j]
                matrix[i][j], rem = divmod(numerator, previous)
                if rem:
                    raise ArithmeticError("Nonexact Bareiss division")
            matrix[i][column] = 0
        previous = pivot
    return sign * matrix[-1][-1]


def pade(n, m):
    """Return integral numerator and denominator for F_n(x), with common scale."""
    h = denominators(n, 2*m)
    rows = [[h[k]//h[k-j] for j in range(m+1)] for k in range(m+1, 2*m+1)]
    minors = [determinant([row[:j] + row[j+1:] for row in rows]) for j in range(m+1)]
    if not all(value > 0 for value in minors):
        raise ArithmeticError("Expected positive Schur minors")
    b = [h[m]*value for value in minors]
    a = [(-1)**k * sum((-1)**j * minors[j] * (h[m]//h[k-j])
                      for j in range(k+1))
         for k in range(m+1)]
    if (-1)**m * a[-1] != determinant(
        [[h[k]//h[k-j] for j in range(m+1)] for k in range(m, 2*m+1)]
    ):
        raise ArithmeticError("Top numerator determinant mismatch")
    if not a[-1]:
        raise ArithmeticError("Numerator lost degree")
    return a, b


def homogenize(coefficients, p, q):
    degree = len(coefficients) - 1
    return sum(c * p**(2*j) * q**(2*(degree-j))
               for j, c in enumerate(coefficients))


def check(p, q, n, m):
    if gcd(p, q) != 1 or not 3*q <= p <= 16*q//5:
        raise ValueError("Need a reduced endpoint in [3, 16/5]")
    a, b = pade(n, m)
    h = denominators(n, 2*m+5)
    e = [Fraction(1, value) for value in h]
    residual = []
    for k in range(2*m+6):
        coeff = sum(Fraction((-1)**(k-j)*b[j], h[k-j])
                    for j in range(min(k,m)+1))
        if k <= m:
            coeff -= a[k]
        residual.append(coeff)
    if any(residual[:2*m+1]) or residual[2*m+1] == 0:
        raise ArithmeticError("Incorrect Padé contact order")
    for left, right in zip(residual[2*m+1:], residual[2*m+2:]):
        if not left*right < 0 or abs(right) > e[1]*abs(left):
            raise ArithmeticError("Signed tail check failed")

    x = Fraction(p*p, q*q)
    p_value, q_value = homogenize(a,p,q), homogenize(b,p,q)
    rational_value = Fraction(p_value,q_value)
    t0 = Fraction(p**(2*n+1)*factorial(n), q**(n+1)*factorial(2*n+1))
    correction = t0 * rational_value / 2
    d = correction.denominator
    # Factorization-free, weaker integer consequence of the valuation theorem.
    lead = abs(a[-1])
    if d*d * lead**(n+3) < q**(2*n+2):
        raise ArithmeticError("Reduced denominator bound failed")

    denominator_at_x = Fraction(q_value, q**(2*m))
    tail_at_x = abs(residual[2*m+1]) * x**(2*m+1)
    error_lower = t0 * tail_at_x * (1-x*e[1]) / (2*denominator_at_x)
    return d, error_lower


def inverse_corrections():
    """Check exact power-of-two cancellation in Padé [0/b]."""
    count = 0
    for n in range(1, 7):
        order = (n+2)//2
        h = denominators(n, order)
        inverse = [Fraction(1)]
        for k in range(1, order+1):
            inverse.append(-sum(Fraction((-1)**j, h[j])*inverse[k-j]
                                for j in range(1, k+1)))
        for exponent in range(3, 8):
            q = 2**exponent
            p = 3*q+1
            x = Fraction(p*p, q*q)
            t0 = Fraction(p**(2*n+1)*factorial(n), q**(n+1)*factorial(2*n+1))
            correction = (t0/2) / sum(term*x**j for j, term in enumerate(inverse))
            numerator = abs(correction.numerator)
            v2_factorial = order-order.bit_count()
            expected = (2*order-n-1)*exponent+3*order+v2_factorial-(n+1)
            if correction.denominator % 2 != 1 or (numerator & -numerator) != 2**expected:
                raise ArithmeticError("Incorrect two-adic power in reduced correction")
            count += 1
    return count


if __name__ == "__main__":
    print(f"Exact power-of-two cancellation confirmed in {inverse_corrections()} cases.")
    cases = ((199, 64), (355, 113), (3553, 1131), (314159, 100000))
    total = 0
    for p, q in cases:
        for n in (1, 2, 4, 9, 16, 25):
            for m in range(1, isqrt(n)+1):
                check(p,q,n,m)
                total += 1
    print(f"Exact tests passed: {total} degree and fraction combinations")
    first_hit = Path(__file__).resolve().parent.parent / "pi-first-hit-window-2026-11-18/certificate.json"
    witness = json.loads(first_hit.read_text())["records"][-1]
    p, q, n = witness["p"], witness["q"], 94
    x = Fraction(p*p, q*q)
    sigma = Fraction(1, 8*(2*n+3))
    t0 = Fraction(p**(2*n+1)*factorial(n), q**(n+1)*factorial(2*n+1))
    height = q**n * factorial(2*n) // factorial(n)
    w_lower = Fraction(witness["cost_floor_1e15"], 10**15)
    drift_lower = height*w_lower/(10*q*q)
    if 100 * t0*(x*sigma)**3 / 2 >= drift_lower:
        raise ArithmeticError("All-order signed dominance was not certified")
    print("All Padé orders m>=1: correction < 1/100 of endpoint drift.")
    for m in (1, 3, 6, 9, 18, 30):
        d, delta = check(witness["p"], witness["q"], 94, m)
        lower = d * delta
        # Find an integer exponent using a bit-length estimate, then verify it exactly.
        power = max(0, (lower.numerator.bit_length() - lower.denominator.bit_length()) * 3 // 10 - 2)
        if lower <= 10**power:
            raise ArithmeticError("Reported lower bound was not certified")
        excess = max(0, (delta.numerator.bit_length() - delta.denominator.bit_length()) * 3 // 10 - 2)
        if delta <= 10**excess:
            raise ArithmeticError("Absolute error lower bound was not certified")
        print(f"155-digit denominator, n=94, m={m}: error > 10^{excess}; D*error > 10^{power}")
