#!/usr/bin/env python3
"""Exact elliptic test of CM-sparse formal-group cancellation weights.

Curve E: y^2=x^3-2, P=(3,5), formal parameter xi=-x/y.
All elliptic group operations and rational cancellation values use Fraction.
Only logarithmic summaries are floating-point.
"""
from fractions import Fraction
from math import comb, log
import json
from pathlib import Path

P0 = (Fraction(3), Fraction(5))
S = 41
SCHEMES = {
    "CM_sparse_J4_order12": {2: -364, 3: 291, 4: -36},
    "CM_sparse_J5_order12": {2: 28, 3: 18, 4: -36, 5: 7},
    "CM_sparse_J6_order12": {2: -19, 3: -4, 4: 13, 5: 1, 6: -1},
    "binomial_order4": {j: (-1)**j * comb(5, j-1) for j in range(2, 7)},
    "binomial_order6": {j: (-1)**j * comb(7, j-1) for j in range(2, 9)},
    "binomial_order7": {j: (-1)**j * comb(8, j-1) for j in range(2, 10)},
}


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and y1 == -y2:
        return None
    if P == Q:
        if y1 == 0:
            return None
        slope = Fraction(3 * x1 * x1, 2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope * slope - x1 - x2
    y3 = slope * (x1 - x3) - y1
    return x3, y3


def mul(n):
    result, base = None, P0
    while n:
        if n & 1:
            result = add(result, base)
        base = add(base, base)
        n >>= 1
    return result


def xi(Q):
    return -Q[0] / Q[1]


def log_int(n):
    n = abs(n)
    if n == 0:
        return float("-inf")
    bits = n.bit_length()
    shift = max(0, bits - 53)
    top = n >> shift
    return log(top) + shift * log(2)


def log10_fraction(x):
    x = abs(x)
    return (log_int(x.numerator) - log_int(x.denominator)) / log(10)


def log_height(x):
    return max(log_int(x.numerator), log_int(x.denominator))


def c7(j):
    return Fraction(6, 7) * (j**7 - j)


def c13(j):
    return Fraction(12, 91) * (4*j**13 - 39*j**7 + 35*j)


def main():
    needed = sorted({S} | {j*S for weights in SCHEMES.values() for j in weights})
    point_xi = {n: xi(mul(n)) for n in needed}
    z = point_xi[S]
    rows = []

    for name, weights in SCHEMES.items():
        S1 = sum(lam*j for j, lam in weights.items())
        S7 = sum(lam*j**7 for j, lam in weights.items())
        S13 = sum(lam*j**13 for j, lam in weights.items())
        coeff_z6 = sum((lam*c7(j) for j, lam in weights.items()), Fraction(0))
        coeff_z12 = sum((lam*c13(j) for j, lam in weights.items()), Fraction(0))
        value = sum((lam * point_xi[j*S] / z for j, lam in weights.items()), Fraction(0))
        error = value - 1
        first_order = 6 if coeff_z6 else 12
        leading_coeff = coeff_z6 if coeff_z6 else coeff_z12

        rows.append({
            "scheme": name,
            "weights": {str(j): lam for j, lam in weights.items()},
            "largest_multiplier": max(weights),
            "S1": S1,
            "S7": S7,
            "S13": S13,
            "C0": sum(abs(lam) for lam in weights.values()),
            "C2": sum(abs(lam)*j*j for j, lam in weights.items()),
            "formal_coefficient_z6": str(coeff_z6),
            "formal_coefficient_z12": str(coeff_z12),
            "leading_error_order": first_order,
            "leading_error_coefficient": str(leading_coeff),
            "predicted_log10_error_from_leading_term":
                log10_fraction(leading_coeff) + first_order*log10_fraction(z),
            "exact_rational_log10_abs_error": log10_fraction(error),
            "exact_rational_log_height": log_height(value),
            "numerator_bits": abs(value.numerator).bit_length(),
            "denominator_bits": value.denominator.bit_length(),
        })

    w4 = SCHEMES["CM_sparse_J4_order12"]
    assert sum(l*j for j, l in w4.items()) == 1
    assert sum(l*j**7 for j, l in w4.items()) == 1
    assert sum((l*c7(j) for j, l in w4.items()), Fraction(0)) == 0
    assert sum((l*c13(j) for j, l in w4.items()), Fraction(0)) == -1031184000

    result = {
        "curve": "E: y^2=x^3-2",
        "point": "P=(3,5)",
        "formal_parameter": "xi=-x/y",
        "s": S,
        "abs_xi_sP_log10": log10_fraction(z),
        "formal_group_expansion": {
            "multiplication_series": "[j](z)=jz+(6/7)(j^7-j)z^7+(12/91)(4j^13-39j^7+35j)z^13+O(z^19)",
            "CM_symmetry": "only powers z^(6k+1) occur in [j](z)"
        },
        "rows": rows,
        "scope": "All elliptic group operations and rational cancellation values are exact. Formal expansions are exact through z^13. This is one curve and one s-value; no theorem about pi follows."
    }
    out = Path(__file__).with_name("cm_sparse_cancellation_experiment.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
