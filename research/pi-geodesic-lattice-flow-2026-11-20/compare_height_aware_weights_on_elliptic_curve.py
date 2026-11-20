#!/usr/bin/env python3
"""Exact comparison of binomial and height-proxy-optimized weights.

Curve: E: y^2=x^3-2, P=(3,5), xi=-x/y.
All elliptic group arithmetic and all B_m values use fractions.Fraction.
Only logarithms of exact integers are floating-point summaries.
"""
from fractions import Fraction
from math import comb, log, log10
import json
from pathlib import Path

P0 = (Fraction(3), Fraction(5))

ALT = {
    4: {2: 4, 3: -5, 5: 5, 6: -4, 7: 1},
    5: {2: 4, 3: -5, 4: 1, 5: 1, 6: 1, 7: 1, 8: -5, 9: 4, 10: -1},
    6: {2: 4, 3: -4, 4: -3, 5: 6, 8: -6, 9: 3, 10: 4, 11: -4, 12: 1},
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
    R, Q = None, P0
    while n:
        if n & 1:
            R = add(R, Q)
        Q = add(Q, Q)
        n >>= 1
    return R


def xi(P):
    return -P[0] / P[1]


def binomial_weights(m):
    return {j: (-1) ** j * comb(m + 1, j - 1)
            for j in range(2, m + 3)}


def moment_vector(weights, m):
    return [sum(lam * j**r for j, lam in weights.items())
            for r in range(m + 1)]


def c2(weights):
    return sum(abs(lam) * j * j for j, lam in weights.items())


def c0(weights):
    return sum(abs(lam) for lam in weights.values())


def log10_abs_fraction(x):
    x = abs(x)
    return log10(x.numerator) - log10(x.denominator)


def log_height(x):
    return max(log(abs(x.numerator)), log(x.denominator))


def evaluate_B(s, weights, point_xi):
    z = point_xi[s]
    return sum((lam * point_xi[j * s] / z
                for j, lam in weights.items()), Fraction(0))


def main():
    s_values = (15, 26, 41)
    needed = set()
    for s in s_values:
        needed.add(s)
        for m in (4, 5, 6):
            needed.update(j * s for j in ALT[m])
            needed.update(j * s for j in binomial_weights(m))
    point_xi = {n: xi(mul(n)) for n in sorted(needed)}

    rows = []
    for s in s_values:
        for m in (4, 5, 6):
            for label, weights in (
                ("binomial", binomial_weights(m)),
                ("height_proxy_candidate", ALT[m]),
            ):
                moments = moment_vector(weights, m)
                assert moments == [1] * (m + 1), (s, m, label, moments)
                B = evaluate_B(s, weights, point_xi)
                error = abs(B - 1)
                rows.append({
                    "s": s,
                    "m": m,
                    "weights": weights,
                    "scheme": label,
                    "C0": c0(weights),
                    "C2": c2(weights),
                    "moments": moments,
                    "log10_abs_B_minus_1": log10_abs_fraction(error),
                    "log_height_B": log_height(B),
                    "numerator_bits": abs(B.numerator).bit_length(),
                    "denominator_bits": B.denominator.bit_length(),
                })

    # Pairwise comparison for each (s,m).
    comparisons = []
    by_key = {(r["s"], r["m"], r["scheme"]): r for r in rows}
    for s in s_values:
        for m in (4, 5, 6):
            base = by_key[(s, m, "binomial")]
            alt = by_key[(s, m, "height_proxy_candidate")]
            comparisons.append({
                "s": s,
                "m": m,
                "C2_reduction_percent":
                    100 * (base["C2"] - alt["C2"]) / base["C2"],
                "height_ratio_candidate_over_binomial":
                    alt["log_height_B"] / base["log_height_B"],
                "log10_error_difference_candidate_minus_binomial":
                    alt["log10_abs_B_minus_1"] - base["log10_abs_B_minus_1"],
                "candidate_has_lower_actual_height":
                    alt["log_height_B"] < base["log_height_B"],
                "candidate_has_smaller_error":
                    alt["log10_abs_B_minus_1"] < base["log10_abs_B_minus_1"],
            })

    result = {
        "curve": "y^2=x^3-2",
        "point": "(3,5)",
        "formal_parameter": "xi=-x/y",
        "arithmetic": "All group-law operations and B_m rational values are exact Fraction arithmetic.",
        "height_definition": "h(a/b)=log(max(|a|,b)) for reduced a/b, b>0.",
        "rows": rows,
        "comparisons": comparisons,
        "scope": "Finite experiment for s=15,26,41 and m=4,5,6. It tests these candidate weights only; it does not establish an asymptotic theorem or a result about pi.",
    }
    path = Path(__file__).with_name("height_aware_weights_elliptic_comparison.json")
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "rows": len(rows),
        "comparisons": comparisons,
        "output": str(path),
    }, indent=2))


if __name__ == "__main__":
    main()
