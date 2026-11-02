#!/usr/bin/env python3
"""Exact rational checks for the all-n balanced-clearing threshold.

Standard library only.  The accompanying note gives the uniform induction;
this script checks its finite premises and the displayed algebraic ratios.
No assertion about odd minimal content is used here.
"""

from fractions import Fraction
from functools import cache
import json
from math import factorial, log10
from pathlib import Path


@cache
def A(n):
    value = 1
    for k in range(n):
        value *= factorial(k)
    return value


T = Fraction(44, 7)
a = T * T / 16


def balanced_denominator(n):
    r, s = n // 2, (n + 1) // 2
    numerator = A(n) ** 2
    denominator = A(r) ** 2 * A(s) ** 2
    assert numerator % denominator == 0
    return numerator // denominator


def W(n):
    r, s = n // 2, (n + 1) // 2
    return T ** (n * n) * Fraction(
        A(n) ** 4, A(2 * n) * A(r) ** 2 * A(s) ** 2
    )


def exact_ratio_formula(n):
    r, s = n // 2, (n + 1) // 2
    return T ** (4 * n + 4) * Fraction(
        (factorial(n) * factorial(n + 1)) ** 4,
        factorial(2 * n)
        * factorial(2 * n + 1)
        * factorial(2 * n + 2)
        * factorial(2 * n + 3)
        * (factorial(r) * factorial(s)) ** 2,
    )


def b(n):
    r, s = n // 2, (n + 1) // 2
    return T ** 4 * (n + 1) ** 2 / 4 * (
        a ** n / (factorial(r) * factorial(s))
    ) ** 2


def b_ratio_formula(n):
    r, s = n // 2, (n + 1) // 2
    return Fraction((n + 3) ** 2, (n + 1) ** 2) * a ** 4 / (
        (r + 1) * (s + 1)
    ) ** 2


def log10_fraction(value):
    return log10(value.numerator) - log10(value.denominator)


def exact_comparison(label, value, bound):
    gap = bound.numerator * value.denominator - value.numerator * bound.denominator
    assert gap > 0
    return {
        "label": label,
        "bound": str(bound),
        "strict_inequality_certified_by_positive_cross_product": True,
        "cross_product_gap_decimal_digits": len(str(gap)),
        "log10_value_diagnostic_only": log10_fraction(value),
    }


def main():
    assert a == Fraction(121, 49) and a < 3
    algebra_checks = []
    for n in range(2, 41):
        ratio = W(n + 2) / W(n)
        assert ratio == exact_ratio_formula(n)
        assert ratio <= Fraction(2 * n + 1, 2 * n + 3) * b(n)
        assert b(n + 2) / b(n) == b_ratio_formula(n)
        balanced_denominator(n)
        algebra_checks.append(n)

    premises = [
        exact_comparison("b_16 < 1/2", b(16), Fraction(1, 2)),
        exact_comparison("b_17 < 1/10", b(17), Fraction(1, 10)),
        exact_comparison("W_23 < 1/60000", W(23), Fraction(1, 60000)),
        exact_comparison("W_24 < 1/10^10", W(24), Fraction(1, 10**10)),
    ]
    # These finite checks audit the formulas; the note proves the all-n step:
    # b_(n+2)/b_n < 4/81 for every n >= 16.
    for n in range(16, 41):
        assert b_ratio_formula(n) < Fraction(4, 81)
        assert b(n) < Fraction(1, 2)
    for n in range(23, 43):
        assert 0 < W(n) < 1

    result = {
        "status": "passed",
        "scope": "Finite exact premises plus sampled audits of identities; uniform induction is in FINITE_THRESHOLD.md.",
        "theorem": "For every n >= 23, the balanced integer polynomial Q_n=D_n P_n satisfies 0 < abs(Q_n(2*pi*i)) < 1.",
        "minimal_denominator_not_needed": True,
        "T": str(T),
        "a": str(a),
        "ratio_identity_audit_n": algebra_checks,
        "exact_premises": premises,
        "uniform_induction": {
            "all_n_at_least": 16,
            "b_two_step_ratio_strict_upper_bound": "4/81",
            "W_two_step_ratio_strict_upper_bound": "1/2",
            "threshold_bases": [23, 24],
        },
        "diagnostic_log10_W": [
            {"n": n, "log10_W": log10_fraction(W(n))}
            for n in [10, 12, 14, 16, 18, 20, 21, 22, 23, 24, 25, 26, 28, 30, 32, 40]
        ],
        "caveat": "Evaluation at an irrational imaginary argument is not an integer; this is not a rational approximation or a bad-approximability theorem.",
    }
    target = Path(__file__).with_name("finite_threshold_verification.json")
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "output": str(target)}))


if __name__ == "__main__":
    main()
