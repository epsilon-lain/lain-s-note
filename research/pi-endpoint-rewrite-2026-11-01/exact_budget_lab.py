#!/usr/bin/env python3
"""Exact checks for adaptive_budget.md; finite examples are not theorem proofs."""
from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import sympy as sp


def weighted_tuples(weights: tuple[int, ...], budget: int):
    if not weights:
        yield ()
        return
    first, *rest = weights
    for value in range(budget // first + 1):
        for tail in weighted_tuples(tuple(rest), budget - first * value):
            yield (value,) + tail


def as_rat(value: Fraction):
    return {"numerator": value.numerator, "denominator": value.denominator}


def main():
    nu, A, chi = sp.symbols("nu A chi", positive=True)
    r = nu**2 / (4 * (nu - 1))
    difference = r / chi - (nu * A - chi) / ((nu - 1) * A**2)
    target = (nu * A - 2 * chi) ** 2 / (4 * chi * (nu - 1) * A**2)
    assert sp.cancel(difference - target) == 0

    rho, x = sp.symbols("rho x", positive=True)
    mean_checks = []
    for m in range(1, 9):
        x0 = (rho / (m + 1)) ** sp.Rational(1, m)
        integral = x0 - sp.Rational(m + 1, m + 1) * x0 ** (m + 1) / rho
        expected = sp.Rational(m, m + 1) * x0
        # x0^(m+1)/rho = x0/(m+1), since rho is positive.
        assert sp.simplify(sp.powdenest(integral - expected, force=True)) == 0
        mean_checks.append(m)

    # These exact checks are illustrative samples; the global inequalities
    # are proved in the accompanying note by ln(1+x)<x and e<3.
    exponential_samples = []
    for m in range(1, 81):
        value = Fraction(m + 2, m) ** m
        assert value < 9
        exponential_samples.append({"m": m, "less_than_9": True})

    row_examples = []
    for m in (1, 2, 3):
        # theta=1/2, v0=1, w_i=1. Strict row inequality s+2Σβ<H.
        limit = Fraction(m, 2 * (m + 2))
        for H in (10, 20, 40):
            rows = list(weighted_tuples((1,) + (2,) * m, H - 1))
            b = Fraction(sum(sum(row[1:]) for row in rows), len(rows) * H)
            row_examples.append({"m": m, "H": H, "row_count": len(rows),
                                 "b_H": as_rat(b), "limit": as_rat(limit)})

    selected_columns = []
    # All columns h+Σα≤H, select floor(N/2) columns of least Σα.
    # The asymptotic lower bound is displayed, not asserted as a finite-H bound.
    for m in (1, 2, 3):
        for H in (10, 20, 30):
            columns = list(weighted_tuples((1,) * (m + 1), H))
            costs = sorted(sum(col[1:]) for col in columns)
            count = len(costs) // 2
            mean = Fraction(sum(costs[:count]), count * H)
            selected_columns.append({"m": m, "H": H, "all_columns": len(columns),
                                     "selected": count, "chi_H": as_rat(mean),
                                     "asymptotic_rho": "1/2",
                                     "asymptotic_lower_bound_formula":
                                     f"{m}/{m+1} * (1/{2*(m+1)})^(1/{m})"})

    dimension_bounds = []
    for value in (Fraction(3), Fraction(5, 2), Fraction(21, 10),
                  Fraction(201, 100), Fraction(2001, 1000)):
        s = value - 2
        bound = 24 * (1 + s) / s**2
        dimension_bounds.append({"nu": as_rat(value),
                                 "necessary_strict_m_lower_bound": as_rat(bound),
                                 "interpretation": "necessary only, not feasible dimension"})

    certificates = []
    value = Fraction(21, 10)
    R = value**2 / (4 * (value - 1))
    for m in (100, 1000, 3000, 10000):
        lhs = R**m
        # eta<1, so this relaxed screen must hold if the stronger necessary
        # condition with eta^(-4) holds. Passing is not a proof of feasibility.
        rhs = Fraction(12800 * (m + 1), 27)
        certificates.append({"nu": as_rat(value), "m": m,
                             "passes_relaxed_necessary_screen": lhs > rhs,
                             "comparison_arithmetic": "exact Fraction",
                             "R": as_rat(R)})

    result = {
        "status": "passed",
        "scope": "exact algebra and finite illustrative checks; theorem proof in adaptive_budget.md",
        "square_identity": True,
        "distribution_mean_symbolic_m": mean_checks,
        "one_plus_2_over_m_power_samples": exponential_samples,
        "row_mean_finite_examples": row_examples,
        "least_cost_column_finite_examples": selected_columns,
        "necessary_dimension_bounds": dimension_bounds,
        "necessary_screen_examples": certificates,
        "slow_kappa_model": {
            "hypotheses": ["h_i=ln(1/kappa_i)<=ln(w_i)",
                           "positive integer w_i strictly increasing", "s=min_i h_i/w_i"],
            "analytic_conclusion": "m*s^2<1, contradicting necessary m*s^2>24*(1+s)",
            "not_a_claim_about_pi": True,
            "pooled_exponent_only": True,
        },
    }
    target_path = Path(__file__).with_name("results_adaptive_budget.json")
    target_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "result": str(target_path),
                      "symbolic_identities": 9, "finite_examples": len(row_examples) + len(selected_columns)},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
