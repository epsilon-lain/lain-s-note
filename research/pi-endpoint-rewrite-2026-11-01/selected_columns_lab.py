#!/usr/bin/env python3
"""Exact rank/weight experiments; these finite matrices do not certify pi017."""
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = spec_from_file_location(
    "pi017_interpolation_lab",
    HERE.parent / "pi-construction-mechanisms-2026-11-01/interpolation_lab.py",
)
LAB = module_from_spec(SPEC)
SPEC.loader.exec_module(LAB)


def cdf(m, x):
    return (m + 1) * x**m - m * x**(m + 1)


def tail_moment(m, x):
    return (Fraction(m, m + 2) - m * x**(m + 1)
            + Fraction(m * (m + 1), m + 2) * x**(m + 2))


def upper_fraction_cost(m, rho, iterations=100):
    lo, hi = Fraction(0), Fraction(1)
    for _ in range(iterations):
        mid = (lo + hi) / 2
        if 1 - cdf(m, mid) > rho:
            lo = mid
        else:
            hi = mid
    # The integrand is nonnegative; tail_moment decreases with x.
    return tail_moment(m, hi) / rho, tail_moment(m, lo) / rho


def finite_rank_profile(H):
    matrix, rows, columns = LAB.jet_matrix(
        H, 1, 1, Fraction(1, 2), 4, [0])
    assert matrix.rank() == matrix.rows
    order = sorted(range(matrix.cols), key=lambda j: (columns[j][1], columns[j][0]))
    sorted_matrix = matrix[:, order]
    pivots = sorted_matrix.rref()[1]
    chosen = [columns[order[j]] for j in pivots]
    M = matrix.rows
    col_cost = sum(a for h, a in chosen)
    row_cost = sum(b for j, s, b in rows)
    profile = []
    for a in range(H + 1):
        indices = [j for j, (h, exponent) in enumerate(columns) if exponent <= a]
        rank = matrix[:, indices].rank()
        profile.append(rank)
    # Discrete layer-cake identity for integer costs.
    assert col_cost == sum(M - profile[a] for a in range(H))
    assert col_cost >= row_cost
    assert len(chosen) == M
    return {
        "H": H, "rows": M, "columns": matrix.cols,
        "rank": matrix.rank(), "minimum_cost_basis": chosen,
        "rank_at_transverse_cutoffs_0_through_H": profile,
        "minimum_basis_transverse_sum": col_cost,
        "row_transverse_sum": row_cost,
        "c_col": str(Fraction(col_cost, M * H)),
        "b_bar": str(Fraction(row_cost, M * H)),
        "c_col_minus_b_bar": str(Fraction(col_cost - row_cost, M * H)),
        "status": "exact finite toy matrix, not an asymptotic source certificate",
    }


def main():
    quantiles = []
    for m in (1, 2, 5, 10, 50, 100):
        lo, hi = upper_fraction_cost(m, Fraction(1, 2))
        assert hi - lo < Fraction(1, 10**25)
        quantiles.append({
            "m": m, "rho": "1/2",
            "upper_half_mean_lower_exact": str(lo),
            "upper_half_mean_upper_exact": str(hi),
            "upper_half_mean_decimal_for_display": float((lo + hi) / 2),
            "greedy_rank_profile_asymptotic_bound": str(Fraction(m + 1, m + 2)),
            "all_columns_mean": str(Fraction(m, m + 2)),
        })
    result = {
        "finite_rank_profiles": [finite_rank_profile(H) for H in range(2, 8)],
        "full_simplex_transverse_distribution": "Beta(m,2)",
        "quantile_bounds": quantiles,
        "checks": ["exact rank", "minimum-cost basis by sorted pivot selection",
                   "rank profile layer-cake identity", "c_col >= b_bar",
                   "rational quantile enclosures"],
    }
    target = HERE / "selected_columns_results.json"
    target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(target), "finite_models": 6,
                      "quantile_enclosures": len(quantiles)}, indent=2))


if __name__ == "__main__":
    main()
