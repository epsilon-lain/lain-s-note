#!/usr/bin/env python3
"""Certify the *parameter budget* in pinned pi017, without assuming its theorem.

Only Python's standard library is needed. Exact rational base parameters and
outward-rounded rational interval powers certify every `certified` result.
No search for rational approximations to pi, interpolation rank computation,
finite-H certificate, or effective denominator bound is performed.

Run: python3 parameter_lab.py > results_parameters.json
"""
from __future__ import annotations

import argparse
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

PIN = "4eaee87548250b2682affda7d05d947c05554bae"
PRECISION = 90
SCALE = 10 ** PRECISION


def ceil_fraction(x: F) -> int:
    return -(-x.numerator // x.denominator)


def decimal(x: F, digits: int = 17, rounding: str | None = None) -> str:
    with localcontext() as ctx:
        ctx.prec = digits
        if rounding is not None:
            ctx.rounding = rounding
        return str(Decimal(x.numerator) / Decimal(x.denominator))


def upper_decimal(x: F, digits: int = 17) -> str:
    return decimal(x, digits, ROUND_CEILING)


def lower_decimal(x: F, digits: int = 17) -> str:
    return decimal(x, digits, ROUND_FLOOR)


def ln_bounds(x: F, terms: int = 128) -> tuple[F, F]:
    """Exact atanh-series bounds, valid for x>=1 (used only near 1 and at 2)."""
    if x < 1:
        lo, hi = ln_bounds(1 / x, terms)
        return -hi, -lo
    z = (x - 1) / (x + 1)
    total = F(0)
    power = z
    for k in range(terms):
        total += 2 * power / (2 * k + 1)
        power *= z * z
    # The first omitted denominator is 2*terms+1; replace every later one
    # by it and sum a positive geometric series.
    return total, total + 2 * power / ((2 * terms + 1) * (1 - z * z))


LN2_LO, LN2_HI = ln_bounds(F(2))
LAMBDA_HI = 4 * LN2_HI
C_LO = LN2_LO / 4


def unit_power_interval(x: F, exponent: int) -> tuple[F, F]:
    """Enclose x**exponent for 0<=x<=1 using bounded-size exact integers.

    Each multiplication rounds down its lower endpoint and up its upper
    endpoint. Repeated squaring requires O(log(exponent)) multiplications.
    Nothing is inferred from floating point, overflow, or underflow.
    """
    if not (0 <= x <= 1) or exponent < 0:
        raise ValueError("unit-interval base and nonnegative exponent required")
    lo = x.numerator * SCALE // x.denominator
    hi = ceil_fraction(x * SCALE)
    result_lo = result_hi = SCALE

    def multiply(a: int, b: int, upward: bool) -> int:
        product = a * b
        return (product + SCALE - 1) // SCALE if upward else product // SCALE

    while exponent:
        if exponent & 1:
            result_lo = multiply(result_lo, lo, False)
            result_hi = multiply(result_hi, hi, True)
        exponent >>= 1
        if exponent:
            lo = multiply(lo, lo, False)
            hi = multiply(hi, hi, True)
    return F(result_lo, SCALE), F(result_hi, SCALE)


def bases(nu: F) -> dict:
    if nu <= 2:
        return {"nu": nu, "status": "structurally_infeasible_in_this_budget"}
    b = (F(1, 2) + 1 - 1 / nu) / 2
    delta = min(F(1, 4), (2 * b - 1) / (2 * b * b))
    theta, A = 1 - delta, 1 - b * delta
    C = (A / theta + 1 / A) / 2
    B = (A + min(1 / C, C * theta)) / 2
    g0 = nu * (A - theta) - (1 - theta)
    eta = g0 / (2 * nu * A)
    g = nu * (A * (1 - eta) - theta) - (1 - theta)
    eps0 = min(g / 2, F(1, 2))
    F0 = max(ceil_fraction(2 / theta), ceil_fraction(3 * nu / eps0)) + 1
    assert 0 < theta < A < B < 1
    assert A * A < theta and 1 < C and B < 1 / C and B < C * theta < 1
    assert 0 < eta < 1 and g == g0 / 2 and 0 < g
    assert F0 > 2 / theta and nu / F0 < eps0 / 3
    return dict(nu=nu, b=b, delta=delta, theta=theta, A=A, B=B, C=C,
                g0=g0, eta=eta, g=g, eps0=eps0, F0=F0,
                status="base_parameters_exact")


def at_dimension(p: dict, m: int) -> dict:
    B, C, theta, A = (p[k] for k in ("B", "C", "theta", "A"))
    # floor(x)>=x/2 for every x>=1: on [1,2) this is immediate, and on
    # [2,infinity) use floor(x)>=x-1>=x/2. Here C>1 and m>=1.
    floor_comparison_valid = C > 1 and m >= 1
    inverse_growth = B / (C * theta)
    error_decay = C * B
    inverse_collision = A / B
    v_inverse_lo, v_inverse_hi = unit_power_interval(inverse_growth, m)
    error_decay_lo, error_decay_hi = unit_power_interval(error_decay, m)
    collision_inverse_lo, collision_inverse_hi = unit_power_interval(inverse_collision, m)
    dimension_error_hi = ((LAMBDA_HI * p["F0"] * m + 2 * LN2_HI)
                          * v_inverse_hi + 100 * error_decay_hi)
    dimension_error_lo = ((LAMBDA_HI * p["F0"] * m + 2 * LN2_HI)
                          * v_inverse_lo + 100 * error_decay_lo)
    # cL >= c_lo*eta**2/[2(m+1)*(A/B)**m]. Strictly positive lower
    # bound even when the exact power is too small to display directly.
    collision_lo = C_LO * p["eta"] ** 2 / (2 * (m + 1) * collision_inverse_hi)
    collision_hi = (C_LO * p["eta"] ** 2 / (2 * (m + 1) * collision_inverse_lo)
                    if collision_inverse_lo > 0 else None)
    dimension_ok = (floor_comparison_valid and
                    dimension_error_hi < p["eps0"] / 3 and collision_lo > 2)
    return dict(m=m, dimension_error_upper=dimension_error_hi,
                dimension_error_lower=dimension_error_lo,
                collision_cL_lower=collision_lo,
                collision_cL_upper=collision_hi,
                dimension_error_ok=dimension_error_hi < p["eps0"] / 3,
                collision_ok=collision_lo > 2,
                conservative_dimension_failure_certified=dimension_error_lo >= p["eps0"] / 3,
                conservative_collision_failure_certified=(collision_hi is not None and collision_hi <= 2),
                floor_comparison_valid=floor_comparison_valid,
                dimension_ok=dimension_ok,
                inverse_collision_power_lower=collision_inverse_lo,
                inverse_collision_power_upper=collision_inverse_hi)


def choose_dimension(p: dict, max_dimension: int) -> dict:
    # From here (a*m+b)*p**m decreases and r**m/(m+1) increases. This
    # permits binary search without assuming global monotonicity near m=1.
    decay1 = p["B"] / (p["C"] * p["theta"])
    decay2 = p["A"] / p["B"]
    start = max(1, ceil_fraction(2 / (1 - decay1)),
                ceil_fraction(2 / (1 - decay2)))
    upper = start
    while upper <= max_dimension and not at_dimension(p, upper)["dimension_ok"]:
        upper *= 2
    if upper > max_dimension:
        return {"status": "search_limit_reached_without_certificate",
                "max_dimension": max_dimension, "first_monotone_dimension": start}
    lower = start - 1
    while upper - lower > 1:
        mid = (lower + upper) // 2
        if at_dimension(p, mid)["dimension_ok"]:
            upper = mid
        else:
            lower = mid
    result = at_dimension(p, upper)
    result["first_monotone_dimension"] = start
    return result


def certificate(nu: F, max_dimension: int) -> dict:
    p = bases(nu)
    if nu <= 2:
        return {"nu": str(nu), "status": p["status"],
                "reason": "nu=2: 2A-theta-1 < -(1-A)^2 < 0 when A^2<theta; nu<2 is even worse"}
    result = choose_dimension(p, max_dimension)
    if result.get("status"):
        return {"nu": str(nu), **result}
    m = result["m"]
    # log(400)<6 and log(C)<=C-1. Thus this rational U dominates the
    # entire weight-error numerator whenever all w_i>=W.
    U = LAMBDA_HI * m + p["theta"] + p["nu"] + 6 + 2 * m * (p["C"] - 1)
    W = ceil_fraction(3 * U / p["eps0"]) + 1
    weight_error_hi = U / W
    total_error_hi = (p["nu"] / p["F0"] + result["dimension_error_upper"]
                      + weight_error_hi)
    margin_lo = min(result["collision_cL_lower"] - 1, p["g"]) - total_error_hi
    normalized_hi = max(total_error_hi / p["g"],
                        (1 + total_error_hi) / result["collision_cL_lower"])
    sigma = min(F(1, 12 * (m + 1)),
                (1 - p["C"] * p["theta"]) / (12 * p["C"] * p["theta"]),
                (1 - p["theta"]) / (12 * p["theta"]))
    assert weight_error_hi < p["eps0"] / 3
    assert total_error_hi < p["eps0"] < p["g"]
    assert margin_lo > 0 and normalized_hi < 1
    ablation_dimensions = sorted(set([1, 10, 100, 1000, m // 10, m // 2,
                                     max(1, m - 1), m]))
    ablations = []
    for test_m in ablation_dimensions:
        test = at_dimension(p, test_m)
        ablations.append({"m": test_m,
                          "conservative_dimension_error_lower": lower_decimal(test["dimension_error_lower"]),
                          "dimension_error_upper": upper_decimal(test["dimension_error_upper"]),
                          "dimension_error_target": decimal(p["eps0"] / 3),
                          "collision_cL_lower": lower_decimal(test["collision_cL_lower"]),
                          "dimension_error_ok": test["dimension_error_ok"],
                          "collision_ok": test["collision_ok"],
                          "conservative_dimension_failure_certified": test["conservative_dimension_failure_certified"],
                          "conservative_collision_failure_certified": test["conservative_collision_failure_certified"],
                          "floor_comparison_valid": test["floor_comparison_valid"],
                          "both_certified": test["dimension_ok"]})
    lnC_lo, lnC_hi = ln_bounds(p["C"])
    return {
        "nu": str(nu),
        "status": "certified_parameter_budget_with_conditional_threshold_policy",
        "base_parameters_exact": {k: str(v) for k, v in p.items() if isinstance(v, F)},
        "F0_integer": p["F0"], "dimension_m": m,
        "first_monotone_dimension": result["first_monotone_dimension"],
        "exact_large_data": {"K": "floor(C**m)", "w0": "B**(-m)",
                             "v0": "2*K*(theta/B)**m", "row_column_volume_ratio": "1/2"},
        "natural_log_C_power_interval": [lower_decimal(m * lnC_lo), upper_decimal(m * lnC_hi)],
        "common_minimum_weight_W": W,
        "sigma_exact": str(sigma),
        "certified_bounds": {
            "nu_over_F0_upper": upper_decimal(p["nu"] / p["F0"]),
            "dimension_error_upper": upper_decimal(result["dimension_error_upper"]),
            "dimension_error_target": decimal(p["eps0"] / 3),
            "weight_error_upper": upper_decimal(weight_error_hi),
            "total_error_upper": upper_decimal(total_error_hi),
            "epsilon0": decimal(p["eps0"]), "gap_g": decimal(p["g"]),
            "collision_cL_lower": lower_decimal(result["collision_cL_lower"]),
            "signed_margin_Phi_lower": lower_decimal(margin_lo),
            "normalized_budget_upper": upper_decimal(normalized_hi)},
        "threshold_policy": {
            "C_star": "2*m!*(((m+2)/sigma)**m)",
            "M_previous": "max(1,v0,w_(i-1)/theta), omitting last entry at i=1",
            "T_i": "max(W, 2*C_star*theta**(-m)*max(1,v0/w0)*(M_previous**m))",
            "choice": "Under the contradictory hypothesis choose a good q_i with ceil(log(q_i))>T_i, then set w_i=ceil(log(q_i)); choose p_i accordingly. No actual q_i is supplied."},
        "ablation_m": ablations,
        "not_certified": ["existence of any actual good pi approximants",
                          "truth of separated-weight interpolation",
                          "a finite H or a nonzero finite determinant",
                          "an effective Q(nu)", "bad approximability of pi"]}


def self_checks() -> dict:
    checks = 0
    for x in [F(0), F(1), F(1, 3), F(99, 100), F(123, 456)]:
        for n in [0, 1, 2, 3, 9, 50, 121]:
            lo, hi = unit_power_interval(x, n)
            assert lo <= x ** n <= hi
            checks += 1
    assert LN2_LO < F(7, 10) and LN2_HI > F(69, 100)
    # Exact rational endpoints of a positive atanh tail certify these.
    assert F(69, 100) < LN2_LO < LN2_HI < F(7, 10)
    return {"exact_power_enclosure_crosschecks": checks,
            "ln2_interval": [lower_decimal(LN2_LO, 40), upper_decimal(LN2_HI, 40)],
            "precision_decimal_places": PRECISION}


def source_audit() -> list:
    from source_support import audit_source_files
    paths = ["build/main.tex", "build/sections/determinant.tex",
             "build/sections/conclusion.tex", "build/sections/series.tex",
             "build/sections/interpolation.tex"]
    result=audit_source_files(paths)
    return result['records'] if result['status']=='passed' else [result]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nu", nargs="*", default=["3", "2.5", "2.1", "2.01", "2.001", "2", "1.99"])
    parser.add_argument("--max-dimension", type=int, default=10 ** 12)
    args = parser.parse_args()
    output = {"scope": "Exact parameter-budget audit only; conditional on interpolation and hypothetical unbounded good denominators",
              "pinned_upstream_commit": PIN,
              "numeric_display": "Exact rational comparisons precede display; upper/lower bound strings are rounded outward. Other decimal strings are diagnostics.",
              "checks": self_checks(),
              "source_byte_audit": source_audit(),
              "certificates": [certificate(F(nu), args.max_dimension) for nu in args.nu],
              "endpoint_model": {
                  "kappa_definition": "q**2 * abs(pi-p/q)",
                  "extra_exponent_definition": "s=log(1/kappa)/log(q)",
                  "fixed_nu_hypothesis": "kappa<=q**(-(nu-2)), so s>=nu-2>0",
                  "slow_decay_example": "kappa=1/log(q) tends to zero but s=log(log(q))/log(q) tends to zero",
                  "obstruction": "A**2<theta and nu=2 imply g0<-(1-A)**2<0 even before eta and error costs"}}
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
