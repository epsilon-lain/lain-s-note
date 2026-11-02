#!/usr/bin/env python3
"""Independently rebuild the saved fractions and check their exact certificates.

This script intentionally does not import run_experiment.  Its integer spinor
evaluation and difference-form Aitken reconstruction differ from the generator.
"""

from fractions import Fraction as Q
from math import gcd
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parent


def spinor_value(t, c):
    a, b = c.numerator, c.denominator
    p = b
    q = (a + b) * t - a * t * t
    dq = a + b - 2 * a * t
    return Q(4 * p) * dq / (p * p + q * q)


def trap_by_intervals(n, c):
    # Sum trapezoidal areas independently interval by interval.
    area = Q(0)
    for j in range(n):
        left, right = Q(j, n), Q(j + 1, n)
        area += (right - left) * (spinor_value(left, c) + spinor_value(right, c)) / 2
    return area


def simpson_by_panels(n, c):
    area = Q(0)
    for j in range(0, n, 2):
        left, middle, right = Q(j, n), Q(j + 1, n), Q(j + 2, n)
        area += (right - left) * (spinor_value(left, c)
                                    + 4 * spinor_value(middle, c)
                                    + spinor_value(right, c)) / 6
    return area


def reconstruct(method, n, c):
    if method == "trapezoid":
        return trap_by_intervals(n, c)
    if method == "simpson":
        return simpson_by_panels(n, c)
    assert method == "aitken"
    s0, s1, s2 = [trap_by_intervals(n * k, c) for k in (1, 2, 4)]
    denominator = (s2 - s1) - (s1 - s0)
    assert denominator != 0
    return s0 - (s1 - s0) ** 2 / denominator


def atan_bounds(k, terms):
    # Build powers recursively, independently of generator's summand formula.
    term_power, partial = Q(1, k), Q(0)
    for j in range(terms):
        partial += term_power / (2 * j + 1)
        term_power /= -(k * k)
    other = partial + term_power / (2 * terms + 1)
    return sorted([partial, other])


def main():
    saved = json.loads((ROOT / "results.json").read_text())
    a, b = atan_bounds(5, 70), atan_bounds(239, 20)
    lo, hi = 16 * a[0] - 4 * b[1], 16 * a[1] - 4 * b[0]
    assert list(map(Q, saved["pi_certificate"]["interval"])) == [lo, hi]
    assert hi - lo < Q(1, 10 ** 95)
    rows = {}
    for row in saved["candidates"]:
        method, n, c = row["method"], row["base_n"], Q(row["c"])
        r = reconstruct(method, n, c)
        assert r == Q(row["fraction"])
        p, q = int(row["p"]), int(row["q"])
        assert r == Q(p, q) and gcd(abs(p), q) == 1 and q > 0
        assert r.numerator == p and r.denominator == q
        rawp, rawq = int(row["cost"]["raw_numerator"]), int(row["cost"]["raw_denominator"])
        g = gcd(abs(rawp), rawq)
        assert g == int(row["cost"]["actual_gcd"])
        assert rawp // g == p and rawq // g == q
        assert abs(p).bit_length() == row["p_bits"] and q.bit_length() == row["q_bits"]
        assert row["integer_spinor"] == {
            "P": [c.denominator], "Q": [0, c.numerator + c.denominator, -c.numerator]}
        assert row["integer_spinor_height"] == max(c.denominator, abs(c.numerator + c.denominator), abs(c.numerator))
        err_lo = Q(0) if lo <= r <= hi else min(abs(lo - r), abs(hi - r))
        err_hi = max(abs(lo - r), abs(hi - r))
        assert list(map(Q, row["absolute_error_interval"])) == [err_lo, err_hi]
        assert list(map(Q, row["quality_interval"])) == [q * q * err_lo, q * q * err_hi]
        max_grid = (4 if method == "aitken" else 1) * n
        assert row["unique_rational_nodes"] == max_grid + 1
        assert row["max_node_denominator"] == max_grid
        rows[(method, n, str(c))] = row
    assert len(rows) == saved["candidate_count"] == 60
    for cmp in saved["comparisons"]:
        old = rows[tuple(cmp["baseline"])]
        new = rows[tuple(cmp["candidate"])]
        assert old["unique_rational_nodes"] == new["unique_rational_nodes"]
        old_lo, old_hi = map(Q, old["quality_interval"])
        new_lo, new_hi = map(Q, new["quality_interval"])
        assert list(map(Q, cmp["quality_improvement_factor_interval"])) == [old_lo / new_hi, old_hi / new_lo]
    # Check the summary's concrete inputs and output exactly.
    assert trap_by_intervals(1, Q(0)) == 3
    assert trap_by_intervals(2, Q(0)) == Q(31, 10)
    assert trap_by_intervals(4, Q(0)) == Q(5323, 1700)
    assert reconstruct("aitken", 1, Q(0)) == Q(368, 117)
    for n in (4, 8, 16):
        assert Q(rows[("simpson", n, "0")]["quality_interval"][1]) < Q(rows[("trapezoid", n, "0")]["quality_interval"][0])
    assert Q(rows[("aitken", 4, "0")]["quality_interval"][0]) > Q(rows[("trapezoid", 16, "0")]["quality_interval"][1])
    ref = saved["posthoc_reference_baselines"]["355/113"]
    ref_hi = Q(ref["quality_interval"][1])
    assert all(Q(row["quality_interval"][0]) > ref_hi for row in rows.values())
    # Uniform integer-fee identity, with rational alpha replacing pi for an
    # independent exact algebra check.  The NOTE proves it for all real alpha.
    identity_cases = 0
    for aa, bb, cc, dd in [(3, 4, 9, 7), (4, 9, 3, 5), (-2, 3, 11, 13)]:
        delta, numerator = aa - 2 * bb + cc, aa * cc - bb * bb
        assert delta != 0
        gg = gcd(abs(numerator), abs(dd * delta))
        rr = Q(numerator, dd * delta)
        for alpha in [Q(-3, 2), Q(0), Q(22, 7)]:
            residuals = [xx - dd * alpha for xx in [aa, bb, cc]]
            mismatch = residuals[0] * residuals[2] - residuals[1] ** 2
            assert rr.denominator ** 2 * abs(rr - alpha) == Q(dd * abs(delta), gg * gg) * abs(mismatch)
            identity_cases += 1
    report = {
        "status": "PASS", "independent_candidate_reconstructions": len(rows),
        "actual_gcd_checks": len(rows), "exact_error_and_quality_certificates": len(rows),
        "same_node_comparisons": len(saved["comparisons"]),
        "all_candidates_worse_quality_than_355_over_113": True,
        "uniform_integer_fee_identity_cases": identity_cases,
        "scope": "Exact finite verification; general identities are proved in NOTE.md.",
    }
    (ROOT / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
