#!/usr/bin/env python3
"""Exact finite certificate for the user's root family, up to q <= 10**100.

All mathematical comparisons use integers and Fraction. Display decimals are
explicitly secondary. The sole number-theoretic theorem used to pass from the
finite list of convergents to *all* reduced fractions is Legendre's criterion.
This is a finite certificate; it makes no assertion for q > 10**100.
"""

from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent


def pack(r):
    return {"numerator": str(r.numerator), "denominator": str(r.denominator)}


def unpack(v):
    return F(int(v["numerator"]), int(v["denominator"]))


def atan_bounds(k, terms):
    # The next term brackets the error of this alternating decreasing series.
    total = F(0)
    for j in range(terms):
        total += F((-1) ** j, (2 * j + 1) * k ** (2 * j + 1))
    next_term = F((-1) ** terms, (2 * terms + 1) * k ** (2 * terms + 1))
    return min(total, total + next_term), max(total, total + next_term)


def quality_bounds(p, q, lower, upper):
    value = F(p, q)
    if value < lower:
        error_lower, error_upper = lower - value, upper - value
        side = "below_pi"
    elif value > upper:
        error_lower, error_upper = value - upper, value - lower
        side = "above_pi"
    else:
        raise AssertionError("pi enclosure insufficient to separate candidate")
    return q * q * error_lower, q * q * error_upper, side


def construct():
    a_lower, a_upper = atan_bounds(5, 200)
    b_lower, b_upper = atan_bounds(239, 70)
    pi_lower = 16 * a_lower - 4 * b_upper
    pi_upper = 16 * a_upper - 4 * b_lower
    assert F(3) < pi_lower < pi_upper < F(4)
    assert pi_upper - pi_lower < F(1, 10 ** 281)

    bound = 10 ** 100
    left, right = pi_lower, pi_upper
    p_older, p_old, q_older, q_old = 0, 1, 1, 0
    convergents, coefficients, records = [], [], []
    running_best = None
    while True:
        digit_lower = left.numerator // left.denominator
        digit_upper = right.numerator // right.denominator
        assert digit_lower == digit_upper, "pi enclosure too wide for next digit"
        digit = digit_lower
        index = len(coefficients)
        coefficients.append(digit)
        p = digit * p_old + p_older
        q = digit * q_old + q_older
        assert gcd(p, q) == 1 and q > 0
        assert 3 * q <= p <= 4 * q
        row = {
            "index": index,
            "common_cf_coefficient": digit,
            "floor_lower_equals_floor_upper": True,
            "p": str(p),
            "q": str(q),
            "root_index_p_minus_3q": str(p - 3 * q),
            "reduced": True,
            "determinant_with_previous": str(p * q_old - p_old * q),
        }
        if q > bound:
            first_beyond = row
            break
        quality_lower, quality_upper, side = quality_bounds(p, q, pi_lower, pi_upper)
        row.update({
            "side": side,
            "quality_lower": pack(quality_lower),
            "quality_upper": pack(quality_upper),
            "quality_midpoint_display_only": float((quality_lower + quality_upper) / 2),
        })
        if running_best is None or quality_upper < running_best[0]:
            row["new_record"] = True
            records.append({"index": index, "p": str(p), "q": str(q),
                            "quality_midpoint_display_only": row["quality_midpoint_display_only"]})
            running_best = (quality_lower, quality_upper, p, q)
        else:
            assert quality_lower > running_best[1]
            row["new_record"] = False
        convergents.append(row)
        p_older, p_old, q_older, q_old = p_old, p, q_old, q
        assert left > digit and right > digit
        # Reciprocal reverses order; a common digit makes this exact interval
        # propagation valid and certifies every retained partial quotient.
        left, right = 1 / (right - digit), 1 / (left - digit)

    assert running_best[2:] == (355, 113)
    best_lower, best_upper = running_best[:2]
    assert best_upper < F(1, 2)
    for row in convergents:
        is_best = (row["p"], row["q"]) == ("355", "113")
        row["quality_strictly_above_355_113"] = not is_best
        if not is_best:
            assert unpack(row["quality_lower"]) > best_upper

    return {
        "description": "Finite all-reduced-fraction minimum of q^2*abs(pi-p/q)",
        "scope": {"positive_reduced_denominator_at_most": str(bound),
                  "numerator_range": "3*q <= p <= 4*q", "infinite_conclusion": False},
        "pi_enclosure": {
            "method": "Machin: pi = 16*atan(1/5)-4*atan(1/239), alternating next-term bounds",
            "atan_5_terms": 200, "atan_239_terms": 70,
            "atan_5_lower": pack(a_lower), "atan_5_upper": pack(a_upper),
            "atan_239_lower": pack(b_lower), "atan_239_upper": pack(b_upper),
            "lower": pack(pi_lower), "upper": pack(pi_upper),
            "width": pack(pi_upper - pi_lower),
            "width_less_than": pack(F(1, 10 ** 281)),
        },
        "cf_certification": {
            "rule": "a=floor(L)=floor(U); next interval is [1/(U-a),1/(L-a)]",
            "initial_interval": "pi_enclosure.lower, pi_enclosure.upper",
            "common_coefficients_including_first_beyond": coefficients,
            "convergent_recurrence": "p=a*p_old+p_older, q=a*q_old+q_older; seeds p=(0,1), q=(1,0)",
            "correct_seeds": {"p_older": 0, "p_old": 1, "q_older": 1, "q_old": 0},
            "first_denominator_beyond_bound_certified": True,
        },
        "convergent_count_within_bound": len(convergents),
        "convergents": convergents,
        "first_convergent_beyond_bound": first_beyond,
        "record_sequence": records,
        "minimum": {"p": "355", "q": "113", "quality_lower": pack(best_lower),
                    "quality_upper": pack(best_upper), "strictly_below_one_half": True,
                    "all_other_listed_convergents_have_disjoint_larger_quality_intervals": True},
        "all_fraction_coverage": {
            "external_theorem": "Legendre: for coprime integers p,q with q>0, abs(alpha-p/q)<1/(2*q*q) implies p/q is a regular continued-fraction convergent of irrational alpha.",
            "application": "Any fraction improving 355/113 would have quality < its upper bound < 1/2; Legendre puts it in the certified list. Every other item has a strictly larger quality interval. Subsequent convergent denominators exceed the denominator cap.",
            "root_family_link": "For every listed fraction, n=p-3*q is an integer in [0,q] and F_q(p/q)=0. Coprimality is certified by adjacent determinant +/-1 (and checked directly).",
            "no_large_scan": "No denominator-by-denominator scan was performed. Exactly 194 convergents were scored.",
        },
    }


def replay_saved(data):
    enclosure = data["pi_enclosure"]
    lo, hi = unpack(enclosure["lower"]), unpack(enclosure["upper"])
    a, b = atan_bounds(5, enclosure["atan_5_terms"])
    c, d = atan_bounds(239, enclosure["atan_239_terms"])
    assert lo == 16 * a - 4 * d and hi == 16 * b - 4 * c
    L, U = lo, hi
    p0, p1, q0, q1 = 0, 1, 1, 0
    best_hi = unpack(data["minimum"]["quality_upper"])
    rows = data["convergents"] + [data["first_convergent_beyond_bound"]]
    for i, row in enumerate(rows):
        aL, aU = L.numerator // L.denominator, U.numerator // U.denominator
        assert aL == aU == row["common_cf_coefficient"]
        p, q = aL * p1 + p0, aL * q1 + q0
        assert (p, q) == (int(row["p"]), int(row["q"]))
        assert gcd(p, q) == 1 and 0 <= p - 3 * q <= q
        assert abs(p * q1 - p1 * q) == 1
        if i < len(data["convergents"]):
            lower, upper, side = quality_bounds(p, q, lo, hi)
            assert lower == unpack(row["quality_lower"])
            assert upper == unpack(row["quality_upper"]) and side == row["side"]
            assert q <= 10 ** 100
            if (p, q) != (355, 113):
                assert lower > best_hi
        else:
            assert q > 10 ** 100
        p0, p1, q0, q1 = p1, p, q1, q
        L, U = 1 / (U - aL), 1 / (L - aL)
    return {"all_exact_checks_passed": True, "convergents_scored": len(data["convergents"]),
            "first_beyond_replayed": True, "all_reduced": True,
            "all_other_quality_intervals_strictly_larger": True}


if __name__ == "__main__":
    ROOT.mkdir(parents=True, exist_ok=True)
    result = construct()
    target = ROOT / "cf_all_denominators_1e100.json"
    target.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    saved = json.loads(target.read_text())
    verification = replay_saved(saved)
    (ROOT / "cf_verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    print(json.dumps({"result_path": str(target), "result_bytes": target.stat().st_size,
                      "verification": verification, "record_sequence": result["record_sequence"],
                      "first_denominator_beyond": result["first_convergent_beyond_bound"]["q"]},
                     ensure_ascii=False, indent=2))
