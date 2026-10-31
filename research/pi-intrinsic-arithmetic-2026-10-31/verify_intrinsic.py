#!/usr/bin/env python3
"""Certified finite checks for the intrinsic pi-digit reciprocal construction.

No continued-fraction data or floating-point arithmetic is used. A rational
Machin enclosure certifies greedy pi-digits of 1. Integer polynomial arithmetic
checks D*C=1 and R*H=1 through a finite order. Ratio certificates retain gcd and
denominator costs. Finite checks do not prove or disprove P2848.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path


def ceil_div(numerator, denominator):
    return -((-numerator) // denominator)


def arctan_interval(inverse, target_denominator):
    """Alternating rational series, enclosed by two successive partial sums."""
    partial = Fraction(0)
    power = inverse
    index = 0
    while True:
        sign = 1 if index % 2 == 0 else -1
        partial += Fraction(sign, (2 * index + 1) * power)
        next_power = power * inverse * inverse
        next_denominator = (2 * index + 3) * next_power
        if next_denominator >= target_denominator:
            adjacent = partial + Fraction(-sign, next_denominator)
            return min(partial, adjacent), max(partial, adjacent), index + 1
        power = next_power
        index += 1


def pi_fixed_interval(bits):
    scale = 1 << bits
    target = scale << 40
    a_lo, a_hi, a_terms = arctan_interval(5, target)
    b_lo, b_hi, b_terms = arctan_interval(239, target)
    lower = 16 * a_lo - 4 * b_hi
    upper = 16 * a_hi - 4 * b_lo
    low = lower.numerator * scale // lower.denominator
    high = ceil_div(upper.numerator * scale, upper.denominator)
    assert 3 * scale < low < high < 4 * scale
    return low, high, scale, {
        "formula": "pi = 16 atan(1/5) - 4 atan(1/239)",
        "atan_1_over_5_terms": a_terms,
        "atan_1_over_239_terms": b_terms,
        "scaled_interval_width": high - low,
    }


def scientific_fraction(numerator, denominator, upper=False, significant=24):
    """Outward decimal display of an exact nonnegative rational number."""
    if numerator == 0:
        return "0"
    assert numerator > 0
    exponent = len(str(numerator)) - len(str(denominator))
    if exponent >= 0:
        if numerator < denominator * 10 ** exponent:
            exponent -= 1
    else:
        if numerator * 10 ** (-exponent) < denominator:
            exponent -= 1
    shift = significant - 1 - exponent
    if shift >= 0:
        mantissa, remainder = divmod(numerator * 10 ** shift, denominator)
    else:
        mantissa, remainder = divmod(numerator, denominator * 10 ** (-shift))
    if upper and remainder:
        mantissa += 1
    mantissa_text = str(mantissa)
    if len(mantissa_text) > significant:
        mantissa_text = mantissa_text[:-1]
        exponent += 1
    return mantissa_text[0] + "." + mantissa_text[1:] + "e" + str(exponent)


def display_interval(lower, upper, scale):
    return [scientific_fraction(lower, scale),
            scientific_fraction(upper, scale, upper=True)]


def greedy_digits(count, pi_lower, pi_upper, scale):
    lower = upper = scale  # r_0 = 1, then all r_n are in [0,1).
    digits = []
    remainders = [(lower, upper)]
    least_margin = None
    for index in range(1, count + 1):
        product_lower = pi_lower * lower
        product_upper = pi_upper * upper
        digit_lower = product_lower // (scale * scale)
        digit_upper = product_upper // (scale * scale)
        if digit_lower != digit_upper:
            raise ArithmeticError("Increase precision: digit %d is ambiguous" % index)
        digit = digit_lower
        assert 0 <= digit <= 3
        lower = product_lower // scale - digit * scale
        upper = ceil_div(product_upper, scale) - digit * scale
        assert 0 <= lower <= upper < scale
        digits.append(digit)
        remainders.append((lower, upper))
        margin = min(lower, scale - upper)
        if least_margin is None or margin < least_margin[1]:
            least_margin = (index, margin)
    assert digits[0] == 3
    return digits, remainders, least_margin


def integer_reciprocal(digits):
    c = [1]
    for n in range(1, len(digits) + 1):
        c.append(sum(digits[k - 1] * c[n - k] for k in range(1, n + 1)))
        assert c[n] >= 3 ** n
    return c


def exact_polynomial_identities(digits, c):
    """Use a formal indeterminate X; do not substitute an approximate pi."""
    count = len(digits)
    remainders = [[1]]
    for n, digit in enumerate(digits, 1):
        remainders.append([-digit] + remainders[-1])
        difference = remainders[n][:]
        for degree, coefficient in enumerate(remainders[n - 1], 1):
            difference[degree] -= coefficient
        assert difference == [-digit] + [0] * n
    reciprocal_coefficients_checked = 0
    for n in range(count + 1):
        coefficient = c[n] - sum(digits[k - 1] * c[n - k] for k in range(1, n + 1))
        assert coefficient == (1 if n == 0 else 0)
        total = [0] * (n + 1)
        # H_j(X)=c_j-X*c_(j-1), with H_0=1.
        for k in range(n + 1):
            j = n - k
            for degree, value in enumerate(remainders[k]):
                total[degree] += c[j] * value
                if j:
                    total[degree + 1] -= c[j - 1] * value
        assert total == ([1] if n == 0 else [0] * (n + 1))
        reciprocal_coefficients_checked += len(total)
    return {
        "D_times_C_checked_through": count,
        "D_equals_1_minus_Xu_times_R_checked_through": count,
        "R_times_H_checked_through": count,
        "R_times_H_integer_polynomial_coefficients_checked": reciprocal_coefficients_checked,
        "H_definition": "h_0(X)=1; h_n(X)=c_n-X*c_(n-1)",
    }


def bounded_linear_forms(c, pi_lower, pi_upper, scale):
    for n in range(1, len(c)):
        lower = c[n] * scale - pi_upper * c[n - 1]
        upper = c[n] * scale - pi_lower * c[n - 1]
        assert max(abs(lower), abs(upper)) <= (1 << (n - 1)) * scale
    return len(c) - 1


def m_interval(remainders, pi_lower, pi_upper, scale):
    rho_lower = scale * scale // pi_upper
    rho_upper = ceil_div(scale * scale, pi_lower)
    power_lower = power_upper = scale
    total_lower = total_upper = 0
    for remainder_lower, remainder_upper in remainders:
        total_lower += remainder_lower * power_lower // scale
        total_upper += ceil_div(remainder_upper * power_upper, scale)
        power_lower = power_lower * rho_lower // scale
        power_upper = ceil_div(power_upper * rho_upper, scale)
    tail_upper = ceil_div(power_upper * scale, scale - rho_upper)
    total_upper += tail_upper
    assert total_lower > scale
    assert total_upper * (pi_upper - scale) < pi_upper * scale
    return total_lower, total_upper, tail_upper


def certify_coefficient_asymptotic(c, pi_lower, pi_upper, scale, m_lower, m_upper):
    # x_n encloses pi^n/m. Each step is rounded outwards to denominator scale.
    lower = scale * scale // m_upper
    upper = ceil_div(scale * scale, m_lower)
    for n, coefficient in enumerate(c):
        error_upper = max(abs(coefficient * scale - lower),
                          abs(upper - coefficient * scale))
        # A LOWER bound for 2^n/(pi-2), sufficient for the asserted inequality.
        bound_lower = (1 << n) * scale * scale // (pi_upper - 2 * scale)
        assert error_upper <= bound_lower
        lower = lower * pi_lower // scale
        upper = ceil_div(upper * pi_upper, scale)
    return len(c)


def ratio_data(n, c, pi_lower, pi_upper, scale):
    divisor = gcd(c[n + 1], c[n])
    numerator, denominator = c[n + 1] // divisor, c[n] // divisor
    error_lower = numerator * scale - denominator * pi_upper
    error_upper = numerator * scale - denominator * pi_lower
    if error_lower >= 0:
        absolute_lower, absolute_upper, side = error_lower, error_upper, "above_pi"
    elif error_upper <= 0:
        absolute_lower, absolute_upper, side = -error_upper, -error_lower, "below_pi"
    else:
        raise ArithmeticError("Increase precision: a ratio's error sign is ambiguous")
    kappa_lower = denominator * absolute_lower
    kappa_upper = denominator * absolute_upper
    # Independent equality of the reduced and unreduced normalization formulas.
    raw_lower = c[n + 1] * scale - c[n] * pi_upper
    raw_upper = c[n + 1] * scale - c[n] * pi_lower
    assert raw_lower == divisor * error_lower
    assert raw_upper == divisor * error_upper
    assert max(abs(raw_lower), abs(raw_upper)) * 3 ** n <= (1 << n) * scale * c[n]
    return {"n": n, "p": numerator, "q": denominator, "gcd": divisor,
            "raw_c_n": c[n], "side": side,
            "kappa_lower_scaled": kappa_lower, "kappa_upper_scaled": kappa_upper}


def compact_ratio(row, scale, bits, full_certificate=False):
    result = {"n": row["n"], "p": str(row["p"]), "q": str(row["q"]),
              "gcd": str(row["gcd"]), "side": row["side"],
              "raw_c_n_decimal_digits": len(str(row["raw_c_n"])),
              "kappa_outward_decimal_interval": display_interval(
                  row["kappa_lower_scaled"], row["kappa_upper_scaled"], scale)}
    if full_certificate:
        result.update(kappa_lower_numerator=str(row["kappa_lower_scaled"]),
                      kappa_upper_numerator=str(row["kappa_upper_scaled"]),
                      denominator_power_of_two=bits,
                      raw_c_n=str(row["raw_c_n"]),
                      normalization="kappa=(c_n/gcd^2)*abs(c_(n+1)-pi*c_n)")
    return result


def intrinsic_score(q, digits):
    """An exact score using q and a finite intrinsic digit prefix only.

    N is the least positive integer with 3**(N-1)>=q**3*2**(N-1).
    The result is q*distance(q*c_N/c_(N-1), Z), computed by one remainder.
    No pi enclosure or real-number operation is used in this function.
    """
    if not isinstance(q, int) or q < 1:
        raise ValueError("q must be a positive integer")
    n, two_power, three_power = 1, 1, 1
    while three_power < q ** 3 * two_power:
        n += 1
        two_power *= 2
        three_power *= 3
    if n > len(digits):
        raise ValueError("Need at least %d intrinsic digits for this q" % n)
    if n > 1:
        assert three_power // 3 < q ** 3 * (two_power // 2)
    prefix_c = integer_reciprocal(digits[:n])
    a, b = prefix_c[n - 1], prefix_c[n]
    residue = (q * b) % a
    distance_numerator = min(residue, a - residue)
    score = Fraction(q * distance_numerator, a)
    bound = Fraction(q * q * two_power, a)
    assert bound <= Fraction(1, q)
    return {"q": q, "N": n, "a": a, "b": b, "residue": residue,
            "score": score, "error_bound": bound}


def intrinsic_score_checks(digits, pi_lower, pi_upper, scale):
    denominators = list(range(1, 129)) + [33102, 10 ** 6]
    sampled = {1, 7, 106, 113, 128, 33102, 10 ** 6}
    results = []
    maximum_prefix = 0
    for q in denominators:
        native = intrinsic_score(q, digits)
        maximum_prefix = max(maximum_prefix, native["N"])
        # pi is used only below, to certify the independently computed score.
        nearest_lower = (2 * q * pi_lower + scale) // (2 * scale)
        nearest_upper = (2 * q * pi_upper + scale) // (2 * scale)
        assert nearest_lower == nearest_upper
        p = nearest_lower
        signed_lower = p * scale - q * pi_upper
        signed_upper = p * scale - q * pi_lower
        if signed_lower >= 0:
            kappa_lower, kappa_upper = q * signed_lower, q * signed_upper
        elif signed_upper <= 0:
            kappa_lower, kappa_upper = -q * signed_upper, -q * signed_lower
        else:
            raise ArithmeticError("Increase interval precision for a score certificate")
        true_lower, true_upper = Fraction(kappa_lower, scale), Fraction(kappa_upper, scale)
        score, bound = native["score"], native["error_bound"]
        assert score - bound <= true_lower <= true_upper <= score + bound
        assert score - Fraction(1, q) <= true_lower <= true_upper <= score + Fraction(1, q)
        # A sharper finite Lipschitz certificate from the actual ratio enclosure.
        ratio_error_lower = native["b"] * scale - native["a"] * pi_upper
        ratio_error_upper = native["b"] * scale - native["a"] * pi_lower
        actual_bound = Fraction(q * q * max(abs(ratio_error_lower), abs(ratio_error_upper)),
                                native["a"] * scale)
        assert score - actual_bound <= true_lower <= true_upper <= score + actual_bound
        if q in sampled:
            results.append({
                "q": q, "prefix_N": native["N"],
                "a_c_N_minus_1": str(native["a"]), "b_c_N": str(native["b"]),
                "residue_qb_mod_a": str(native["residue"]),
                "score_numerator": str(score.numerator), "score_denominator": str(score.denominator),
                "score_outward_decimal_interval": [scientific_fraction(score.numerator, score.denominator),
                    scientific_fraction(score.numerator, score.denominator, upper=True)],
                "universal_error_bound_numerator": str(bound.numerator),
                "universal_error_bound_denominator": str(bound.denominator),
                "universal_error_bound_le_1_over_q": True,
                "nearest_integer_to_q_pi": str(p),
                "true_q_distance_outward_decimal_interval": display_interval(kappa_lower, kappa_upper, scale),
                "actual_ratio_lipschitz_bound_upper": scientific_fraction(
                    actual_bound.numerator, actual_bound.denominator, upper=True),
                "whole_true_interval_inside_score_plus_minus_bound": True,
            })
    return {
        "definition": "S(q)=q/a*min((q*b) mod a, a-((q*b) mod a)); a=c_(N-1), b=c_N",
        "prefix_rule": "N is least >=1 with 3^(N-1)>=q^3*2^(N-1)",
        "theorem_bound": "abs(S(q)-q*distance(q*pi,Z))<=q^2*2^(N-1)/c_(N-1)<=1/q",
        "scope": "Finite checks of the score and error certificates. The asymptotic equivalence with P2848 requires the general proof; finite scores alone do not decide it.",
        "q_range": "1..128, 33102, 1000000", "case_count": len(denominators),
        "maximum_prefix_N": maximum_prefix,
        "all_certificates_passed": True, "selected_samples": results,
    }


def tail_polynomial_checks(digits, c):
    results = []
    for n in [0, 1, 2, 4, 8, 16, 32, 64]:
        if n + 1 > len(digits):
            continue
        coefficients = []
        for degree in range(len(digits) + 1):
            product = c[degree] if degree <= n else 0
            for j in range(min(n, degree - 1) + 1):
                product -= c[j] * digits[degree - j - 1]
            coefficients.append((1 if degree == 0 else 0) - product)
        assert all(value == 0 for value in coefficients[:n + 1])
        assert coefficients[n + 1] == c[n + 1]
        assert all(value >= 0 for value in coefficients[n + 1:])
        results.append({"N": n, "checked_through_degree": len(digits),
                        "first_nonzero_degree": n + 1,
                        "first_nonzero_coefficient": str(c[n + 1]),
                        "all_inspected_later_coefficients_nonnegative": True})
    return {
        "definition": "B_N=1-D*C_N, C_N=sum_(j=0)^N c_j u^j",
        "cases": results,
        "scope": "Finite coefficient checks only. Evaluation B_N(1/pi)=1 is an infinite-series identity proved separately, not established by this truncation.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits", type=int, default=256)
    parser.add_argument("--ratios", type=int, default=128)
    parser.add_argument("--bits", type=int, default=768)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("intrinsic_results.json"))
    args = parser.parse_args()
    if args.digits < 2 or not 1 <= args.ratios < args.digits:
        raise ValueError("Require digits >= 2 and 1 <= ratios < digits")
    if args.bits < 100:
        raise ValueError("Use at least 100 bits of interval precision")
    pi_lower, pi_upper, scale, machin = pi_fixed_interval(args.bits)
    digits, remainders, least_margin = greedy_digits(args.digits, pi_lower, pi_upper, scale)
    c = integer_reciprocal(digits)
    identities = exact_polynomial_identities(digits, c)
    h_count = bounded_linear_forms(c, pi_lower, pi_upper, scale)
    m_lower, m_upper, m_tail = m_interval(remainders, pi_lower, pi_upper, scale)
    asymptotic_count = certify_coefficient_asymptotic(c, pi_lower, pi_upper, scale, m_lower, m_upper)
    rows = [ratio_data(n, c, pi_lower, pi_upper, scale) for n in range(1, args.ratios + 1)]
    best = min(rows, key=lambda row: row["kappa_upper_scaled"])
    best_fraction = (best["p"], best["q"])
    different = [row for row in rows if (row["p"], row["q"]) != best_fraction]
    minimum_certified = all(best["kappa_upper_scaled"] < row["kappa_lower_scaled"] for row in different)
    assert minimum_certified
    records = []
    running = None
    for row in rows:
        if running is None or row["kappa_upper_scaled"] < running["kappa_lower_scaled"]:
            records.append(row)
            running = row
    sample_indices = sorted(({1, 2, 3, 4, 5, 6, 7, 8, 16, 32, 64, args.ratios, best["n"]}
                             & set(range(1, args.ratios + 1))))
    largest_gcd = max(rows, key=lambda row: row["gcd"])
    tail_checks = tail_polynomial_checks(digits, c)
    score_checks = intrinsic_score_checks(digits, pi_lower, pi_upper, scale)
    report = {
        "status": "all_exact_checks_passed",
        "scope": "Certified finite intrinsic pi-digit construction. No continued-fraction input, no claim about the full sequence's approximation constants, and no proof or disproof of P2848.",
        "arithmetic": "Python standard-library integers and fractions.Fraction; outward rational interval arithmetic only.",
        "pi_enclosure": {**machin, "denominator_power_of_two": args.bits,
                         "lower_numerator": str(pi_lower), "upper_numerator": str(pi_upper),
                         "outward_decimal_interval": display_interval(pi_lower, pi_upper, scale)},
        "greedy_digits": {"count": args.digits, "definition": "r_0=1, d_n=floor(pi*r_(n-1)), r_n=pi*r_(n-1)-d_n",
                          "digits": digits,
                          "sha256_of_digit_bytes": hashlib.sha256(bytes(digits)).hexdigest(),
                          "every_digit_interval_unambiguous": True,
                          "every_remainder_certified_in_0_1": True,
                          "least_boundary_margin_at_n": least_margin[0],
                          "least_boundary_margin_certified_lower": scientific_fraction(least_margin[1], scale),
                          "last_remainder_interval_width_upper": scientific_fraction(
                              remainders[-1][1] - remainders[-1][0], scale, upper=True)},
        "integer_reciprocal": {"definition": "c_0=1, c_n=sum_(k=1)^n d_k*c_(n-k)",
                               "count": len(c), "first_17_coefficients": [str(x) for x in c[:17]],
                               "sha256_of_comma_joined_decimal_coefficients": hashlib.sha256(
                                   ",".join(map(str, c)).encode()).hexdigest(),
                               "c_n_ge_3_to_n_checked_through": args.digits},
        "formal_identities": identities,
        "linear_form_bound": {"definition": "h_n=c_n-pi*c_(n-1)",
                              "bound": "abs(h_n)<=2^(n-1)", "certified_through": h_count},
        "renewal_constant": {"definition": "m=sum_(n>=0) r_n*pi^(-n)=sum_(k>=1) k*d_k*pi^(-k)",
                             "lower_numerator": str(m_lower), "upper_numerator": str(m_upper),
                             "denominator_power_of_two": args.bits,
                             "outward_decimal_interval": display_interval(m_lower, m_upper, scale),
                             "tail_upper_after_remainder_n": args.digits,
                             "tail_upper": scientific_fraction(m_tail, scale, upper=True),
                             "one_lt_m_lt_pi_over_pi_minus_one_certified": True},
        "coefficient_asymptotic_check": {"bound": "abs(c_n-pi^n/m)<=2^n/(pi-2)",
                                         "checked_n_zero_through": asymptotic_count - 1},
        "ratio_candidates": {"definition": "p=c_(n+1)/g, q=c_n/g, g=gcd(c_n,c_(n+1))",
                             "tested_n_min": 1, "tested_n_max": args.ratios,
                             "distinct_reduced_fraction_count": len({(r["p"], r["q"]) for r in rows}),
                             "ratio_error_bound_certified": "abs(c_(n+1)/c_n-pi)<=(2/3)^n",
                             "minimum_certified_over_tested_range": minimum_certified,
                             "minimum": compact_ratio(best, scale, args.bits, full_certificate=True),
                             "minimum_fraction_occurs_at_n": [r["n"] for r in rows if (r["p"], r["q"]) == best_fraction],
                             "strict_record_improvements": [compact_ratio(r, scale, args.bits) for r in records],
                             "largest_observed_gcd": {"n": largest_gcd["n"], "gcd": str(largest_gcd["gcd"])},
                             "selected_samples": [compact_ratio(rows[n - 1], scale, args.bits) for n in sample_indices]},
        "B_N_checks": tail_checks,
        "intrinsic_score_checks": score_checks,
    }
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "certified_digits": args.digits,
                      "ratio_count": args.ratios, "minimum": compact_ratio(best, scale, args.bits),
                      "intrinsic_score_case_count": score_checks["case_count"],
                      "maximum_score_prefix_N": score_checks["maximum_prefix_N"],
                      "m_interval": report["renewal_constant"]["outward_decimal_interval"],
                      "output": str(args.output)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
