#!/usr/bin/env python3
"""Certified intrinsic pi-base addition, with no decimal expansion as input.

The inputs are the finite pi-base words 0.02 and 0.02. The output is the
greedy pi-base expansion of their sum 4/pi^2. Rational enclosures for pi
come from Machin's identity and alternating arctangent sums. Every floor
operation is certified by an exact interval. Finite computations do not
prove aperiodicity: the separate transcendence argument does that.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from pathlib import Path


def arctangent_interval(denominator: int, terms: int):
    """Closed rational bounds for arctan(1/denominator)."""
    if denominator <= 1 or terms < 1:
        raise ValueError("need denominator > 1 and terms >= 1")
    value = sum((F((-1) ** k, (2 * k + 1) * denominator ** (2 * k + 1))
                 for k in range(terms)), F(0))
    next_value = value + F((-1) ** terms,
                          (2 * terms + 1) * denominator ** (2 * terms + 1))
    return min(value, next_value), max(value, next_value)


def pi_interval(terms5: int, terms239: int):
    lower5, upper5 = arctangent_interval(5, terms5)
    lower239, upper239 = arctangent_interval(239, terms239)
    lower = 16 * lower5 - 4 * upper239
    upper = 16 * upper5 - 4 * lower239
    assert F(3) < lower < upper < F(4)
    return lower, upper


def floor_fraction(value: F):
    return value.numerator // value.denominator


def floor_from_interval(lower: F, upper: F):
    if lower > upper:
        raise ValueError("reversed interval")
    digit = floor_fraction(lower)
    if not upper < digit + 1:
        raise ArithmeticError("pi enclosure is too wide to certify this digit")
    return digit


def dyadic_lower(value: F, bits: int = 64):
    return F(floor_fraction(value * (1 << bits)), 1 << bits)


def dyadic_upper(value: F, bits: int = 64):
    scaled = value * (1 << bits)
    ceiling = -((-scaled.numerator) // scaled.denominator)
    return F(ceiling, 1 << bits)


def short_interval(lower, upper):
    return {"lower": str(dyadic_lower(lower)),
            "upper": str(dyadic_upper(upper))}


def exact_fraction(value: F):
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def certify_unit_digits(count, beta_lower, beta_upper):
    """Independently certify the native greedy expansion of the unit 1."""
    remainder_lower = remainder_upper = F(1)
    digits = []
    for _ in range(count):
        product_lower = beta_lower * remainder_lower
        product_upper = beta_upper * remainder_upper
        digit = floor_from_interval(product_lower, product_upper)
        assert 0 <= digit <= 3
        remainder_lower = product_lower - digit
        remainder_upper = product_upper - digit
        assert 0 < remainder_lower < remainder_upper < 1
        digits.append(digit)
    return digits


def check_principal_zero_relation(output_digits, beta_lower, beta_upper):
    """Check the concrete zero relation F=(4u^2-E(u))=D(u)A(u).

    The unit digits d_n and output digits e_n are certified independently.
    Infinite evaluation vanishes by the exact greedy expansion of 4/pi^2;
    the finite products below are identities modulo u^(N+1), not assertions
    that truncated F has zero value.
    """
    count = len(output_digits)
    unit_digits = certify_unit_digits(count, beta_lower, beta_upper)
    f_coefficients = [0] + [4 * (n == 2) - output_digits[n - 1]
                            for n in range(1, count + 1)]
    assert max(map(abs, f_coefficients)) <= 4
    a_coefficients = []
    for n in range(count + 1):
        next_a = f_coefficients[n]
        for k in range(1, n + 1):
            next_a += unit_digits[k - 1] * a_coefficients[n - k]
        a_coefficients.append(next_a)

    d_coefficients = [1] + [-digit for digit in unit_digits]
    # Dense convolution is separate from the coefficient-by-coefficient solve.
    product = [0] * (2 * count + 1)
    for i, left in enumerate(d_coefficients):
        for j, right in enumerate(a_coefficients):
            product[i + j] += left * right
    assert product[:count + 1] == f_coefficients

    # beta_upper > pi makes this a LOWER bound for 4*2^n/(pi-1).
    # Satisfying this stricter rational comparison certifies the desired bound.
    for n, value in enumerate(a_coefficients):
        assert abs(value) * (beta_upper - 1) <= 4 * (1 << n)
    max_abs = max(map(abs, a_coefficients))
    max_indices = [n for n, value in enumerate(a_coefficients) if abs(value) == max_abs]
    normalized = [F(abs(value), 1 << n) for n, value in enumerate(a_coefficients)]
    return {
        "status": "pass",
        "F": "4*u^2 - sum_{n>=1} e_n*u^n, where e_n are this addition's canonical output digits",
        "D": "1 - sum_{n>=1} d_n*u^n, where d_n are the independently certified greedy pi-base digits of 1",
        "infinite_identity": "F=D*A in the principal relation kernel (proved separately)",
        "finite_identity_checked": "D*A=F modulo u^(N+1)",
        "N": count,
        "unit_floor_tests": count,
        "product_coefficients_checked": count + 1,
        "quantitative_bound": "abs(a_n)<=4*2^n/(pi-1)",
        "bound_coefficients_checked": count + 1,
        "rational_certification": "abs(a_n)*(pi_upper-1)<=4*2^n; pi_upper is an exact rational upper enclosure",
        "unit_digits": unit_digits,
        "quotient_coefficients_a0_through_aN": a_coefficients,
        "first_coefficients_a0_through_a32": a_coefficients[:33],
        "largest_absolute_coefficient": max_abs,
        "indices_attaining_largest_absolute_coefficient": max_indices,
        "largest_abs_a_n_divided_by_2_power_n": str(max(normalized)),
        "limits": "finite coefficient checks do not prove the whole principal-kernel theorem; they check its use on one genuine intrinsic addition relation",
    }


def certify(length: int, terms5: int, terms239: int):
    if length < 2:
        raise ValueError("length must be at least two")
    beta_lower, beta_upper = pi_interval(terms5, terms239)

    # The inputs 2/beta^2 each have greedy digits 0,2,0,... because beta>3.
    # For their sum: floor(4/beta)=1 and floor(4-beta)=0.
    first = floor_from_interval(4 / beta_upper, 4 / beta_lower)
    second = floor_from_interval(4 - beta_upper, 4 - beta_lower)
    assert (first, second) == (1, 0)
    digits = [first, second]
    carry_polynomial = [4, -1]  # Coefficients of C_2=4-beta, in increasing degree.
    carry_lower, carry_upper = 4 - beta_upper, 4 - beta_lower
    assert 0 < carry_lower < carry_upper < 1
    early_states = [
        {"n": 1, "polynomial_coefficients_in_increasing_degree": [-1],
         "interval": {"lower": "-1", "upper": "-1"}},
        {"n": 2, "polynomial_coefficients_in_increasing_degree": carry_polynomial[:],
         "interval": short_interval(carry_lower, carry_upper)},
    ]
    minimum_floor_margin = F(1)
    for n in range(3, length + 1):
        candidate_lower = beta_lower * carry_lower
        candidate_upper = beta_upper * carry_upper
        digit = floor_from_interval(candidate_lower, candidate_upper)
        assert 0 <= digit <= 3
        minimum_floor_margin = min(minimum_floor_margin,
                                   candidate_lower - digit,
                                   digit + 1 - candidate_upper)
        carry_lower, carry_upper = candidate_lower - digit, candidate_upper - digit
        assert 0 < carry_lower < carry_upper < 1
        carry_polynomial = [-digit] + carry_polynomial
        assert len(carry_polynomial) == n and carry_polynomial[-1] == -1
        digits.append(digit)
        if n <= 10:
            early_states.append({"n": n,
                                 "polynomial_coefficients_in_increasing_degree": carry_polynomial[:],
                                 "interval": short_interval(carry_lower, carry_upper)})

    # Independently check the formal coefficient recurrence from the input and
    # output digits, including the carry at the first position.
    direct_polynomial = []
    for n, digit in enumerate(digits, start=1):
        input_digit_sum = 4 if n == 2 else 0
        direct_polynomial = [input_digit_sum - digit] + direct_polynomial
        assert len(direct_polynomial) == n and direct_polynomial[-1] == -1
    assert direct_polynomial == carry_polynomial
    return {
        "status": "pass",
        "definition": "greedy beta expansion with beta=pi, digits {0,1,2,3}",
        "inputs": ["0.02 followed by zeros", "0.02 followed by zeros"],
        "input_value_each": "2*pi^(-2)",
        "sum_value": "4*pi^(-2)",
        "certified_digit_count": length,
        "sum_prefix": "0." + "".join(map(str, digits)),
        "pi_source": "pi=16*arctan(1/5)-4*arctan(1/239); alternating rational bounds",
        "arctangent_terms": {"denominator5": terms5, "denominator239": terms239},
        "pi_enclosure": {"lower": exact_fraction(beta_lower), "upper": exact_fraction(beta_upper)},
        "certification": {
            "floor_tests": length,
            "all_floors_unambiguous": True,
            "formal_recurrence_checks": length,
            "carry_degree_at_n": "n-1",
            "carry_leading_coefficient": -1,
            "floor_margin_lower_dyadic": str(dyadic_lower(minimum_floor_margin)),
            "end_carry_enclosure": short_interval(carry_lower, carry_upper),
            "end_carry_width_upper_dyadic_256bits": str(dyadic_upper(carry_upper - carry_lower, 256)),
        },
        "early_carry_states": early_states,
        "principal_zero_relation": check_principal_zero_relation(digits, beta_lower, beta_upper),
        "exact_remainder_identity": "4/pi^2 = sum_{k=1}^N digit[k]*pi^(-k) + C_N*pi^(-N), 0<C_N<1",
        "proved_separately_not_by_finite_scan": [
            "All nonzero numerical carry states are distinct because their integer polynomials have different degrees and pi is transcendental.",
            "The canonical output is not eventually periodic: otherwise its rational-function expression would give a nonzero polynomial relation for pi.",
        ],
        "limitations": [
            "No statistical randomness, equidistribution, or complexity growth rate is inferred.",
            "No conclusion about P2848 follows from this finite computation or the generic transcendental-base carry obstruction.",
            "The carry-state obstruction concerns exact numerical carry states; it is not a claim that all algorithms for pi-base arithmetic are impossible.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits", type=int, default=128)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("beta_shift_carry_results.json"))
    args = parser.parse_args()
    terms5, terms239 = args.digits // 2 + 32, args.digits // 8 + 16
    for attempt in range(5):
        try:
            result = certify(args.digits, terms5, terms239)
            break
        except ArithmeticError:
            terms5 *= 2
            terms239 *= 2
    else:
        raise RuntimeError("could not certify all requested digits")
    result["precision_attempts"] = attempt + 1
    sibling_report = Path(__file__).with_name("intrinsic_results.json")
    if sibling_report.exists():
        certified_elsewhere = json.loads(sibling_report.read_text(encoding="utf-8"))["greedy_digits"]["digits"]
        independently_certified = result["principal_zero_relation"]["unit_digits"]
        overlap = min(len(certified_elsewhere), len(independently_certified))
        assert certified_elsewhere[:overlap] == independently_certified[:overlap]
        result["principal_zero_relation"]["sibling_unit_digit_cross_check"] = {
            "file": sibling_report.name, "count": overlap, "status": "pass"}
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"],
                      "digits": result["certified_digit_count"],
                      "sum_prefix": result["sum_prefix"],
                      "zero_relation_coefficients_checked": result["principal_zero_relation"]["product_coefficients_checked"],
                      "largest_zero_relation_quotient_coefficient": result["principal_zero_relation"]["largest_absolute_coefficient"],
                      "precision_attempts": result["precision_attempts"],
                      "output": str(args.output)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
