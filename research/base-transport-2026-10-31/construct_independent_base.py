"""Finite certificates for a computable algebraically independent transport pair.

Run: python3 construct_independent_base.py --stages 16
Increase --stages to continue the same deterministic infinite construction.
Every decision below uses integers/Fraction; decimals are outward display only.
The infinite conclusion needs the termination proof in the research note.
No rational midpoint or finite decimal prefix is asserted to be transcendental.
"""

import argparse
import importlib.util
import itertools
import json
import math
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
HELPER = HERE.parent / "pi-native-space-2026-10-31" / "verify_native.py"
SPEC = importlib.util.spec_from_file_location("native_pi_bounds", HELPER)
NATIVE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(NATIVE)


def fraction_text(value):
    return f"{value.numerator}/{value.denominator}"


def interval_text(bounds):
    return [fraction_text(value) for value in bounds]


def display_decimal(value, rounding):
    with localcontext() as context:
        context.prec = 24
        context.rounding = rounding
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def interval_display(bounds):
    return {"lower": display_decimal(bounds[0], ROUND_FLOOR),
            "upper": display_decimal(bounds[1], ROUND_CEILING)}


def floor(value):
    return value.numerator // value.denominator


def ceil(value):
    return -floor(-value)


class CertifiedDigits:
    """Adaptive Machin bounds with dyadic outward remainder enclosures.

    Each retry doubles the working precision. No fixed cap on the number of
    terms, bits, digits, or retries is imposed. Pi's transcendence ensures that
    each finite greedy prefix eventually avoids all integer boundaries.
    """

    def __init__(self):
        self.digits = []
        self.certificate = None

    def get(self, count):
        if count <= len(self.digits):
            return self.digits[:count]
        bits = max(64, 4 * count + 32)
        retries = 0
        while True:
            terms = (bits + 3) // 4 + 8
            machin_low, machin_high = NATIVE.pi_bounds(terms)
            grid = 1 << bits
            pi_low = Q(floor(machin_low * grid), grid)
            pi_high = Q(ceil(machin_high * grid), grid)
            assert 3 < pi_low <= machin_low < machin_high <= pi_high < 4
            low = high = Q(1)
            digits = []
            for _ in range(count):
                raw_low, raw_high = pi_low * low, pi_high * high
                digit = floor(raw_low)
                if digit != floor(raw_high):
                    break
                assert digit in range(4)
                digits.append(digit)
                low = Q(floor((raw_low - digit) * grid), grid)
                high = Q(ceil((raw_high - digit) * grid), grid)
                assert 0 <= low <= high <= 1
            if len(digits) == count:
                assert digits[:len(self.digits)] == self.digits
                self.digits = digits
                self.certificate = {
                    "method": "Machin alternating rational bounds, adaptive precision",
                    "machin_terms_per_arctangent": terms,
                    "dyadic_bits": bits,
                    "retries_for_final_prefix": retries,
                    "pi_interval": interval_text((pi_low, pi_high)),
                    "digit_count": count,
                    "digits": "".join(map(str, digits)),
                }
                return digits
            bits *= 2
            retries += 1


def monomials(height):
    # Total degree first; within each degree, the exponent of T increases.
    return [(t_degree, degree - t_degree)
            for degree in range(height + 1)
            for t_degree in range(degree + 1)]


def polynomial_enumeration():
    """Primitive sign representatives, one height layer at a time.

    Height = max(total degree, largest absolute coefficient). A polynomial is
    retained iff its coefficient gcd is one and its first nonzero coefficient
    is positive. This covers the zero locus of every nonzero integer polynomial.
    Within each height, coefficient vectors are in lexicographic order.
    """
    height = 1
    while True:
        powers = monomials(height)
        length = len(powers)
        # This directly enumerates the vectors with positive first nonzero
        # coefficient in lexicographic order, skipping all negative-sign copies.
        for first in range(length - 1, -1, -1):
            for leading in range(1, height + 1):
                for rest in itertools.product(range(-height, height + 1),
                                              repeat=length - first - 1):
                    coeffs = (0,) * first + (leading,) + rest
                    if math.gcd(*coeffs) != 1:
                        continue
                    degree = max(sum(power) for coefficient, power in zip(coeffs, powers)
                                 if coefficient)
                    coefficient_height = max(map(abs, coeffs))
                    if max(degree, coefficient_height) != height:
                        continue
                    terms = [(coefficient, *power)
                             for coefficient, power in zip(coeffs, powers)
                             if coefficient]
                    yield {"height": height,
                           "dense_coefficients": list(coeffs),
                           "monomial_order": [list(power) for power in powers],
                           "terms": [list(term) for term in terms]}
        height += 1


def phi_enclosure(interval, digits):
    low_t, high_t = interval
    count = len(digits)
    low = sum((Q(digit) / high_t**n
               for n, digit in enumerate(digits, 1)), Q(0))
    high = sum((Q(digit) / low_t**n
                for n, digit in enumerate(digits, 1)), Q(0))
    high += 3 * low_t**(-count) / (low_t - 1)
    assert 0 < low < high
    return low, high


def polynomial_enclosure(polynomial, t_interval, y_interval):
    # Both input coordinates are positive, so each monomial is monotone in
    # each coordinate. Coefficients of either sign are then handled outward.
    t_low, t_high = t_interval
    y_low, y_high = y_interval
    assert 0 < t_low <= t_high and 0 < y_low <= y_high
    result_low = result_high = Q(0)
    for coefficient, t_degree, y_degree in polynomial["terms"]:
        monomial_low = t_low**t_degree * y_low**y_degree
        monomial_high = t_high**t_degree * y_high**y_degree
        if coefficient > 0:
            result_low += coefficient * monomial_low
            result_high += coefficient * monomial_high
        else:
            result_low += coefficient * monomial_high
            result_high += coefficient * monomial_low
    assert result_low <= result_high
    return result_low, result_high


def excludes_zero(interval):
    return interval[1] < 0 or interval[0] > 0


def choose_interval(stage, parent, polynomial, digit_oracle):
    precision = stage + 2
    candidates_tested = 0
    precisions_with_cells = 0
    while True:
        grid = 1 << precision
        # Require both endpoints strictly inside the parent.
        first_k = floor(parent[0] * grid) + 1
        last_k = ceil(parent[1] * grid) - 2
        if first_k <= last_k:
            digits = digit_oracle.get(precision)
            precisions_with_cells += 1
            for k in range(first_k, last_k + 1):
                candidate = (Q(k, grid), Q(k + 1, grid))
                image = phi_enclosure(candidate, digits)
                polynomial_image = polynomial_enclosure(polynomial, candidate, image)
                candidates_tested += 1
                if excludes_zero(polynomial_image):
                    return candidate, {
                        "stage": stage,
                        "polynomial": polynomial,
                        "precision_N": precision,
                        "dyadic_cell_index": k,
                        "interval": interval_text(candidate),
                        "Phi_interval": interval_text(image),
                        "P_interval": interval_text(polynomial_image),
                        "interval_decimal_display": interval_display(candidate),
                        "candidates_tested": candidates_tested,
                        "precisions_with_cells": precisions_with_cells,
                    }
        precision += 1


def verify_certificate(certificate):
    """Recheck the exact stored inequalities, prefix, and polynomial enumeration."""
    digit_cert = certificate["digit_certificate"]
    pi_low, pi_high = map(Q, digit_cert["pi_interval"])
    machin_low, machin_high = NATIVE.pi_bounds(digit_cert["machin_terms_per_arctangent"])
    assert pi_low <= machin_low < machin_high <= pi_high
    grid = 1 << digit_cert["dyadic_bits"]
    low = high = Q(1)
    digits = []
    for char in digit_cert["digits"]:
        raw_low, raw_high = pi_low * low, pi_high * high
        digit = int(char)
        assert floor(raw_low) == floor(raw_high) == digit
        low = Q(floor((raw_low - digit) * grid), grid)
        high = Q(ceil((raw_high - digit) * grid), grid)
        assert 0 <= low <= high <= 1
        digits.append(digit)
    assert len(digits) == digit_cert["digit_count"]

    parent = tuple(map(Q, certificate["initial_interval"]))
    polynomials = polynomial_enumeration()
    for stage, item in enumerate(certificate["stages"], 1):
        assert item["stage"] == stage
        assert item["polynomial"] == next(polynomials)
        interval = tuple(map(Q, item["interval"]))
        assert parent[0] < interval[0] < interval[1] < parent[1]
        assert interval[1] - interval[0] <= Q(1, 1 << (stage + 2))
        precision = item["precision_N"]
        index = item["dyadic_cell_index"]
        assert interval == (Q(index, 1 << precision), Q(index + 1, 1 << precision))
        assert precision <= len(digits)
        image = phi_enclosure(interval, digits[:precision])
        assert interval_text(image) == item["Phi_interval"]
        polynomial_image = polynomial_enclosure(item["polynomial"], interval, image)
        assert interval_text(polynomial_image) == item["P_interval"]
        assert excludes_zero(polynomial_image)
        parent = interval
    assert interval_text(parent) == certificate["final_interval"]
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stages", type=int, default=16)
    parser.add_argument("--output", type=Path,
                        default=HERE / "independent_base_certificate.json")
    args = parser.parse_args()
    if args.stages < 1:
        parser.error("--stages must be positive")

    initial = (Q(41, 4), Q(43, 4))
    current = initial
    oracle = CertifiedDigits()
    stages = []
    enumeration = polynomial_enumeration()
    for stage in range(1, args.stages + 1):
        current, record = choose_interval(stage, current, next(enumeration), oracle)
        stages.append(record)

    result = {
        "construction": "Deterministic nested dyadic intervals avoiding P_j(t,Phi(t))=0",
        "initial_interval": interval_text(initial),
        "enumeration": {
            "height": "max(total degree, maximum absolute coefficient)",
            "monomials": "increasing total degree; within a degree, increasing T exponent",
            "normalization": "coefficient gcd 1; first nonzero coefficient positive",
            "order": "increasing height, then lexicographic coefficient vector",
            "coverage": "Every nonzero integer polynomial has the same zero locus as one primitive sign representative enumerated here",
        },
        "search": {
            "precision": "N starts at j+2 and increases by one",
            "candidates": "consecutive 2^-N dyadic cells strictly inside the preceding interval, left to right",
            "acceptance": "rational interval evaluation of P_j excludes zero",
        },
        "digit_certificate": oracle.certificate,
        "stages": stages,
        "final_interval": interval_text(current),
        "final_interval_decimal_display": interval_display(current),
        "scope": "A finite prefix certificate for the defined infinite construction. The rational midpoint is not asserted to be transcendental. Algebraic independence of the infinite limit pair relies on the proved termination of every stage and enumeration of all polynomial zero loci; these finite stages alone do not establish it. No conclusion about Phi(10) or P2848.",
    }
    result["exact_certificate_recheck_passed"] = verify_certificate(result)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stages": args.stages,
                      "digits_certified": len(oracle.digits),
                      "final_interval": result["final_interval"],
                      "final_interval_decimal_display": result["final_interval_decimal_display"],
                      "exact_certificate_recheck_passed": result["exact_certificate_recheck_passed"],
                      "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
