#!/usr/bin/env python3
"""Exact EDS and formal-parameter denominator experiment for E:y^2=x^3-2.

The EDS recurrence computes large multiples with integers only. Small n are
cross-checked against the exact rational elliptic group law. The finite
experiment is not a proof of the asymptotic primitive-part theorem.
"""
from fractions import Fraction
from itertools import combinations
from math import comb, gcd, isqrt, log
import json
from pathlib import Path

NMAX = 720
P = (Fraction(3), Fraction(5))

ALT = {
    4: {2: 4, 3: -5, 5: 5, 6: -4, 7: 1},
    5: {2: 4, 3: -5, 4: 1, 5: 1, 6: 1, 7: 1, 8: -5, 9: 4, 10: -1},
    6: {2: 4, 3: -4, 4: -3, 5: 6, 8: -6, 9: 3, 10: 4, 11: -4, 12: 1},
}


def point_add(P1, P2):
    if P1 is None:
        return P2
    if P2 is None:
        return P1
    x1, y1 = P1
    x2, y2 = P2
    if x1 == x2 and y1 == -y2:
        return None
    if P1 == P2:
        if y1 == 0:
            return None
        slope = Fraction(3 * x1 * x1, 2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope * slope - x1 - x2
    y3 = slope * (x1 - x3) - y1
    return x3, y3


def point_mul(n):
    result, base = None, P
    while n:
        if n & 1:
            result = point_add(result, base)
        base = point_add(base, base)
        n >>= 1
    return result


def xi_point(Q):
    return -Q[0] / Q[1]


def build_eds(nmax):
    W = {0: 0, 1: 1, 2: 10, 3: 171, 4: -7660}
    for m in range(3, nmax - 1):
        numerator = (
            W[m + 1] * W[m - 1] * W[2] ** 2
            - W[3] * W[1] * W[m] ** 2
        )
        assert numerator % W[m - 2] == 0
        W[m + 2] = numerator // W[m - 2]
    return W


def phi(n, W):
    return 3 * W[n] ** 2 - W[n + 1] * W[n - 1]


def xi_eds(n, W):
    return Fraction(-2 * phi(n, W) * W[n] ** 2, W[2 * n])


def binomial_weights(m):
    return {
        j: (-1) ** j * comb(m + 1, j - 1)
        for j in range(2, m + 3)
    }


def formal_combination(s, weights, W):
    base = xi_eds(s, W)
    return sum(
        (lam * xi_eds(j * s, W) / base for j, lam in weights.items()),
        Fraction(0),
    )


def strip_bad_primes(n):
    while n % 2 == 0:
        n //= 2
    while n % 3 == 0:
        n //= 3
    return n


def main():
    W = build_eds(NMAX)

    # Cross-check the EDS coordinate formula against exact group arithmetic.
    coordinate_checks = []
    for n in range(1, 13):
        Q = point_mul(n)
        direct = xi_point(Q)
        via_eds = xi_eds(n, W)
        assert direct == via_eds
        x, y = Q
        B = isqrt(x.denominator)
        assert B * B == x.denominator
        assert y.denominator == B ** 3
        assert abs(W[n]) == B
        coordinate_checks.append({"n": n, "B_n_bits": B.bit_length()})

    # Strong divisibility on a finite grid.
    strong_divisibility_checks = 0
    for m in range(1, 41):
        for n in range(1, 41):
            assert gcd(abs(W[m]), abs(W[n])) == abs(W[gcd(m, n)])
            strong_divisibility_checks += 1

    # Exact denominator data for the ratios R_j(15).
    s = 15
    ratio_denominators = {}
    for j in range(2, 13):
        ratio = xi_eds(j * s, W) / xi_eds(s, W)
        ratio_denominators[j] = ratio.denominator

    coprime_pair_checks = []
    for j, k in combinations(range(2, 13), 2):
        if gcd(j, k) == 1:
            common = gcd(ratio_denominators[j], ratio_denominators[k])
            assert strip_bad_primes(common) == 1
            coprime_pair_checks.append({
                "j": j, "k": k, "common_gcd_bits": common.bit_length(),
                "common_gcd": str(common),
            })

    ratio_bits = {
        str(j): d.bit_length()
        for j, d in ratio_denominators.items()
    }

    # Compare exact reduced rational heights for the old and alternative weights.
    comparisons = []
    for m in (4, 5, 6):
        for label, weights in (
            ("binomial", binomial_weights(m)),
            ("height_proxy_alternative", ALT[m]),
        ):
            assert all(
                sum(lam * j ** r for j, lam in weights.items()) == 1
                for r in range(m + 1)
            )
            value = formal_combination(s, weights, W)
            comparisons.append({
                "m": m,
                "scheme": label,
                "largest_multiplier": max(weights),
                "support": weights,
                "denominator_bits": value.denominator.bit_length(),
                "numerator_bits": abs(value.numerator).bit_length(),
                "log_denominator": log(value.denominator),
                "log_height": max(log(abs(value.numerator)), log(value.denominator)),
                "log10_abs_error": (
                    log(abs((value - 1).numerator))
                    - log(abs((value - 1).denominator))
                ) / log(10),
            })

    h_estimates = {
        str(n): log(abs(W[n])) / (n * n)
        for n in (100, 200, 300, 400, 500, 600, 720)
    }

    result = {
        "curve": "y^2=x^3-2",
        "point": "(3,5)",
        "formal_parameter": "xi=-x/y",
        "eds_initial_values": [W[n] for n in range(1, 5)],
        "eds_recurrence": (
            "W_(m+2) W_(m-2) = W_(m+1) W_(m-1) W_2^2 "
            "- W_3 W_1 W_m^2"
        ),
        "eds_cutoff": NMAX,
        "coordinate_formula_checks": coordinate_checks,
        "strong_divisibility_checks_for_1_to_40": strong_divisibility_checks,
        "ratio_denominator_bits_s_15": ratio_bits,
        "coprime_multiplier_pair_checks": coprime_pair_checks,
        "exact_weight_comparisons_s_15": comparisons,
        "log_Bn_over_n_squared_estimates": h_estimates,
        "scope": (
            "All EDS integer recurrences, small-point cross-checks, and "
            "rational combinations are exact. The primitive-part asymptotic "
            "bound is a literature theorem; this script does not factor "
            "the large EDS terms or independently certify their primitive parts."
        ),
    }
    out = Path(__file__).with_name("eds_primitive_denominator_experiment.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "eds_cutoff": NMAX,
        "strong_divisibility_checks": strong_divisibility_checks,
        "ratio_denominator_bits_s15": ratio_bits,
        "comparisons": [
            (r["m"], r["scheme"], r["largest_multiplier"],
             r["denominator_bits"], r["log10_abs_error"])
            for r in comparisons
        ],
        "output": str(out),
    }, indent=2))


if __name__ == "__main__":
    main()
