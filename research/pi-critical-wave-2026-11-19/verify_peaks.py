"""Enclose the quadratic wave energy at exact analytic peak times.

Uses the independently checked pi interval and witnesses in the preceding study.
All comparisons in this script use integers and fractions.
"""
import argparse
import json
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path


def decimal_floor(value, digits):
    scale = 10**digits
    return value.numerator * scale // value.denominator


def decimal_ceil(value, digits):
    return -decimal_floor(-value, digits)


def sqrt_floor(value, digits):
    scale = 10**digits
    return isqrt((value.numerator * scale**2) // value.denominator)


def show(numerator, digits):
    if digits == 0:
        return str(numerator)
    whole, tail = divmod(numerator, 10**digits)
    return f"{whole}.{tail:0{digits}d}"


def arctan_bounds(inverse, terms):
    """Alternating series with an even number of retained terms."""
    assert terms % 2 == 0
    partial, power = Fraction(0), inverse
    for j in range(terms):
        partial += Fraction((-1)**j, (2*j+1)*power)
        power *= inverse**2
    return partial, partial + Fraction(1, (2*terms+1)*power)


def independent_pi_bounds():
    lower5, upper5 = arctan_bounds(5, 268)
    lower239, upper239 = arctan_bounds(239, 80)
    # Machin's identity. Both endpoints are rational.
    return 16*lower5 - 4*upper239, 16*upper5 - 4*lower239


def certify_peak(p, q, pi_lower, pi_upper, expected_energy, time_digits):
    divisor = gcd(p, 4*q)
    u, k = p//divisor, (4*q)//divisor
    assert gcd(u, k) == 1 and u > 0
    delta_lower, delta_upper = sorted((4*u-k*pi_lower, 4*u-k*pi_upper))
    assert delta_lower*delta_upper > 0
    error_lower, error_upper = sorted((abs(delta_lower), abs(delta_upper)))
    lambda2_lower, lambda2_upper = 2*u*error_lower, 2*u*error_upper
    assert lambda2_upper < pi_lower  # the primitive vector is the shortest
    cutoff = pi_lower**2 // lambda2_upper
    # exp(-pi^2/lambda^2) < 2^(-cutoff), from e>2.
    correction = Fraction(1, 2**cutoff)
    factor = ((1+correction)/(1-correction))**2
    lower_energy_squared = pi_lower/lambda2_upper
    upper_energy_squared = (pi_upper/lambda2_lower)*factor**2
    lo = sqrt_floor(lower_energy_squared, 5)
    hi = sqrt_floor(upper_energy_squared, 5) + 1
    assert (lo, hi) == expected_energy, (lo, hi)
    time_lo = decimal_floor(error_lower/(2*u), time_digits)
    time_hi = decimal_ceil(error_upper/(2*u), time_digits)
    assert time_lo > 0 and time_lo < time_hi
    print(f"q has {len(str(q))} digits, gcd(p,4q)={divisor}, "
          f"Z(t*) in [{show(lo,5)}, {show(hi,5)}]")
    print(f"  t* in [{time_lo}, {time_hi}] * 10^(-{time_digits})")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", nargs="?", type=Path,
                        default=Path(__file__).resolve().parents[1] /
                        "pi-first-hit-window-2026-11-18" / "certificate.json")
    args = parser.parse_args()
    data = json.loads(args.certificate.read_text())
    scale = 10**data["pi_grid_digits"]
    pi_lower = Fraction(data["pi_lower_integer"], scale)
    pi_upper = Fraction(data["pi_upper_integer"], scale)
    actual_lower, actual_upper = independent_pi_bounds()
    assert pi_lower <= actual_lower < actual_upper <= pi_upper
    previous, latest = data["records"][-2:]
    for record, bounds, digits in [(previous, (605777, 605778), 14),
                                    (latest, (2956084, 2956085), 318)]:
        certify_peak(record["p"], record["q"], pi_lower, pi_upper, bounds, digits)


if __name__ == "__main__":
    main()
