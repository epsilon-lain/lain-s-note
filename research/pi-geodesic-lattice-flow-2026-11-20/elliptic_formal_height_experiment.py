#!/usr/bin/env python3
"""Exact rational experiment for high-order elliptic formal-group cancellation.

Curve: E: y^2 = x^3 - 2, point P=(3,5), formal parameter xi=-x/y.
Uses exact Fraction arithmetic for elliptic multiples and B_m(P_s).
mpmath is used only for the real elliptic-log integral and decimal diagnostics.
"""
from fractions import Fraction
from math import comb, log, log10
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 70
E_P = (Fraction(3), Fraction(5))


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and y1 == -y2:
        return None
    if P == Q:
        if y1 == 0:
            return None
        slope = Fraction(3 * x1 * x1, 2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope * slope - x1 - x2
    y3 = slope * (x1 - x3) - y1
    return x3, y3


def mul(n, P=E_P):
    R = None
    Q = P
    while n:
        if n & 1:
            R = add(R, Q)
        Q = add(Q, Q)
        n >>= 1
    return R


def xi(P):
    return -P[0] / P[1]


def log10_fraction(x):
    return log10(x.numerator) - log10(x.denominator)


def log_height_fraction(x):
    return max(log(x.numerator), log(x.denominator))


def B_m(s, m, cache):
    z = cache[s]
    total = Fraction(0)
    for j in range(2, m + 3):
        lam = (-1) ** j * comb(m + 1, j - 1)
        total += lam * cache[j * s] / z
    return total


def main():
    # Exact points needed for s=15,26 (m<=4), and s=41 (m<=6).
    max_multiple = 8 * 41
    cache = {n: xi(mul(n)) for n in range(1, max_multiple + 1)}
    selected = []
    for s, max_m in ((15, 4), (26, 4), (41, 6)):
        z = cache[s]
        for m in range(1, max_m + 1):
            B = B_m(s, m, cache)
            err = abs(B - 1)
            selected.append({
                "s": s,
                "m": m,
                "abs_xi_sP": mp.nstr(abs(mp.mpf(z.numerator) / z.denominator), 30),
                "log10_abs_Bm_minus_1": log10_fraction(err),
                "log_height_Bm": log_height_fraction(B),
                "numerator_bits": B.numerator.bit_length(),
                "denominator_bits": B.denominator.bit_length(),
                "log10_abs_error_over_xi_power_m":
                    log10_fraction(err) - m * log10_fraction(abs(z)),
                "weight_l1_norm": 2 ** (m + 1) - 1,
            })

    # Real elliptic logarithm for P on E(R): u(P)=integral_x(P)^infinity dx/sqrt(x^3-2).
    root = mp.root(2, 3)
    integrand = lambda x: 1 / mp.sqrt(x ** 3 - 2)
    period = mp.quad(integrand, [root, mp.inf])
    uP = mp.quad(integrand, [mp.mpf(3), mp.inf])
    theta = uP / period
    s = 41
    nearest = int(mp.nint(s * theta))
    delta = abs(s * theta - nearest)
    z41 = abs(mp.mpf(cache[41].numerator) / cache[41].denominator)

    # Naive logarithmic height/canonical-height estimate at 41P and 82P.
    hvals = {}
    for n in (41, 82):
        x = mul(n)[0]
        h_x = max(log(x.numerator), log(x.denominator))
        hvals[str(n)] = {
            "log_height_x_nP": h_x,
            "h_x_over_2n2": h_x / (2 * n * n),
        }

    result = {
        "curve": "y^2 = x^3 - 2",
        "point": "(3,5)",
        "point_non_torsion_reason": "Lutz-Nagell: y^2=25 does not divide discriminant -1728",
        "formal_parameter": "xi=-x/y",
        "real_uniformization": {
            "period": mp.nstr(period, 35),
            "uP": mp.nstr(uP, 35),
            "theta_u_over_period": mp.nstr(theta, 35),
            "s": 41,
            "nearest_integer_to_s_theta": nearest,
            "distance_s_theta_to_integer": mp.nstr(delta, 35),
            "period_times_distance": mp.nstr(period * delta, 35),
            "abs_xi_41P": mp.nstr(z41, 35),
            "ratio_period_distance_over_2xi": mp.nstr(period * delta / (2 * z41), 25),
        },
        "height_diagnostics": hvals,
        "exact_fraction_experiments": selected,
        "scope": "Exact rational elliptic arithmetic for B_m; real elliptic-log values are high-precision numerical diagnostics. No D-H-K explicit constant has been evaluated.",
    }
    out = Path(__file__).with_name("elliptic_formal_height_experiment.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "curve": result["curve"],
        "real_uniformization": result["real_uniformization"],
        "experiments": [
            (r["s"], r["m"], r["log10_abs_Bm_minus_1"], r["log_height_Bm"])
            for r in selected
        ],
        "output": str(out),
    }, indent=2))


if __name__ == "__main__":
    main()
