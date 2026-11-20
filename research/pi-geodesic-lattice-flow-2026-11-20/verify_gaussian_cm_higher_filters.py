#!/usr/bin/env python3
"""Exact higher-order Gaussian-CM trace filters for the square-CM pi bridge.

Imports exact elliptic arithmetic and Machin interval helpers from
verify_gaussian_cm_trace_filter.py. All weights, moments, point evaluations,
filters, and rational pi approximants are exact Fraction computations.
"""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path

from verify_gaussian_cm_trace_filter import (
    rational_point, xi, trace_multiplier, xi_rational,
    local_zeta_interval, error_log_interval, rational_cocycles,
    machin_pi_interval, pi_error_log_interval,
)


def gmul(x, y):
    a, b = x
    c, d = y
    return a*c-b*d, a*d+b*c


def gadd(x, y):
    return x[0]+y[0], x[1]+y[1]


def gsub(x, y):
    return x[0]-y[0], x[1]-y[1]


def gscale(k, x):
    return k*x[0], k*x[1]


def gpow(x, n):
    out = (1, 0)
    for _ in range(n):
        out = gmul(out, x)
    return out


def real_part_polynomial(a, b, degree):
    alpha = (a, b)
    a4 = gpow(alpha, 4)
    a8 = gpow(alpha, 8)
    if degree == 3:
        factor = gadd(gsub(gscale(11, a8), gscale(106, a4)), (375, 0))
    elif degree == 4:
        a12 = gpow(alpha, 12)
        factor = gadd(
            gsub(gadd(gsub(gscale(633, a12), gscale(9091, a8)),
                      gscale(55883, a4)), (0, 0)),
            (-170625, 0),
        )
    elif degree == 5:
        a12 = gpow(alpha, 12)
        a16 = gpow(alpha, 16)
        factor = gadd(
            gsub(gadd(gsub(gscale(1607, a16), gscale(30676, a12)),
                      gscale(265906, a8)), gscale(1326212, a4)),
            (3729375, 0),
        )
    else:
        raise ValueError(degree)
    value = gmul(gmul(alpha, gsub(a4, (1, 0))), factor)
    return F(value[0], a)


def r1(a, b):
    return a**4-10*a*a*b*b+5*b**4-1


def r2(a, b):
    return (a**8-36*a**6*b*b+126*a**4*b**4-84*a*a*b**6+9*b**8
            -6*a**4+60*a*a*b*b-30*b**4+5)


def r3(a, b):
    return real_part_polynomial(a, b, 3)


def r4(a, b):
    return real_part_polynomial(a, b, 4)


def r5(a, b):
    return real_part_polynomial(a, b, 5)


def filter_value(n, entries):
    Q = rational_point(n)
    z = xi(Q).r
    F_value = sum((weight*trace_multiplier(Q, a, b)
                   for a, b, weight in entries), F(0))
    return F_value/z


def L_trunc(z, K):
    return sum((F((-1)**r * __import__("math").comb(2*r, r)*36**r, 4*r+1)
                * z**(4*r+1) for r in range(K+1)), F(0))


def pi_approximant(entries, K):
    z7, z11 = xi_rational(7), xi_rational(11)
    R = rational_cocycles(11)
    omega = 7*L_trunc(z11, K)-11*L_trunc(z7, K)
    H = (11*R[7]-7*R[11]
         -11*filter_value(7, entries)+7*filter_value(11, entries))
    return omega*H


ORDER19 = [
    (1, 0, F(202121, 178717)),
    (1, 1, F(-72883, 536151)),
    (2, 0, F(8893, 1072302)),
    (1, 2, F(-661, 1072302)),
    (2, 1, F(-1445, 536151)),
]

ORDER23 = [
    (1, 0, F(12486023931, 10913532322)),
    (1, 1, F(-7249221697, 49110895449)),
    (2, 0, F(31564637, 4677228138)),
    (1, 2, F(-56405183, 49110895449)),
    (2, 1, F(-103699966, 49110895449)),
    (2, 2, F(1685917, 49110895449)),
]


def verify_moments(entries, max_order):
    moment_functions = [lambda a,b: 1, r1, r2, r3, r4, r5]
    targets = [F(1), F(5,6), F(-45,14), F(4875,22),
               F(-193375,2), F(78316875,38)]
    for k in range(max_order+1):
        value = sum((w*moment_functions[k](a,b)
                     for a,b,w in entries), F(0))
        assert value == targets[k], (k, value, targets[k])


def summarize_filter(name, entries, max_order, K):
    verify_moments(entries, max_order)
    local = []
    for n in (7, 11):
        z = xi_rational(n)
        value = filter_value(n, entries)
        true_interval = local_zeta_interval(z, 30)
        local.append({
            "n": n,
            "numerator_bits": abs(value.numerator).bit_length(),
            "denominator_bits": value.denominator.bit_length(),
            "log10_local_error_interval":
                error_log_interval(value, true_interval),
        })
    pi_hat = pi_approximant(entries, K)
    pi_err, side = pi_error_log_interval(pi_hat, machin_pi_interval(100))
    return {
        "name": name,
        "weights": [
            {"a":a,"b":b,"weight":str(w)} for a,b,w in entries
        ],
        "largest_gaussian_norm": max(a*a+b*b for a,b,_ in entries),
        "moment_order": max_order,
        "local_zeta_error_order": 4*max_order+3,
        "local_results": local,
        "pi_approximant": {
            "L_truncation_K": K,
            "numerator_bits": abs(pi_hat.numerator).bit_length(),
            "denominator_bits": pi_hat.denominator.bit_length(),
            "log10_abs_error_interval": pi_err,
            "side_of_pi": side,
        },
    }


def main():
    order19 = summarize_filter("order19", ORDER19, 4, 3)
    order23 = summarize_filter("order23", ORDER23, 5, 4)
    result = {
        "curve": "E:y^2=x^3-36x",
        "formal_parameter": "xi=-x/y",
        "order19_target_r4": "-193375/2",
        "order23_target_r5": "78316875/38",
        "filters": [order19, order23],
        "scope": (
            "All Gaussian moment identities, elliptic evaluations, and "
            "rational pi approximants are exact. Local zeta errors use "
            "alternating-series intervals; pi errors use a 100-term exact "
            "Machin interval."
        ),
    }
    out = Path(__file__).with_name("gaussian_cm_higher_filter_results.json")
    out.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
