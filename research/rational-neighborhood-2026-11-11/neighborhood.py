#!/usr/bin/env python3
"""Exact charts for the q^-2 ball of a reduced rational; stdlib only.

This computes classical continued-fraction/Farey coordinates, not new pi
approximations. Run without arguments for the stored 22/7 and 355/113 examples.
"""
from fractions import Fraction as F
from math import gcd
import argparse
import json


def fs(x):
    x = F(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def parents(p, q):
    if q < 2 or gcd(p, q) != 1:
        raise ValueError("Use a reduced p/q with q >= 2")
    s = pow(p, -1, q)
    r = (p * s - 1) // q
    return ((r, s), (p - r, q - s))


def chart(p, q, r, s, tail):
    tail = F(tail)
    return (p * tail + r) / (q * tail + s)


def describe(p, q):
    center = F(p, q)
    p, q = center.numerator, center.denominator
    record = {"center": fs(center), "p": p, "q": q,
              "ball": [fs(center - F(1, q*q)), fs(center + F(1, q*q))]}
    if q == 1:
        record["integer_case"] = {
            "representation": "x=p+u, irrational -1<u<1",
            "convergent_interval": [fs(center-F(1, 2)), fs(center+1)],
            "note": "The q>=2 two-parent chart is not used for integers."}
        return record
    par = parents(p, q)
    endpoints = []
    branches = []
    for side, (r, s) in zip(("left", "right"), par):
        endpoint = F(p+r, q+s)
        endpoints.append(endpoint)
        branches.append({"side": side, "parent": fs(F(r,s)), "r": r, "s": s,
                         "determinant_p_s_minus_q_r": p*s-q*r,
                         "tail_lower_bound_strict": fs(1-F(s,q)),
                         "convergent_tail_condition": "t > 1",
                         "convergent_outer_endpoint": fs(endpoint),
                         "convergent_radius": fs(F(1,q*(q+s)))})
    length = endpoints[1]-endpoints[0]
    record.update({"branches": branches,
                   "convergent_interval": list(map(fs,endpoints)),
                   "sharp_symmetric_convergent_radius": fs(F(1,q*(q+max(s for r,s in par)))),
                   "convergent_length_fraction": fs(length/F(2,q*q))})
    return record


def arctan_interval(inv, terms=80):
    total = sum((F((-1)**k, (2*k+1)*inv**(2*k+1)) for k in range(terms)), F(0))
    other = total + F((-1)**terms, (2*terms+1)*inv**(2*terms+1))
    return min(total, other), max(total, other)


def pi_interval():
    l5, u5 = arctan_interval(5)
    l239, u239 = arctan_interval(239)
    return 16*l5-4*u239, 16*u5-4*l239


def decimal_enclosure(lo, hi, digits=24):
    """Outward decimal bounds, using exact integers, including negatives."""
    scale = 10**digits
    a = (lo*scale).numerator//(lo*scale).denominator
    b = -((-hi*scale).numerator//(-hi*scale).denominator)
    def dec(n):
        sign = "-" if n < 0 else ""
        n = abs(n)
        return f"{sign}{n//scale}.{n%scale:0{digits}d}"
    return [dec(a), dec(b)]


def examples():
    lo, hi = pi_interval()
    p, q, r, s = 22, 7, 3, 1
    # On this lower branch x(t) increases, so the inverse also increases.
    def inverse(x):
        return (r-s*x)/(q*x-p)
    tlo, thi = inverse(lo), inverse(hi)
    assert F(15) < tlo < thi < F(16)
    candidate = chart(p,q,r,s,16)
    assert candidate == F(355,113)
    q22 = (q*q*(F(p,q)-hi), q*q*(F(p,q)-lo))
    q355 = (113**2*(candidate-hi), 113**2*(candidate-lo))
    tail_error = (16-thi,16-tlo)
    # Verify interval-wise that the transformed tail identity gives the
    # same bounds as direct pi subtraction (the expression decreases in t).
    def transferred(t):
        return (16-t)*F(113)/(7*t+1)
    assert transferred(thi) == q355[0]
    assert transferred(tlo) == q355[1]
    return {
        "status": "PASS",
        "scope": "Classical coordinates and a known-fraction example; no new record or infinite pi theorem",
        "centers": [describe(22,7),describe(355,113),describe(1,2),describe(-5,7),describe(3,1)],
        "pi_example": {
            "interval_method": "16 arctan(1/5)-4 arctan(1/239); 80 alternating terms each",
            "pi_interval_exact": [fs(lo),fs(hi)],
            "pi_interval_decimal": decimal_enclosure(lo,hi),
            "lower_22_over_7_tail_exact": [fs(tlo),fs(thi)],
            "tail_decimal": decimal_enclosure(tlo,thi),
            "tail_minus_16_abs_decimal": decimal_enclosure(*tail_error),
            "tail_16_image": fs(candidate),
            "quality_22_over_7_decimal": decimal_enclosure(*q22),
            "quality_355_over_113_decimal": decimal_enclosure(*q355)}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("p", type=int, nargs="?")
    parser.add_argument("q", type=int, nargs="?")
    args = parser.parse_args()
    if (args.p is None) != (args.q is None):
        parser.error("Provide both p and q, or neither")
    result = examples() if args.p is None else describe(args.p,args.q)
    print(json.dumps(result, ensure_ascii=False, indent=2))
