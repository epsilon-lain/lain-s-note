#!/usr/bin/env python3
"""Certified finite brute-force search for rational approximations to pi.

Standard library only. No float pi, continued fractions, prime filter or network.
Run python3 search.py --max-q 10000000 --output results.json
"""
import argparse
import json
import sys
import time
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
from math import gcd
from pathlib import Path


class NeedPrecision(Exception):
    pass


def atan_bounds(a, terms):
    power, total = a, Fraction(0)
    for k in range(terms):
        total += Fraction((-1)**k, (2*k+1)*power)
        power *= a*a
    other = total + Fraction((-1)**terms, (2*terms+1)*power)
    return min(total, other), max(total, other)


def pi_bounds(bits):
    """Return integers L,U,D with L/D < pi < U/D, D=2**bits."""
    terms = bits//4 + 16
    a, b = atan_bounds(5, terms)
    c, d = atan_bounds(239, terms)
    lo, hi = 16*a-4*d, 16*b-4*c  # Machin's identity
    scale = 1 << bits
    lower = (lo.numerator*scale)//lo.denominator
    upper = -((-hi.numerator*scale)//hi.denominator)
    assert 3*scale < lower < upper < 4*scale
    return lower, upper, scale


def clipped(p, pmin, pmax):
    return max(pmin, p) if pmax is None else min(pmax, max(pmin, p))


def nearest(q, lower, upper, scale, pmin, pmax):
    a = clipped((q*lower + scale//2)//scale, pmin, pmax)
    b = clipped((q*upper + scale//2)//scale, pmin, pmax)
    if a != b:
        raise NeedPrecision(f"nearest numerator ambiguous at q={q}")
    return a


def error_numerators(p, q, lower, upper, scale):
    """Enclose |p-q*pi| by [a/scale,b/scale]."""
    a, b = p*scale-q*upper, p*scale-q*lower
    if a >= 0:
        return a, b
    if b <= 0:
        return -b, -a
    return 0, max(-a, b)


def improves(candidate, best, metric):
    """Tuples (p,q,a,b); common scale cancels in exact comparisons."""
    if best is None:
        return True
    p, q, a, b = candidate
    bp, bq, ba, bb = best
    if metric == "absolute":
        if p*bq == bp*q:  # same rational, exact tie
            return False
        if b*bq < ba*q:
            return True
        if a*bq > bb*q:
            return False
    else:
        if p == bp and q == bq:
            return False
        if b*q < ba*bq:
            return True
        if a*q > bb*bq:
            return False
    raise NeedPrecision(f"{metric} comparison ambiguous at {p}/{q}")


def scan_once(qmin, qmax, pmin, pmax, bits, exhaustive, progress):
    lower, upper, scale = pi_bounds(bits)
    absolute, normalized = [], []
    best_abs = best_norm = None
    evaluated = 0
    for q in range(qmin, qmax+1):
        if exhaustive:
            # Explicit baseline: evaluate EVERY pair, then retain the row winner.
            row = None
            for p in range(pmin, pmax+1):
                a, b = error_numerators(p, q, lower, upper, scale)
                candidate = (p, q, a, b)
                evaluated += 1
                if improves(candidate, row, "absolute"):
                    row = candidate
            candidate = row
        else:
            p = nearest(q, lower, upper, scale, pmin, pmax)
            a, b = error_numerators(p, q, lower, upper, scale)
            candidate = (p, q, a, b)
            evaluated += 1
        if improves(candidate, best_abs, "absolute"):
            best_abs = candidate
            absolute.append(candidate)
        if improves(candidate, best_norm, "normalized"):
            best_norm = candidate
            normalized.append(candidate)
        if progress and (q-qmin+1) % progress == 0:
            print(f"certified q={q}; best normalized {best_norm[0]}/{best_norm[1]}",
                  file=sys.stderr, flush=True)
    return {"absolute": absolute, "normalized": normalized,
            "pi": (lower, upper, scale), "bits": bits,
            "pairs_evaluated": evaluated}


def scan(qmax, qmin=1, pmin=1, pmax=None, bits=160,
         exhaustive=False, progress=0):
    if not (1 <= qmin <= qmax and pmin >= 1 and bits >= 8 and progress >= 0):
        raise ValueError("Require 1 <= min-q <= max-q, min-p >= 1, bits >= 8.")
    if pmax is not None and pmax < pmin:
        raise ValueError("max-p must be at least min-p")
    if exhaustive and pmax is None:
        raise ValueError("--exhaustive requires --max-p")
    start, retries = time.perf_counter(), 0
    while True:
        try:
            result = scan_once(qmin, qmax, pmin, pmax, bits, exhaustive, progress)
            result.update({"precision_retries": retries,
                           "elapsed_seconds": time.perf_counter()-start,
                           "range": {"qmin": qmin, "qmax": qmax,
                                     "pmin": pmin, "pmax": pmax},
                           "mode": "exhaustive_pairs" if exhaustive else "nearest_per_denominator"})
            return result
        except NeedPrecision as error:
            retries += 1
            bits *= 2
            if progress:
                print(f"{error}; restarting with {bits} bits", file=sys.stderr, flush=True)


def decimal_bound(n, d, rounding):
    with localcontext() as ctx:
        ctx.prec = 30
        ctx.rounding = rounding
        return str(Decimal(n)/Decimal(d))


def interval_json(a, b, d):
    return {"lower_numerator": str(a), "upper_numerator": str(b),
            "denominator": str(d),
            "lower_decimal": decimal_bound(a, d, ROUND_FLOOR),
            "upper_decimal": decimal_bound(b, d, ROUND_CEILING)}


def record_json(candidate, scale):
    p, q, a, b = candidate
    common = gcd(p, q)
    return {"p": p, "q": q, "reduced_p": p//common, "reduced_q": q//common,
            "absolute_error": interval_json(a, b, scale*q),
            "normalized_error": interval_json(a*q, b*q, scale)}


def native_digits(count):
    """Certify the greedy pi-unit digits with outward fixed-point intervals."""
    bits = max(160, 2*count+64)
    while True:
        lower, upper, scale = pi_bounds(bits)
        lo = hi = scale
        digits = []
        for _ in range(count):
            lo = (lo*lower)//scale
            hi = -((-hi*upper)//scale)
            digit = lo//scale
            if digit != hi//scale:
                break
            assert 0 <= digit <= 3
            digits.append(digit)
            lo -= digit*scale
            hi -= digit*scale
        if len(digits) == count:
            return digits
        bits *= 2


def migration_interval(p, q, digits):
    """Enclose F(p/q)-1, F(t)=sum d_j t^-j; only valid for p/q>1."""
    if p <= q:
        return {"not_evaluated": "The migration series requires p/q > 1."}
    base = Fraction(p, q)
    value, power = Fraction(0), Fraction(1)
    for digit in digits:
        power /= base
        value += digit*power
    lo = value-1
    hi = lo + 3*power/(base-1)
    denominator = lo.denominator*hi.denominator
    return interval_json(lo.numerator*hi.denominator,
                         hi.numerator*lo.denominator, denominator)


def report(result, digit_count=64):
    lower, upper, scale = result["pi"]
    data = {key: value for key, value in result.items()
            if key not in ("absolute", "normalized", "pi")}
    data.update({"certificate": "Machin alternating-series rational enclosure; exact integer comparisons",
                 "pi_interval": interval_json(lower, upper, scale),
                 "absolute_records": [record_json(x, scale) for x in result["absolute"]],
                 "normalized_records": [record_json(x, scale) for x in result["normalized"]],
                 "scope": "Finite range only; no new global approximation record or proof of P2848 is claimed. Normalized scores use the original q, not a denominator after reduction."})
    if digit_count:
        digits = native_digits(digit_count)
        data["native_unit_digit_prefix"] = ''.join(map(str, digits))
        data["winner_migrations"] = {}
        for metric in ("absolute", "normalized"):
            p, q = result[metric][-1][:2]
            data["winner_migrations"][metric] = {
                "p": p, "q": q, "F_p_over_q_minus_1": migration_interval(p, q, digits)}
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-q", type=int, required=True)
    parser.add_argument("--min-q", type=int, default=1)
    parser.add_argument("--min-p", type=int, default=1)
    parser.add_argument("--max-p", type=int)
    parser.add_argument("--bits", type=int, default=160)
    parser.add_argument("--native-digits", type=int, default=64)
    parser.add_argument("--exhaustive", action="store_true")
    parser.add_argument("--progress-every", type=int, default=1000000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.native_digits < 0:
        parser.error("--native-digits must be nonnegative")
    try:
        result = scan(args.max_q, args.min_q, args.min_p, args.max_p,
                      args.bits, args.exhaustive, args.progress_every)
    except ValueError as error:
        parser.error(str(error))
    data = report(result, args.native_digits)
    encoded = json.dumps(data, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")
    for label in ("absolute", "normalized"):
        item = data[label+"_records"][-1]
        bounds = item[label+"_error"]
        print(f"best {label}: {item['p']}/{item['q']}; "
              f"[{bounds['lower_decimal']}, {bounds['upper_decimal']}]", file=sys.stderr)
    print(f"{result['pairs_evaluated']} pairs; {result['bits']} bits; "
          f"{result['elapsed_seconds']:.3f} seconds", file=sys.stderr)


if __name__ == "__main__":
    main()
