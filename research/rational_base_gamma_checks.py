"""Exact checks for rational-base fixed lines and denominator-scale waves.

Standard library only. Pi decisions use a Machin-series rational enclosure.
Gamma reflection is checked numerically only; its proof is not this script.
No infinite Diophantine conclusion is claimed by finite checks.
"""

from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from math import gcd, lgamma, exp, sin, pi, log
import json
import random


def atan_bounds(denominator, terms):
    total = sum((F((-1) ** k, (2 * k + 1) * denominator ** (2 * k + 1))
                 for k in range(terms)), F(0))
    following = F((-1) ** terms,
                  (2 * terms + 1) * denominator ** (2 * terms + 1))
    return min(total, total + following), max(total, total + following)


def pi_bounds():
    a, b = atan_bounds(5, 160)
    c, d = atan_bounds(239, 50)
    lower, upper = 16 * a - 4 * d, 16 * b - 4 * c
    scale = 10 ** 200
    lower = F(lower.numerator * scale // lower.denominator, scale)
    upper = F(-(-upper.numerator * scale // upper.denominator), scale)
    assert upper - lower == F(1, scale)
    return lower, upper


def decimal(number, precision=22):
    with localcontext() as context:
        context.prec = precision
        return str(Decimal(number.numerator) / Decimal(number.denominator))


def decimal_interval(lower, upper, precision=22):
    out = []
    with localcontext() as context:
        context.prec = precision
        for value, rounding in ((lower, ROUND_FLOOR), (upper, ROUND_CEILING)):
            context.rounding = rounding
            out.append(str(Decimal(value.numerator) / Decimal(value.denominator)))
    return out


def abs_interval(lower, upper):
    assert lower * upper > 0
    return min(abs(lower), abs(upper)), max(abs(lower), abs(upper))


def distance_to_integer(number):
    return abs(number - (number + F(1, 2)).__floor__())


def identities():
    rng = random.Random(20261028)
    for _ in range(8000):
        q0 = rng.randrange(2, 1000)
        y = F(rng.randrange(1, q0), q0)
        a, q = y.numerator, y.denominator
        base = 1 + 1 / y
        assert base == F(a + q, a)
        assert (base.numerator, base.denominator) == (a + q, a)
        assert 1 / (base - 1) == y
        assert (base * y).__floor__() == 1 and base * y - 1 == y
        x = F(rng.randrange(1, 999), 1000)
        target_base = 1 + 1 / x
        assert x - y == -x * a / q * (target_base - base)
        assert q * q * abs(x - y) == x * a * q * abs(target_base - base)
        n = rng.randrange(1, 100)
        assert distance_to_integer(n * base) == distance_to_integer(n * (base + F(1, n)))
        assert distance_to_integer(n * F(rng.randrange(1, 100), n)) == 0
        p, k = rng.randrange(1, 10000), rng.randrange(1, 10000)
        g = gcd(p, k)
        assert F(k * k, g * g) * abs(x - F(p, k)) == (k // g) ** 2 * abs(x - F(p, k))
    return 8000


def rational_aperiodic_example():
    remainder = F(1, 2)
    digits = []
    for step in range(1, 81):
        scaled = F(3, 2) * remainder
        digit = scaled.__floor__()
        remainder = scaled - digit
        assert remainder.denominator == 2 ** (step + 1)
        assert remainder.numerator % 2 == 1
        digits.append(str(digit))
    return ''.join(digits)


def analyze_pi():
    plo, phi = pi_bounds()
    xlo, xhi = plo - 3, phi - 3
    blo, bhi = 1 + 1 / xhi, 1 + 1 / xlo
    selected = [(22, 7), (333, 106), (355, 113),
                (103993, 33102), (104348, 33215),
                (245850922, 78256779)]
    rows = []
    for p, q in selected:
        assert gcd(p, q) == 1
        a = p - 3 * q
        y = F(a, q)
        assert F(1, 8) <= y <= F(1, 6)
        base = 1 + 1 / y
        lo, hi = xlo, xhi
        digits = []
        for _ in range(100):
            scaled_lo, scaled_hi = base * lo, base * hi
            digit = scaled_lo.__floor__()
            assert digit == scaled_hi.__floor__(), 'Insufficient pi precision'
            digits.append(digit)
            lo, hi = scaled_lo - digit, scaled_hi - digit
            if digit != 1:
                break
        else:
            raise AssertionError('No first departure found')
        length = len(digits) - 1
        elo, ehi = abs_interval(xlo - y, xhi - y)
        threshold = 1 - y if xlo > y else y
        assert base ** length * ehi < threshold
        assert base ** (length + 1) * elo > threshold
        # The uniform, deliberately coarse compact-interval shadow bound.
        assert F(1, 72) / base ** length < elo
        assert ehi < F(7, 8) / base ** length
        klo, khi = q * q * elo, q * q * ehi
        phase_lo, phase_hi = abs_interval(blo - base, bhi - base)
        rows.append({
            'rational': f'{p}/{q}', 'small_fraction': str(y),
            'rational_base': str(base), 'base_decimal': decimal(base),
            'leading_ones': length, 'next_digit': digits[-1],
            'kappa_interval': decimal_interval(klo, khi),
            'base_normalized_interval': decimal_interval(a * a * phase_lo, a * a * phase_hi),
            'shadow_excess_approx': length - 2 * log(q) / log(float(base)),
        })
    return {'target_base_interval': decimal_interval(blo, bhi, 35), 'rows': rows}


def gamma_sanity():
    maximum_relative_error = 0.0
    for n in range(2, 100):
        for a in range(1, n):
            phase = a / n
            gamma_wave = exp(-lgamma(phase) - lgamma(1 - phase))
            sine_wave = sin(pi * phase) / pi
            maximum_relative_error = max(maximum_relative_error,
                                         abs(gamma_wave / sine_wave - 1))
            assert abs(gamma_wave - sine_wave) < 2e-14
            distance = min(phase, 1 - phase)
            assert (2 / pi) * distance - 2e-14 <= gamma_wave <= distance + 2e-14
    return maximum_relative_error


if __name__ == '__main__':
    result = {'exact_identity_samples': identities(),
              'rational_aperiodic_prefix_base_3_over_2': rational_aperiodic_example(),
              'pi': analyze_pi(),
              'gamma_reflection_float_max_relative_error': gamma_sanity()}
    print(json.dumps(result, ensure_ascii=False, indent=2))
