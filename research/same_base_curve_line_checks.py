"""Exact finite checks for the same-base digit curve and its all-one line.

Uses the sibling rational_base_gamma_checks module for a certified pi interval.
All digit and area decisions use fractions. Pattern-root computations are
finite consistency checks, not proofs of asymptotic Diophantine claims.
"""

from fractions import Fraction as F
import json
import random

from rational_base_gamma_checks import pi_bounds, decimal_interval


def exact_digits(x, base, length):
    remainder = x
    digits, remainders = [], [x]
    for _ in range(length):
        product = base * remainder
        digit = product.__floor__()
        remainder = product - digit
        digits.append(digit)
        remainders.append(remainder)
    return digits, remainders


def area_interval(xlo, xhi, base, length=80):
    line = 1 / (base - 1)
    lo, hi = xlo, xhi
    partial, zero_partial = F(0), F(0)
    weight = F(1)
    digits = []
    for _ in range(length):
        a, b = base * lo, base * hi
        digit = a.__floor__()
        assert digit == b.__floor__(), 'Digit interval needs more precision'
        lo, hi = a - digit, b - digit
        digits.append(digit)
        weight /= base
        partial += abs(digit - 1) * weight
        if digit == 0:
            zero_partial += weight
    signed_lo, signed_hi = xlo - line, xhi - line
    lower = max(partial, signed_lo + 2 * zero_partial)
    upper = min(partial + weight * (hi + line),
                signed_hi + 2 * zero_partial + 2 * line * weight)
    assert 0 <= lower <= upper
    return lower, upper, digits


def exact_finite_checks():
    rng = random.Random(284820261028)
    samples = 800
    for _ in range(samples):
        denominator = rng.randrange(1, 41)
        base = 7 + F(rng.randrange(0, 2 * denominator + 1), denominator)
        x = F(rng.randrange(1, 1000), 1000)
        line = 1 / (base - 1)
        digits, remainders = exact_digits(x, base, 60)
        signed = x - line
        prefix = 0
        while prefix < len(digits) and digits[prefix] == 1:
            prefix += 1
        partial = sum((abs(d - 1) / base ** j for j, d in enumerate(digits, 1)), F())
        zeros = sum((1 / base ** j for j, d in enumerate(digits, 1) if d == 0), F())
        assert x == sum((d / base ** j for j, d in enumerate(digits, 1)), F()) + remainders[-1] / base ** 60
        assert partial == signed + 2 * zeros + (line - remainders[-1]) / base ** 60
        assert remainders[prefix] - line == base ** prefix * signed
        tail_area = sum((abs(d - 1) / base ** j
                         for j, d in enumerate(digits[prefix:], 1)), F())
        assert partial == tail_area / base ** prefix
        # The weighted-median estimate for arbitrary horizontal lines.
        best_line = F(digits[0]) / (base - 1)
        test_line = F(rng.randrange(-100, 201), 100)
        first_shell = 1 - 1 / base
        tail_shell = 1 / base
        assert first_shell * abs(test_line - best_line) - tail_shell * abs(test_line - best_line) == (1 - 2 / base) * abs(test_line - best_line)
        if not signed:
            assert all(d == 1 for d in digits)
            continue
        assert prefix < 60
        first = digits[prefix]
        assert (first == 0) == (signed < 0)
        constant = 2 * base - 1 if signed < 0 else base / (base - 2)
        lower, upper, _ = area_interval(x, x, base, 60)
        assert lower >= abs(signed) or upper >= abs(signed)
        assert upper <= constant * abs(signed)
    return samples


def pi_areas():
    plo, phi = pi_bounds()
    xlo, xhi = plo - 3, phi - 3
    bases = [F(8), F(121, 15), F(129, 16), F(37789, 4687), F(37918, 4703)]
    rows = []
    for base in bases:
        line = 1 / (base - 1)
        q = base.numerator - base.denominator
        p = 3 * q + base.denominator
        lower, upper, digits = area_interval(xlo, xhi, base)
        error_lo = min(abs(xlo - line), abs(xhi - line))
        error_hi = max(abs(xlo - line), abs(xhi - line))
        assert (xlo - line) * (xhi - line) > 0
        assert error_lo <= lower <= upper <= 17 * error_hi
        assert digits[0] == 1
        rows.append({
            'base': str(base), 'line': str(line), 'pi_fraction': f'{p}/{q}',
            'first_digits': digits[:14],
            'first_departure': next(j for j, digit in enumerate(digits, 1) if digit != 1),
            'q_squared_error': decimal_interval(q * q * error_lo, q * q * error_hi),
            'q_squared_area': decimal_interval(q * q * lower, q * q * upper),
            'area_over_error': decimal_interval(lower / error_hi, upper / error_lo),
        })
    return rows


def pattern_value(base, position, sign, pair):
    correction = (1 - 1 / base) if pair else F(1)
    return 1 / (base - 1) + sign * correction / base ** position


def pattern_root(position, sign, pair, xlo, xhi):
    lo, hi = F(7), F(9)
    assert pattern_value(lo, position, sign, pair) > xhi
    assert pattern_value(hi, position, sign, pair) < xlo
    for _ in range(200):
        mid = (lo + hi) / 2
        value = pattern_value(mid, position, sign, pair)
        if value > xhi:
            lo = mid
        elif value < xlo:
            hi = mid
        else:
            raise ArithmeticError('Root bisection reached the pi uncertainty')
    return lo, hi


def verify_pattern_prefix(blo, bhi, position, sign, pair, xlo, xhi):
    lo, hi = xlo, xhi
    for index in range(1, position + 13):
        product_lo, product_hi = blo * lo, bhi * hi
        digit = product_lo.__floor__()
        assert digit == product_hi.__floor__()
        expected = 1
        if index == position:
            expected += sign
        elif pair and index == position + 1:
            expected -= sign
        assert digit == expected, (position, sign, pair, index, digit, expected)
        lo, hi = product_lo - digit, product_hi - digit


def pattern_families():
    plo, phi = pi_bounds()
    xlo, xhi = plo - 3, phi - 3
    B_lo, B_hi = 1 + 1 / xhi, 1 + 1 / xlo
    rows = []
    for position in (2, 4, 8, 16):
        for sign in (-1, 1):
            for pair in (False, True):
                blo, bhi = pattern_root(position, sign, pair, xlo, xhi)
                verify_pattern_prefix(blo, bhi, position, sign, pair, xlo, xhi)
                if sign > 0:
                    gap_lo, gap_hi = blo - B_hi, bhi - B_lo
                else:
                    gap_lo, gap_hi = B_lo - bhi, B_hi - blo
                assert gap_lo > 0
                area_lo = ((1 + 1 / bhi) if pair else 1) / bhi ** position
                area_hi = ((1 + 1 / blo) if pair else 1) / blo ** position
                slope_lo, slope_hi = area_lo / gap_hi, area_hi / gap_lo
                expected_lo = xlo * xlo * ((1 + 2 * xlo) if pair else 1)
                expected_hi = xhi * xhi * ((1 + 2 * xhi) if pair else 1)
                if position == 16:
                    assert abs(slope_lo - expected_hi) < F(1, 10 ** 9)
                    assert abs(slope_hi - expected_lo) < F(1, 10 ** 9)
                rows.append({
                    'position': position, 'side': 'above' if sign > 0 else 'below',
                    'pattern': 'compensating_pair' if pair else 'single_defect',
                    'base_interval': decimal_interval(blo, bhi, 40),
                    'absolute_base_gap': decimal_interval(gap_lo, gap_hi),
                    'area_over_absolute_base_gap': decimal_interval(slope_lo, slope_hi),
                    'predicted_slope_limit': decimal_interval(expected_lo, expected_hi),
                })
    return rows


if __name__ == '__main__':
    result = {'exact_rational_samples': exact_finite_checks(),
              'pi_area_certificates': pi_areas(),
              'real_base_pattern_certificates': pattern_families()}
    print(json.dumps(result, ensure_ascii=False, indent=2))
