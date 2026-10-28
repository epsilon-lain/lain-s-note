"""Exact regression coefficients of the normalized beta-digit coordinate cloud.

Standard library only; uses a sibling certified Machin enclosure for pi.
No neural network or asymptotic theorem is inferred from these finite checks.
"""

from fractions import Fraction as F
import json
import random

from rational_base_gamma_checks import pi_bounds, decimal_interval


def interval_product(first, second):
    values = [a * b for a in first for b in second]
    return min(values), max(values)


def interval_ratio(first, second):
    assert second[0] * second[1] > 0
    return interval_product(first, (1 / second[1], 1 / second[0]))


def fit_interval(base, xlo, xhi, length=70):
    lo, hi = xlo, xhi
    squared_weight_sum, height_square_sum = F(0), F(0)
    weight = F(1)
    digits = []
    D = -(-base.numerator // base.denominator) - 1
    for _ in range(length):
        product_lo, product_hi = base * lo, base * hi
        digit = product_lo.__floor__()
        assert digit == product_hi.__floor__(), 'More source precision needed'
        lo, hi = product_lo - digit, product_hi - digit
        digits.append(digit)
        weight /= base
        squared_weight_sum += digit * weight * weight
        height_square_sum += F(digit * digit) * weight / (base - 1)
    compressed_tail = min(F(D) * weight * weight / (base * base - 1),
                          weight * weight * hi / base)
    moment_square_tail = F(D) * weight * hi / (base - 1)
    slo = 6 * (base + 1) * squared_weight_sum - 6 * xhi
    shi = 6 * (base + 1) * (squared_weight_sum + compressed_tail) - 6 * xlo
    ilo, ihi = xlo - shi / 2, xhi - slo / 2
    sq_min = F(0) if slo <= 0 <= shi else min(slo * slo, shi * shi)
    sq_max = max(slo * slo, shi * shi)
    mse_lo = max(F(0), height_square_sum - xhi * xhi - sq_max / 12)
    mse_hi = height_square_sum + moment_square_tail - xlo * xlo - sq_min / 12
    assert mse_hi >= mse_lo
    correction_lo = 6 * (base + 1) * (squared_weight_sum - 1 / (base * base - 1))
    correction_hi = correction_lo + 6 * (base + 1) * compressed_tail
    return {'slope': (slo, shi), 'intercept': (ilo, ihi),
            'mse': (mse_lo, mse_hi), 'digits': digits,
            'quadratic_correction': (correction_lo, correction_hi)}


def finite_moment_checks():
    rng = random.Random(2026102817)
    count = 1000
    for _ in range(count):
        denominator = rng.randrange(1, 50)
        base = 8 + F(rng.randrange(denominator + 1), denominator)
        D = -(-base.numerator // base.denominator) - 1
        digits = [rng.randrange(D + 1) for _ in range(rng.randrange(1, 14))]
        mean, cross_moment, second_moment, compressed = F(0), F(0), F(0), F(0)
        for j, d in enumerate(digits, 1):
            left, right = base ** (-j), base ** (1 - j)
            height = F(d) / (base - 1)
            mean += height * (right - left)
            cross_moment += height * (right * right - left * left) / 2
            second_moment += height * height * (right - left)
            compressed += d * base ** (-2 * j)
        slope = 12 * cross_moment - 6 * mean
        intercept = mean - slope / 2
        assert slope == 6 * (base + 1) * compressed - 6 * mean
        assert mean - slope / 2 - intercept == 0
        assert cross_moment - slope / 3 - intercept / 2 == 0
        mse = second_moment - mean * mean - slope * slope / 12
        assert mse >= 0
        # Exact midpoint-cloud variance, after adding within-shell variance.
        midpoint_variance = base / (4 * (base * base + base + 1))
        within_variance = (base - 1) ** 2 / (12 * (base * base + base + 1))
        assert midpoint_variance + within_variance == F(1, 12)
    return count


def pi_certificate(base, xlo, xhi):
    fit = fit_interval(base, xlo, xhi)
    slo, shi = fit['slope']
    error = (xlo - 1 / (base - 1), xhi - 1 / (base - 1))
    assert error[0] * error[1] > 0
    abs_error = min(abs(error[0]), abs(error[1])), max(abs(error[0]), abs(error[1]))
    if error[0] > 0:
        assert shi < 0
        abs_slope = -shi, -slo
    else:
        assert slo > 0
        abs_slope = slo, shi
        assert fit['digits'][:2] == [1, 1]
    assert fit['digits'][0] == 1
    assert abs_slope[0] > 4 * abs_error[1]
    assert abs_slope[1] < 8 * abs_error[0]
    assert max(abs(v) for v in fit['quadratic_correction']) < 8160 * abs_error[0] ** 2
    q = base.numerator - base.denominator
    p = 3 * base.numerator - 2 * base.denominator
    return {
        'base': str(base), 'pi_fraction': f'{p}/{q}',
        'slope_interval': decimal_interval(*fit['slope']),
        'intercept_interval': decimal_interval(*fit['intercept']),
        'mse_interval': decimal_interval(*fit['mse']),
        'slope_over_signed_error': decimal_interval(*interval_ratio(fit['slope'], error)),
        'normalized_slope_q2_abs_s_over_6': decimal_interval(q * q * abs_slope[0] / 6,
                                                           q * q * abs_slope[1] / 6),
        'normalized_error_q2': decimal_interval(q * q * abs_error[0], q * q * abs_error[1]),
        'first_digits': fit['digits'][:14],
    }


def second_level_fits(xlo, xhi):
    center = F(8062513305931, 10 ** 12)
    rows = []
    for width in (F(1, 100), F(1, 1000), F(1, 10000), F(1, 1000000)):
        bases = [center + width * F(i, 10) for i in range(-10, 11)]
        values = [fit_interval(base, xlo, xhi)['slope'] for base in bases]
        mean_base = sum(bases, F()) / len(bases)
        denominator = sum(((base - mean_base) ** 2 for base in bases), F())
        low, high = F(0), F(0)
        for base, interval in zip(bases, values):
            weight = (base - mean_base) / denominator
            contribution = interval_product((weight, weight), interval)
            low += contribution[0]
            high += contribution[1]
        rows.append({'half_width': str(width), 'sample_count': len(bases),
                     'fitted_slope_of_base_to_slope': decimal_interval(low, high),
                     'predicted_derivative_at_special_base': decimal_interval(-6 * xhi * xhi,
                                                                             -6 * xlo * xlo)})
    return rows


if __name__ == '__main__':
    plo, phi = pi_bounds()
    xlo, xhi = plo - 3, phi - 3
    selected = [F(8), F(121, 15), F(129, 16), F(37789, 4687),
                F(37918, 4703), F(89337364, 11080585), F(9)]
    sample_certificates = [pi_certificate(base, xlo, xhi) for base in selected]
    grid_count = 0
    for i in range(101):
        pi_certificate(8 + F(i, 100), xlo, xhi)
        grid_count += 1
    output = {'finite_moment_checks': finite_moment_checks(),
              'pi_grid_certificates': grid_count,
              'selected_bases': sample_certificates,
              'base_to_slope_linear_fits': second_level_fits(xlo, xhi)}
    print(json.dumps(output, ensure_ascii=False, indent=2))
