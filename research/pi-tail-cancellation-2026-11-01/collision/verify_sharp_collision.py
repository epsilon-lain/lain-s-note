#!/usr/bin/env python3
"""Exact checks for the sharp distinct-Taylor-degree summation estimate."""
from fractions import Fraction as F
from math import factorial, prod
from pathlib import Path
import json


def finite_elementary(r, nmax, last_degree):
    values = [F(1)] + [F(0)] * nmax
    for d in range(last_degree + 1):
        z = r ** d
        for n in range(nmax, 0, -1):
            values[n] += z * values[n - 1]
    return values


def compositions(n):
    if n == 0:
        yield ()
    else:
        for first in range(1, n + 1):
            for tail in compositions(n - first):
                yield (first,) + tail


def polynomial_division(a, b):
    a = list(map(F, a))
    b = list(map(F, b))
    q = [F(0)] * (len(a) - len(b) + 1)
    for k in range(len(q) - 1, -1, -1):
        q[k] = a[k + len(b) - 1] / b[-1]
        for j, coefficient in enumerate(b):
            a[k + j] -= q[k] * coefficient
    assert all(x == 0 for x in a[len(b) - 1:])
    return q, a[:len(b) - 1]


def main():
    max_n, max_degree = 8, 70
    # Coefficients of the elementary symmetric functions of 1,x,...,x^70.
    # Degrees above 70 cannot affect any compared coefficient.
    coefficients = [[0] * (max_degree + 1) for _ in range(max_n + 1)]
    coefficients[0][0] = 1
    for j in range(max_degree + 1):
        for n in range(max_n, 0, -1):
            for d in range(j, max_degree + 1):
                coefficients[n][d] += coefficients[n - 1][d - j]
    checked_coefficients = 0
    for n in range(max_n + 1):
        partition = [1] + [0] * max_degree
        for k in range(1, n + 1):
            for d in range(k, max_degree + 1):
                partition[d] += partition[d - k]
        shift = n * (n - 1) // 2
        for d in range(max_degree + 1):
            rhs = 0 if d < shift else factorial(n) * partition[d - shift]
            lhs = factorial(n) * coefficients[n][d]
            assert lhs == rhs, (n, d, lhs, rhs)
            checked_coefficients += 1

    checked_truncations = 0
    rational_radii = [F(1, 2), F(7, 100), F(2, 3)]
    for r in rational_radii:
        for last in [8, 20, 40]:
            finite = finite_elementary(r, max_n, last)
            for n in range(max_n + 1):
                closed = F(factorial(n)) * r ** (n * (n - 1) // 2)
                closed /= prod((1 - r ** k for k in range(1, n + 1)), start=F(1))
                partial = factorial(n) * finite[n]
                gap = closed - partial
                # Union bound: at least one ordered coordinate is >last;
                # ignore distinctness for all coordinates in this upper bound.
                tail = F(0) if n == 0 else n * r ** (last + 1) / (1 - r) ** n
                assert F(0) <= gap <= tail, (r, n, last)
                checked_truncations += 1

    checked_compositions = 0
    for total in range(1, 13):
        for parts in compositions(total):
            assert prod(factorial(n) for n in parts) <= factorial(total)
            for r in rational_radii:
                denominator = prod(
                    (1 - r ** k for n in parts for k in range(1, n + 1)),
                    start=F(1),
                )
                assert denominator >= (1 - r) ** total
            assert sum(n * (n - 1) for n in parts) == sum(n * n for n in parts) - total
            checked_compositions += 1

    # x^4(1-x)^4/(1+x^2): polynomial division proves the positive
    # integral equals 22/7 - 4 atan(1) = 22/7 - pi.
    numerator = [0, 0, 0, 0, 1, -4, 6, -4, 1]
    quotient, remainder = polynomial_division(numerator, [1, 0, 1])
    assert remainder == [-4, 0]
    integral_polynomial = sum((a / (i + 1) for i, a in enumerate(quotient)), F(0))
    assert integral_polynomial == F(22, 7)
    # ln2 = 2 sum_{k>=0} (1/3)^(2k+1)/(2k+1)
    # < 2/3 + (2/3) sum_{k>=1}(1/3)^(2k+1) = 25/36.
    log2_bound = F(2, 3) + F(2, 3) * F(1, 27) / (1 - F(1, 9))
    assert log2_bound == F(25, 36)
    for centers in range(1, 101):
        assert 2 * F(22, 7) * (centers - 1) + log2_bound < 7 * centers
        assert F(7 * centers, 14 * centers) == F(1, 2)

    # The radius objective is g(u)=d ln u-Au, u>1. Its derivative
    # and stationary point are rational for rational d,A,u.
    derivative_checks = 0
    for d, A in [(F(1), F(1)), (F(1, 2), F(3)),
                 (F(3), F(1)), (F(7, 5), F(2, 7))]:
        if d <= A:
            for u in [F(1001, 1000), F(3, 2), F(2), F(10)]:
                assert d / u - A < 0
                assert -d / u ** 2 < 0
                derivative_checks += 1
        else:
            stationary = d / A
            assert stationary > 1
            assert d / stationary - A == 0
            assert d / ((1 + stationary) / 2) - A > 0
            assert d / (2 * stationary) - A < 0
            assert -d / stationary ** 2 < 0
            derivative_checks += 4

    family_checks = []
    for m in [8, 10, 16, 24, 32, 40, 64, 80, 100]:
        theta, B, C = F(1, 2), F(3, 5), F(13, 10)
        power = C ** m
        K = power.numerator // power.denominator
        w0 = B ** (-m)
        v0 = 2 * K * theta ** m * w0
        d, A = 1 / (2 * v0), 7 * K / w0
        assert power >= 8 and K >= F(7, 8) * power
        assert v0 >= F(7, 4) * F(13, 12) ** m
        assert d / A == 1 / (28 * K ** 2 * theta ** m)
        assert d / A <= F(16, 343) * F(200, 169) ** m
        assert A <= 7 * F(39, 50) ** m
        # This rational bound follows from log(1+x)<=x.
        upper_gain = F(62 * m, 1183) * F(12, 13) ** m
        assert d * F(31 * m, 169) <= upper_gain
        family_checks.append({'m': m, 'K': K, 'stationary_ratio_d_over_A': str(d/A),
                              'gain_upper_bound': str(upper_gain)})

    result = {
        'status': 'passed',
        'integer_generating_coefficients_checked': checked_coefficients,
        'exact_rational_truncations_with_tail_bound': checked_truncations,
        'positive_integer_compositions_checked': checked_compositions,
        'rational_radius_ratios': [str(r) for r in rational_radii],
        'pi_upper_bound_positive_integral': 'integral_0^1 x^4(1-x)^4/(1+x^2) dx = 22/7 - pi > 0',
        'log_two_rational_upper': str(log2_bound),
        'disk_bound_centers_checked': 100,
        'old_collision_coefficient': 'ln(2)/4',
        'new_collision_coefficient_at_R_14K': 'ln(2)/2',
        'new_holomorphic_error': 'max(14K/w0, ln(28K)/w_min) + ln(2)/v0',
        'radius_objective_derivative_checks': derivative_checks,
        'radius_optimized_family_checks': family_checks,
        'radius_gain_uniform_bound': '-7(39/50)^m <= G_m <= (62m/1183)(12/13)^m',
        'scope': 'Exact combinatorial checks, not a validation of manuscript interpolation thresholds or a result on pi bad approximability.'
    }
    output = Path(__file__).with_name('sharp_collision_verification.json')
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
