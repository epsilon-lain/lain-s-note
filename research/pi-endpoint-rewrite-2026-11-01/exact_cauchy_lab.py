#!/usr/bin/env python3
"""Exact checks for the grouped distinct-Taylor-degree sum.

These checks concern a finite combinatorial identity and bounds. The general
identity, collision estimate and endpoint obstruction are proved in the note;
they do not certify the pi paper or a bad-approximability conclusion.
"""
from fractions import Fraction
from math import factorial
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def ordered_distinct_coefficients(n, cutoff):
    # Coefficient of u^n in product_{d=0}^cutoff (1 + u q^d), then n!.
    dp = [[0] * (cutoff + 1) for _ in range(n + 1)]
    dp[0][0] = 1
    for d in range(cutoff + 1):
        for k in range(n, 0, -1):
            for total in range(cutoff, d - 1, -1):
                dp[k][total] += dp[k - 1][total - d]
    return [factorial(n) * x for x in dp[n]]


def predicted_coefficients(n, cutoff):
    # n! q^{n(n-1)/2} / product_{j=1}^n (1-q^j), as a formal series.
    shift = n * (n - 1) // 2
    p = [0] * (cutoff + 1)
    p[0] = 1
    for j in range(1, n + 1):
        for total in range(j, cutoff + 1):
            p[total] += p[total - j]
    return [0 if k < shift else factorial(n) * p[k - shift]
            for k in range(cutoff + 1)]


def exact_group_sum(n, gamma):
    product = Fraction(1)
    for j in range(1, n + 1):
        product *= 1 - gamma ** j
    return factorial(n) * gamma ** (n * (n - 1) // 2) / product


def truncated_group_sum(n, gamma, max_degree):
    # Coefficient of u^n in product_{d=0}^{max_degree}(1+u gamma^d).
    p = [Fraction(0)] * (n + 1)
    p[0] = Fraction(1)
    for d in range(max_degree + 1):
        for k in range(n, 0, -1):
            p[k] += p[k - 1] * gamma ** d
    return factorial(n) * p[n]


def serial(x):
    return str(x.numerator) if x.denominator == 1 else f'{x.numerator}/{x.denominator}'


def main():
    cutoff = 80
    formal = []
    for n in range(9):
        observed = ordered_distinct_coefficients(n, cutoff)
        predicted = predicted_coefficients(n, cutoff)
        assert observed == predicted
        formal.append({'n': n, 'coefficients_checked': cutoff + 1, 'passed': True})
    rational = []
    for gamma in (Fraction(1, 10), Fraction(1, 2), Fraction(3, 5)):
        for n in range(1, 7):
            full = exact_group_sum(n, gamma)
            partial = truncated_group_sum(n, gamma, 16)
            # Union bound for at least one degree >16, dropping distinctness.
            tail_bound = n * gamma ** 17 / (1 - gamma) ** n
            assert 0 <= full - partial <= tail_bound
            # Bound used to absorb denominators into an o_H(1) term.
            product = Fraction(1)
            for j in range(1, n + 1):
                product *= 1 - gamma ** j
            assert product >= (1 - gamma) ** n
            rational.append({'gamma': serial(gamma), 'n': n,
                             'full_sum': serial(full),
                             'finite_degree_sum': serial(partial),
                             'tail_union_bound': serial(tail_bound),
                             'exact_checks_passed': True})
    # Exact elementary bounds supporting the radius-optimized no-go estimates.
    # ln(2)>2/3 follows from 2*atanh(1/3)=2/3+positive terms.
    # e<3 is proved in the note from its positive factorial series.
    # sqrt(m+2)/(8 ln(2)(m+1)) < 3 sqrt(m+2)/(16(m+1)) < 1.
    endpoint = []
    for m in range(1, 501):
        assert 9 * (m + 2) < 256 * (m + 1) ** 2
        assert (m + 2) ** m < 3 * (m + 1) ** m
        endpoint.append({'m': m, 'integer_inequalities_passed': True})
    result = {
        'scope': 'Exact combinatorial and finite supporting checks; not a pi theorem certificate.',
        'formal_series_checks': formal,
        'rational_sum_checks': rational,
        'finite_supporting_inequalities': endpoint,
        'all_passed': True,
    }
    out = HERE / 'results_cauchy.json'
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'all_passed': True, 'formal_coefficients': 9 * (cutoff + 1),
                      'rational_cases': len(rational), 'supporting_dimensions': 500,
                      'output': str(out)}))


if __name__ == '__main__':
    main()
