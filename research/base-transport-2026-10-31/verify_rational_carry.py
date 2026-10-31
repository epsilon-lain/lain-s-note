"""Exact integer certificates for rational-base carry survival.

Run: python3 verify_rational_carry.py
All decisions use integers and Fraction; decimal strings are outward displays.
The accompanying mathematical note proves the infinite statements.
"""

import importlib.util
import json
from fractions import Fraction as F
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    'cf_migration', HERE/'verify_continued_fraction_migration.py')
cf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cf)


def inspect_base(base, digits, pi_interval, slope_interval, test_taylor):
    base = F(base)
    p, q = base.numerator, base.denominator
    assert 3 <= base <= 4 and gcd(p, q) == 1
    critical = 0
    while p**critical < q**(critical+2):
        critical += 1
    assert p**critical >= q**(critical+2)
    assert critical == 0 or p**(critical-1) < q**(critical+1)

    numerator = 1
    numerators = [numerator]
    escape = None
    suffix_checks = head_checks = gcd_checks = truncation_checks = 0
    last_nonzero = 0
    for n, digit in enumerate(digits, 1):
        numerator = p*numerator-digit*q**n
        numerators.append(numerator)
        assert gcd(numerator, q) == 1
        assert F(numerator, q**n).denominator == q**n
        gcd_checks += 1
        if digit:
            last_nonzero = n
        if gcd(p, 6) == 1:
            assert F(numerator, p**n).denominator == p**last_nonzero
            truncation_checks += 1
        if p >= 4:
            for k in range(1, min(n, 5)+1):
                modulus = p**k
                q_inverse = pow(q, -1, modulus)
                actual = numerator*pow(q_inverse, n, modulus) % modulus
                expected = -sum(digits[n-1-i]*pow(p*q_inverse, i, modulus)
                                for i in range(k)) % modulus
                assert actual == expected
                suffix_checks += 1
        if q > 1:
            for k in range(1, min(n+1, 5)+1):
                modulus = q**k
                p_inverse = pow(p, -1, modulus)
                actual = numerator*pow(p_inverse, n, modulus) % modulus
                expected = (1-sum(digits[j-1]*pow(q*p_inverse, j, modulus)
                                  for j in range(1, k))) % modulus
                assert actual == expected
                head_checks += 1
        if escape is None and abs(numerator) >= 2*q**n:
            escape = n
        if escape is not None:
            assert abs(numerator) >= 2*q**n
            assert (numerator > 0) == (numerators[escape] > 0)
            if n == escape+3:
                break
    assert escape is not None, 'Increase certified prefix for this base'
    assert len(numerators) == escape+4
    assert all(abs(value) < 2*q**n for n, value in enumerate(numerators[:escape]))

    phi_low, phi_high = cf.layers.phi_bounds(base, digits)
    displacement = (1-phi_high, 1-phi_low)
    magnitude = cf.absolute_interval(*displacement)
    error = cf.absolute_interval(base-pi_interval[1], base-pi_interval[0])
    assert magnitude[0] > 0 and error[0] > 0
    assert (displacement[0] > 0) == (numerators[escape] > 0)
    assert (base > pi_interval[1]) == (numerators[escape] > 0)
    assert F(1, 2)*base**(-escape) <= magnitude[0]
    assert magnitude[1] < 14*base**(-escape)
    kappa = (q*q*error[0], q*q*error[1])
    excess = escape-critical
    assert F(1, 4)*base**(-excess) <= kappa[0]
    assert kappa[1] < 56*base**(-excess)

    survival_n = escape-1
    survival_numerator = numerators[survival_n]
    height_ratio = F(p**survival_n, q**(survival_n+2))
    max_integer_M = height_ratio.__floor__()
    assert abs(survival_numerator) < 2*q**survival_n
    if max_integer_M >= 1:
        assert p**survival_n >= max_integer_M*q**(survival_n+2)
        assert kappa[1] < F(14, max_integer_M)

    taylor = None
    if test_taylor:
        defect = (q*q*magnitude[0]-slope_interval[1]*kappa[1],
                  q*q*magnitude[1]-slope_interval[0]*kappa[0])
        rigorous_upper = max(abs(defect[0]), abs(defect[1]))
        permitted_lower = kappa[0]**2/F(4*q*q)
        assert rigorous_upper <= permitted_lower
        taylor = {'verified': True,
                  'absolute_defect_upper': cf.layers.interval_record(rigorous_upper, rigorous_upper),
                  'bound_formula': 'kappa^2/(4*q^2)',
                  'permitted_bound_lower': cf.layers.interval_record(permitted_lower, permitted_lower)}

    return {
        'base': str(base), 'p': p, 'q': q,
        'escape_depth_E': escape, 'critical_depth_B': critical,
        'excess_depth_H': excess,
        'escape_sign': 'positive' if numerators[escape] > 0 else 'negative',
        'last_surviving_K': str(F(survival_numerator, q**survival_n)),
        'escape_K': str(F(numerators[escape], q**escape)),
        'displacement_1_minus_F': cf.layers.interval_record(*displacement),
        'kappa_q_squared_abs_base_minus_pi': cf.layers.interval_record(*kappa),
        'normalized_displacement_q_squared_abs_1_minus_F':
            cf.layers.interval_record(q*q*magnitude[0], q*q*magnitude[1]),
        'integer_strip_certificate': {
            'N': survival_n, 'A_N': str(survival_numerator),
            'max_integer_M': max_integer_M,
            'height_ratio_p_to_N_over_q_to_N_plus_2': str(height_ratio),
            'all_required_inequalities_verified': max_integer_M >= 1},
        'checks': {'coprime_q_and_exact_carry_denominator': gcd_checks,
                   'suffix_mod_p_power': suffix_checks,
                   'prefix_mod_q_power': head_checks,
                   'truncation_denominator_when_gcd_p_6_is_1': truncation_checks,
                   'absorbing_steps_after_escape': 3},
        'optional_taylor_refinement': taylor}


def main():
    pi_low, pi_high, digits = cf.layers.transport.certify_digits(192)
    partials, status = cf.certified_cf(pi_low, pi_high, 40)
    assert len(partials) == 40
    approximants = cf.convergents(partials)
    extra_bases = [F(4), F(19, 6), F(25, 8)]
    bases = list(dict.fromkeys(approximants+extra_bases))
    m_low = cf.moment_bounds(pi_high, digits)[0]
    m_high = cf.moment_bounds(pi_low, digits)[1]
    slope = (m_low/pi_high, m_high/pi_low)
    rows = [inspect_base(base, digits, (pi_low, pi_high), slope,
                         test_taylor=(index < 20 or base in extra_bases))
            for index, base in enumerate(bases)]
    result = {
        'native_digit_certificate': {'count': len(digits),
            'method': 'Machin rational bounds, 600-bit outward enclosure'},
        'pi_cf_certificate': {'partial_quotients': partials, 'status': status},
        'definitions': {
            'A_0': 1, 'A_N': 'p*A_(N-1)-d_N*q^N',
            'E': 'min {N>=1: |A_N|>=2*q^N}',
            'B': 'min {N>=0: p^N>=q^(N+2)}', 'H': 'E-B'},
        'universal_bounds_finitely_checked': [
            '(1/2)*t^(-E) <= |1-F(t)| < 14*t^(-E)',
            '(1/4)*t^(-H) <= q^2*|t-pi| < 56*t^(-H)'],
        'counts': {'pi_convergents': 40, 'control_bases': 3,
                   'total_bases': len(rows),
                   'largest_escape_depth': max(row['escape_depth_E'] for row in rows),
                   'taylor_refinements': sum(row['optional_taylor_refinement'] is not None for row in rows)},
        'records': rows,
        'scope': 'Exact finite checks and integer certificates; accompanying proofs justify universal bounds. No conclusion on P2848.'}
    (HERE/'rational_carry_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    wanted = {'22/7', '333/106', '355/113', '103993/33102'}
    print(json.dumps({'counts': result['counts'], 'selected': [
        {key: row[key] for key in ['base', 'escape_depth_E', 'critical_depth_B',
                                   'excess_depth_H', 'escape_sign']}
        | {'M': row['integer_strip_certificate']['max_integer_M']}
        for row in rows if row['base'] in wanted]}, indent=2))


if __name__ == '__main__':
    main()
