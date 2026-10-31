"""Exact checks for rational-grid transport and chord-root certificates.

Run: python3 verify_rational_grid_chords.py
All decisions use integers/Fraction and independent Machin bounds for pi.
Finite experiments do not establish the P2848 liminf statement.
"""

import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    'carry_certificates', HERE / 'verify_rational_carry.py')
carry = importlib.util.module_from_spec(spec)
spec.loader.exec_module(carry)
cf = carry.cf
record = cf.layers.interval_record


def nearest_integer_distance(value):
    remainder = value - value.__floor__()
    return min(remainder, 1-remainder)


def unit_distance_interval(interval):
    low, high = interval
    assert 0 <= low <= high <= 1
    lower = min(low, 1-high)
    upper = (F(1, 2) if low <= F(1, 2) <= high
             else max(min(low, 1-low), min(high, 1-high)))
    return lower, upper


def difference_upper(first, second):
    return max(abs(first[0]-second[1]), abs(first[1]-second[0]))


def finite_rational_grid_certificate(n, digits):
    """No pi value or true-pi grid bracket is an input to this algorithm."""
    length, power = 3, 27
    while power < 6*n**3:
        length += 1
        power *= 3
    assert length <= len(digits)
    assert power >= 6*n**3
    assert length == 3 or power//3 < 6*n**3
    prefix = digits[:length]
    queries = 0

    def homogeneous_numerator(m):
        nonlocal queries
        queries += 1
        result, n_power = 1, 1
        for digit in prefix:
            n_power *= n
            result = m*result-digit*n_power
        return result

    low, high = 3*n, 4*n
    assert homogeneous_numerator(low) < 0 < homogeneous_numerator(high)
    while high-low > 1:
        middle = (low+high)//2
        sign = homogeneous_numerator(middle)
        assert sign != 0  # A rational root in (3,4) would be an integer.
        if sign < 0:
            low = middle
        else:
            high = middle
    assert high == low+1
    x, y = F(low, n), F(high, n)
    dx = F(homogeneous_numerator(low), low**length)
    dy = F(homogeneous_numerator(high), high**length)
    assert dx < 0 < dy
    assert dx == 1-cf.layers.phi_bounds(x, prefix)[0]
    assert dy == 1-cf.layers.phi_bounds(y, prefix)[0]
    assert n*(dy-dx) >= F(3, 16)
    root = x-dx/(n*(dy-dx))
    assert x < root < y
    score = n*nearest_integer_distance(n*root)
    return {'L': length, 'power_3_L': power, 'a_floor_n_beta_L': low,
            'bracket': [x, y], 'chord_root': root, 'J': score,
            'integer_sign_queries': queries}


def inspect_grid(n, digits, pi_interval, slope):
    pi_low, pi_high = pi_interval
    a = (n*pi_low).__floor__()
    assert a == (n*pi_high).__floor__()
    x, y = F(a, n), F(a+1, n)
    assert 3 <= x < pi_low < pi_high < y <= 4
    theta = (n*pi_low-a, n*pi_high-a)
    dist = unit_distance_interval(theta)
    kappa = tuple(n*value for value in dist)

    # U=F(x)-1>0, V=1-F(y)>0. The chord phase is U/(U+V).
    fx, fy = cf.layers.phi_bounds(x, digits), cf.layers.phi_bounds(y, digits)
    u, v = (fx[0]-1, fx[1]-1), (1-fy[1], 1-fy[0])
    assert u[0] > 0 and v[0] > 0
    tau = (u[0]/(u[0]+v[1]), u[1]/(u[1]+v[0]))
    assert tau[0] >= theta[1]
    phase_bound_lower = min(theta[0]*(1-theta[0]),
                            theta[1]*(1-theta[1]))/n
    assert tau[1]-theta[0] <= phase_bound_lower
    chord_root = (x+tau[0]/n, x+tau[1]/n)
    assert chord_root[0] >= pi_high
    assert chord_root[1]-pi_low <= kappa[0]/n**3
    chord_dist = unit_distance_interval(tau)
    chord_score = tuple(n*value for value in chord_dist)
    assert difference_upper(chord_score, kappa) <= kappa[0]/n

    gamma = (min(u[0], v[0]), min(u[1], v[1]))
    normalized_gamma = tuple(n*n*value for value in gamma)
    linear_term = (slope[0]*kappa[0], slope[1]*kappa[1])
    assert difference_upper(normalized_gamma, linear_term) <= kappa[0]**2/(4*n*n)

    # The prefix length is chosen by exact integer comparisons.
    length, power = 0, 1
    while power < 8*n**3:
        length += 1
        power *= 3
    assert 1 <= length <= len(digits)
    assert power >= 8*n**3 and power//3 < 8*n**3
    dx = 1-cf.layers.phi_bounds(x, digits[:length])[0]
    dy = 1-cf.layers.phi_bounds(y, digits[:length])[0]
    assert n*(dy-dx) >= F(3, 16)
    finite_root = x-dx/(n*(dy-dx))
    finite_score = n*nearest_integer_distance(n*finite_root)
    finite_tail_bound = F(8*n*n, power)
    assert finite_tail_bound <= F(1, n)
    assert difference_upper((finite_root, finite_root), chord_root) <= F(8, power)
    assert difference_upper((finite_score, finite_score), chord_score) <= finite_tail_bound
    total_bound_lower = kappa[0]/n + finite_tail_bound
    assert difference_upper((finite_score, finite_score), kappa) <= total_bound_lower

    # Primary certificate: bracket the truncated root using only integers.
    primary = finite_rational_grid_certificate(n, digits)
    primary_score = primary['J']
    primary_tail = F(6*n*n, primary['power_3_L'])
    assert primary_tail <= F(1, n)
    primary_error = difference_upper((primary_score, primary_score), kappa)
    assert primary_error <= kappa[0]/n+primary_tail

    return {
        'n': n,
        'grid_bracket': [str(x), str(y)],
        'a_floor_n_pi': a,
        'theta_fractional_part_n_pi': record(*theta),
        'distance_pi_to_grid': record(dist[0]/n, dist[1]/n),
        'kappa_n_times_distance_n_pi_to_Z': record(*kappa),
        'normalized_transport_minimum_n_squared_Gamma': record(*normalized_gamma),
        'chord_phase': record(*tau),
        'chord_root_minus_pi': record(chord_root[0]-pi_high, chord_root[1]-pi_low),
        'chord_score': record(*chord_score),
        'finite_prefix_length_L': length,
        'finite_score_rational': str(finite_score),
        'finite_score': record(finite_score, finite_score),
        'finite_score_error_upper': record(
            difference_upper((finite_score, finite_score), kappa),
            difference_upper((finite_score, finite_score), kappa)),
        'finite_score_error_bound': record(kappa[0]/n+finite_tail_bound,
                                           kappa[1]/n+finite_tail_bound),
        'primary_integer_only_certificate': {
            'L': primary['L'],
            'a_floor_n_beta_L': primary['a_floor_n_beta_L'],
            'bracket': list(map(str, primary['bracket'])),
            'same_bracket_as_pi': primary['a_floor_n_beta_L'] == a,
            'integer_sign_queries': primary['integer_sign_queries'],
            'J_exact_rational': str(primary_score),
            'J': record(primary_score, primary_score),
            'absolute_error_upper': record(primary_error, primary_error),
            'bound_kappa_over_n_plus_6n_squared_3_to_minus_L':
                record(kappa[0]/n+primary_tail, kappa[1]/n+primary_tail),
        },
        'all_exact_checks_passed': True,
    }


def main():
    pi_low, pi_high, digits = cf.layers.transport.certify_digits(192)
    partials, status = cf.certified_cf(pi_low, pi_high, 20)
    assert len(partials) == 20
    cf_denominators = [value.denominator for value in cf.convergents(partials)]
    all_n = sorted(set(range(1, 257)) | set(cf_denominators) | {33102})
    m_low = cf.moment_bounds(pi_high, digits)[0]
    m_high = cf.moment_bounds(pi_low, digits)[1]
    slope = (m_low/pi_high, m_high/pi_low)
    rows = [inspect_grid(n, digits, (pi_low, pi_high), slope) for n in all_n]
    selected = {1, 7, 106, 113, 33102}
    results = {
        'status': 'All finite exact checks passed; no proof of the P2848 liminf.',
        'pi_certificate': 'Machin rational enclosure, 600-bit outward grid',
        'native_digits_certified': len(digits),
        'native_digit_prefix': ''.join(map(str, digits[:64])),
        'slope_c_m_over_pi': record(*slope),
        'grid_sample': {
            'full_range': [1, 256],
            'extra_cf_denominators_from_20_terms': cf_denominators,
            'cf_status': status,
            'distinct_n_count': len(all_n),
            'maximum_n': max(all_n),
            'maximum_finite_prefix_length': max(row['finite_prefix_length_L'] for row in rows),
            'maximum_primary_prefix_length': max(row['primary_integer_only_certificate']['L'] for row in rows),
            'primary_integer_sign_queries': sum(row['primary_integer_only_certificate']['integer_sign_queries'] for row in rows),
            'different_truncated_and_pi_grid_brackets': sum(not row['primary_integer_only_certificate']['same_bracket_as_pi'] for row in rows),
        },
        'checks_for_every_n': [
            'floor(n*pi) certified independently of the transported values',
            'F(x)>1>F(y)',
            '0 <= tau-theta <= theta*(1-theta)/n',
            '0 <= chord_root-pi <= kappa/n^3',
            '|chord_score-kappa| <= kappa/n',
            '|n^2*Gamma-c*kappa| <= kappa^2/(4*n^2)',
            'minimal L satisfies 3^L >= 8*n^3',
            'finite interpolation slope >= 3/16',
            '|finite_root-chord_root| <= 8*3^(-L)',
            '|finite_score-chord_score| <= 8*n^2*3^(-L)',
            '|finite_score-kappa| <= kappa/n+8*n^2*3^(-L)',
            'primary certificate obtains its own bracket by integer sign binary search, with no pi-value input',
            'primary minimal L>=3 satisfies 3^L>=6*n^3',
            '|primary_J-kappa| <= kappa/n+6*n^2*3^(-L) <= (kappa+1)/n',
        ],
        'selected_rows': [row for row in rows if row['n'] in selected],
        'limits': [
            'No claim that all infinite transport values at rational bases are rational.',
            'The primary algorithm uses only n and L native digits and brackets the finite root beta_L, not pi.',
            'The separate true-pi-bracket diagnostic uses an independent pi enclosure and is not an input to the primary algorithm.',
            'Grid n is retained even when its endpoint fractions reduce.',
        ],
    }
    destination = HERE/'rational_grid_chord_results.json'
    destination.write_text(json.dumps(results, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'distinct_n': len(all_n),
                      'max_n': max(all_n), 'max_L': results['grid_sample']['maximum_finite_prefix_length'],
                      'output': str(destination)}, indent=2))


if __name__ == '__main__':
    main()
