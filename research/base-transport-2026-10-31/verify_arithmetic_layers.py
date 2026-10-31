"""Exact finite checks for arithmetic layers in frozen-digit base transport.

Run: python3 verify_arithmetic_layers.py
Reuses certified native digits from verify_transport.py. Integer/Fraction
arithmetic controls every decision; outward decimals are display only.
The inverse asymptotic estimate and infinite counting theorem require proofs.
"""

import importlib.util
import json
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('transport_certificate', HERE / 'verify_transport.py')
transport = importlib.util.module_from_spec(spec)
spec.loader.exec_module(transport)


def outward(value, rounding):
    with localcontext() as context:
        context.prec = 40
        context.rounding = rounding
        return str(Decimal(value.numerator)/Decimal(value.denominator))


def interval_record(low, high):
    assert low <= high
    return {'lower': outward(low, ROUND_FLOOR),
            'upper': outward(high, ROUND_CEILING)}


def phi_bounds(base, digits):
    """Evaluate the polynomial using integer homogeneous coordinates.

    For base=a/b, the partial sum has denominator a**N. Its integer
    numerator is built without repeated rational simplification. The unknown
    tail is bounded by 3*(b/a)**(N+1)/(1-b/a).
    """
    base = F(base)
    numerator, denominator = base.numerator, base.denominator
    assert numerator > denominator > 0
    partial_numerator = 0
    denominator_power = 1
    for digit in digits:
        denominator_power *= denominator
        partial_numerator = numerator*partial_numerator + digit*denominator_power
    numerator_power = numerator**len(digits)
    low = F(partial_numerator, numerator_power)
    high = F(partial_numerator*(numerator-denominator)
             + 3*denominator_power*denominator,
             numerator_power*(numerator-denominator))
    return low, high


def inverse_interval(target, left, right, digits, requested_width=F(1, 2**80)):
    """Certified bisection with a derivative fallback for ambiguous probes.

    Phi is strictly decreasing, and -Phi'(t)>=3/right**2 throughout a
    bracket, since the first digit equals 3. If a probe interval includes
    the target, this derivative lower bound gives a valid root enclosure;
    no assumption about strict inequality at the midpoint is needed.
    """
    left, right, target = F(left), F(right), F(target)
    assert digits[0] == 3
    assert phi_bounds(left, digits)[0] > target
    assert phi_bounds(right, digits)[1] < target
    bisections = fallback_steps = 0
    while right-left > requested_width:
        middle = (left+right)/2
        value_low, value_high = phi_bounds(middle, digits)
        if value_low > target:
            left = middle
            bisections += 1
        elif value_high < target:
            right = middle
            bisections += 1
        else:
            derivative_lower = F(3)/right**2
            radius = max(target-value_low, value_high-target)/derivative_lower
            new_left = max(left, middle-radius)
            new_right = min(right, middle+radius)
            assert new_right-new_left < right-left, 'Increase certified digit count'
            left, right = new_left, new_right
            fallback_steps += 1
    assert 0 < right-left <= requested_width
    value_low = phi_bounds(right, digits)[0]
    value_high = phi_bounds(left, digits)[1]
    assert value_low <= target <= value_high
    return (left, right), (value_low, value_high), {
        'bisections': bisections, 'derivative_fallback_steps': fallback_steps,
        'certified_width_upper': str(requested_width),
        'actual_width': str(right-left)}


def arithmetic_layers(digits, pi_interval):
    results = []
    for denominator in [2, 3, 10, 100, 1000]:
        target = 1-F(1, denominator)
        root, value, algorithm = inverse_interval(target, 3, 10, digits)
        left, right = root
        assert left > pi_interval[1]
        low_value, high_value = value
        assert high_value < 1
        inverse_D = (1/(1-low_value), 1/(1-high_value))
        assert inverse_D[0] <= denominator <= inverse_D[1]
        results.append({
            'Q': denominator, 'definition': 'Phi(t_Q)=1-1/Q; C(1/t_Q)=Q',
            'target_Phi': str(target),
            't_Q': interval_record(left, right),
            't_Q_exact_rational_interval': [str(left), str(right)],
            't_Q_minus_pi': interval_record(left-pi_interval[1], right-pi_interval[0]),
            'Phi_t_Q_minus_target': interval_record(low_value-target, high_value-target),
            'C_1_over_t_Q_minus_Q': interval_record(inverse_D[0]-denominator,
                                                  inverse_D[1]-denominator),
            'algorithm': algorithm
        })
    return results


def inverse_asymptotics(digits):
    results = []
    for denominator in [3, 10, 100]:
        root, value, algorithm = inverse_interval(F(1, denominator),
                                                 3*denominator, 4*denominator, digits)
        approximation = (F(3*denominator) + F(1, 9*denominator)
                         + F(1, 27*denominator**2) - F(2, 243*denominator**3))
        error = (root[0]-approximation, root[1]-approximation)
        proposed_bound = F(24, denominator**4)
        assert -proposed_bound <= error[0] <= error[1] <= proposed_bound
        results.append({
            'q': denominator, 'definition': 'tau(1/q)=Phi^(-1)(1/q)',
            'root': interval_record(*root),
            'root_exact_rational_interval': [str(root[0]), str(root[1])],
            'approximation_exact': str(approximation),
            'approximation_formula': '3*q+1/(9*q)+1/(27*q^2)-2/(243*q^3)',
            'error': interval_record(*error),
            'absolute_error_bound_tested': str(proposed_bound),
            'within_24_over_q4': True, 'algorithm': algorithm
        })
    return results


def reduced_rational_counts(digits, pi_interval):
    lower_endpoint = phi_bounds(F(10), digits)
    upper_endpoint = phi_bounds(F(4), digits)
    target_counts = [10, 30, 100, 300, 1000]
    count = 0
    output = []
    for denominator in range(1, max(target_counts)+1):
        lower_numerator = (denominator*lower_endpoint[0]).__ceil__()
        assert lower_numerator == (denominator*lower_endpoint[1]).__ceil__()
        upper_numerator = (denominator*upper_endpoint[0]).__floor__()
        assert upper_numerator == (denominator*upper_endpoint[1]).__floor__()
        count += sum(gcd(numerator, denominator) == 1
                     for numerator in range(lower_numerator, upper_numerator+1))
        if denominator in target_counts:
            quotient = F(count, denominator**2)
            output.append({'Q': denominator, 'count': count,
                           'count_over_Q_squared_exact': str(quotient),
                           'count_over_Q_squared': interval_record(quotient, quotient)})
    expected = (F(3)/pi_interval[1]**2*(upper_endpoint[0]-lower_endpoint[1]),
                F(3)/pi_interval[0]**2*(upper_endpoint[1]-lower_endpoint[0]))
    return {
        'base_interval_J': '[4,10]',
        'rational_value_interval': '[Phi(10),Phi(4)]',
        'Phi_10': interval_record(*lower_endpoint),
        'Phi_4': interval_record(*upper_endpoint),
        'endpoint_decisions_certified_for_every_denominator_through': max(target_counts),
        'counts': output,
        'asymptotic_coefficient_formula': '(3/pi^2)*(Phi(4)-Phi(10))',
        'asymptotic_coefficient': interval_record(*expected)
    }


def main():
    degree = 192
    pi_low, pi_high, digits = transport.certify_digits(degree)
    # Deliberately place a rational target inside the midpoint's enclosure,
    # so the derivative fallback is actually exercised rather than only read.
    midpoint_value = phi_bounds(F(4), digits)
    fallback_target = sum(midpoint_value)/2
    fallback_root, _, fallback_algorithm = inverse_interval(
        fallback_target, 3, 5, digits)
    assert fallback_algorithm['derivative_fallback_steps'] >= 1
    result = {
        'native_digit_certificate': {
            'count': degree, 'pi_source': 'Machin rational bounds, 160 alternating terms',
            'working_intervals': '600-bit dyadic outward enclosures',
            'all_decisions': 'integer and Fraction arithmetic',
            'display': '40-significant-digit Decimal outward bounds'},
        'arithmetic_layers_near_pi': arithmetic_layers(digits, (pi_low, pi_high)),
        'inverse_asymptotic_finite_checks': inverse_asymptotics(digits),
        'ambiguous_midpoint_fallback_check': {
            'purpose': 'A rational target inside the Phi(4) enclosure exercises the derivative branch.',
            'enclosed_root': interval_record(*fallback_root),
            'algorithm': fallback_algorithm},
        'reduced_rational_layer_counts': reduced_rational_counts(digits, (pi_low, pi_high)),
        'scope': 'Finite certificates only; the inverse expansion bound for all q and the infinite asymptotic counting law require the accompanying proofs. No P2848 conclusion.'
    }
    (HERE/'arithmetic_layer_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({
        'arithmetic_layers': len(result['arithmetic_layers_near_pi']),
        'inverse_asymptotic_cases': len(result['inverse_asymptotic_finite_checks']),
        'rational_count_max_denominator': 1000,
        'output': 'arithmetic_layer_results.json'
    }, indent=2))


if __name__ == '__main__':
    main()
