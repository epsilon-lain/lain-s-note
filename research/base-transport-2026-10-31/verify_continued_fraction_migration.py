"""Exact finite certificates for continued fractions of Psi(t)=t*Phi(t).

Run: python3 verify_continued_fraction_migration.py
All certification decisions use Fraction. The optional SVG is a floating-point
display of certified intervals and series, not an independent proof.
"""

import importlib.util
import json
import math
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('arithmetic_layers', HERE/'verify_arithmetic_layers.py')
layers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(layers)


def psi_bounds(base, digits):
    base = F(base)
    low, high = layers.phi_bounds(base, digits)
    return base*low, base*high


def certified_cf(low, high, requested_terms=12):
    digits = []
    for _ in range(requested_terms):
        quotient = low.__floor__()
        if quotient != high.__floor__():
            return digits, 'ambiguous floor'
        digits.append(quotient)
        low, high = low-quotient, high-quotient
        if low <= 0:
            return digits, 'interval touches a reciprocal pole'
        low, high = 1/high, 1/low
    return digits, 'requested prefix certified'


def convergents(partials):
    old_p, p, old_q, q = 0, 1, 1, 0
    result = []
    for partial in partials:
        old_p, p = p, partial*p+old_p
        old_q, q = q, partial*q+old_q
        result.append(F(p, q))
    return result


def inverse_psi(target, digits, width=F(1, 2**112)):
    left, right = F(3), F(4)
    target = F(target)
    assert digits[0:3] == [3, 0, 1]
    assert psi_bounds(left, digits)[0] > target
    assert psi_bounds(right, digits)[1] < target
    bisections = fallback_steps = 0
    while right-left > width:
        midpoint = (left+right)/2
        low, high = psi_bounds(midpoint, digits)
        if low > target:
            left = midpoint
            bisections += 1
        elif high < target:
            right = midpoint
            bisections += 1
        else:
            # d_3=1 gives -Psi'(t) >= 2/t^3 >= 2/right^3.
            derivative_lower = F(2)/right**3
            radius = max(target-low, high-target)/derivative_lower
            new_left = max(left, midpoint-radius)
            new_right = min(right, midpoint+radius)
            assert new_right-new_left < right-left, 'Increase digit precision'
            left, right = new_left, new_right
            fallback_steps += 1
    residual = (psi_bounds(right, digits)[0]-target,
                psi_bounds(left, digits)[1]-target)
    assert residual[0] <= 0 <= residual[1]
    assert 0 < right-left <= F(1, 2**110)
    return (left, right), residual, {
        'bisections': bisections, 'derivative_fallback_steps': fallback_steps,
        'width_exact': str(right-left), 'requested_maximum_width': '1/2^110'}


def moment_bounds(base, digits):
    base = F(base)
    numerator, denominator = base.numerator, base.denominator
    partial_numerator = 0
    denominator_power = 1
    for index, digit in enumerate(digits, 1):
        denominator_power *= denominator
        partial_numerator = numerator*partial_numerator + index*digit*denominator_power
    degree = len(digits)
    numerator_power = numerator**degree
    partial = F(partial_numerator, numerator_power)
    tail = F(3*denominator_power*denominator
             *((degree+1)*numerator-degree*denominator),
             numerator_power*(numerator-denominator)**2)
    return partial, partial+tail


def absolute_interval(low, high):
    assert low <= high
    if low > 0:
        return low, high
    if high < 0:
        return -high, -low
    return F(0), max(-low, high)


def figure(digits, root_intervals, pi_interval):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return {'generated': False, 'reason': 'matplotlib unavailable'}
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'svg.fonttype': 'none'})
    fig, (curve, chamber) = plt.subplots(1, 2, figsize=(12.4, 5.1), layout='constrained')

    def display_psi(base):
        total = 0.0
        for digit in reversed(digits):
            total = (total+digit)/base
        return base*total

    bases = [math.pi+(10-math.pi)*index/500 for index in range(501)]
    curve.plot(bases, [display_psi(base) for base in bases], color='#356b83', lw=2)
    curve.scatter([math.pi], [math.pi], color='#9b3e70', s=35, zorder=3)
    curve.annotate(r'$\Psi(\pi)=\pi$', (math.pi, math.pi),
                   xytext=(14, -7), textcoords='offset points')
    integers = list(range(4, 11))
    curve.scatter(integers, [display_psi(base) for base in integers],
                  color='#356b83', s=24)
    curve.set(xlim=(3.0, 10.3), ylim=(3.0, 3.155),
              xlabel='Target base t', ylabel=r'$\Psi(t)=t\Phi(t)$',
              title='The frozen pi word moves through CF cylinders')
    curve.text(0.5, 0.79, r'$\Psi(b)=[3;\ b^2-b,\ 1,\ \lfloor(b-2)/3\rfloor,\ldots]$',
               transform=curve.transAxes, fontsize=11, ha='center')
    curve.text(0.5, 0.69, 'for every integer b >= 5',
               transform=curve.transAxes, fontsize=10, ha='center', color='#555555')
    curve.grid(alpha=.2)
    curve.spines[['top', 'right']].set_visible(False)

    midpoint_pi = sum(pi_interval)/2
    # Subtract exact rational intervals before the display conversion;
    # never subtract two machine-precision approximations to pi.
    offsets = {key: float((sum(interval)/2-midpoint_pi)*10**6)
               for key, interval in root_intervals.items()}
    left = offsets['355/113']
    right = offsets['688/219']
    chamber.axvspan(left, right, color='#e4d6ef', alpha=.9)
    chamber.axvline(left, color='#87599e', lw=1.5)
    chamber.axvline(right, color='#87599e', lw=1.5)
    chamber.axvline(0, color='#222222', lw=1.3)
    chamber.set(xlim=(left-28, right+24), ylim=(0, 1), yticks=[],
                xlabel=r'$(t-\pi)\,10^6$',
                title='A chamber retaining the prefix [3; 7, 15, 1]')
    chamber.text(.5, .90, r'$688/219<\Psi(t)<355/113$',
                 ha='center', transform=chamber.transAxes, fontsize=12)
    chamber.annotate('Psi(t) = 355/113', (left, .18), xytext=(35, .07),
                     textcoords='data', arrowprops={'arrowstyle': '->', 'color': '#87599e'},
                     fontsize=10, color='#67417b')
    chamber.annotate('Psi(t) = 688/219', (right, .35), xytext=(right-135, .30),
                     textcoords='data', arrowprops={'arrowstyle': '->', 'color': '#87599e'},
                     fontsize=10, color='#67417b')
    chamber.spines[['top', 'right', 'left']].set_visible(False)
    chamber.grid(axis='x', alpha=.2)
    inset = chamber.inset_axes([.30, .49, .62, .28])
    inset.axvspan(left, .8, color='#e4d6ef', alpha=.9)
    inset.axvline(left, color='#87599e', lw=1.5)
    inset.axvline(0, color='#222222', lw=1.2)
    inset.set(xlim=(left-.7, .8), ylim=(0, 1), yticks=[],
              xticks=[-3, -2, -1, 0], title='Left edge, magnified')
    inset.tick_params(labelsize=8)
    inset.title.set_fontsize(9)
    inset.text(.12, .55, '355/113', transform=inset.transAxes, fontsize=8, color='#67417b')
    inset.text(.81, .55, 'pi', transform=inset.transAxes, fontsize=8)
    inset.spines[['top', 'right', 'left']].set_visible(False)
    fig.suptitle('Continued fractions along a continuously changing radix', fontsize=15)
    fig.savefig(HERE/'cf_migration_geometry.svg', bbox_inches='tight', pad_inches=.18)
    fig.savefig('/tmp/cf_migration_geometry_preview.png', dpi=140,
                bbox_inches='tight', pad_inches=.18)
    plt.close(fig)
    return {'generated': True, 'file': 'cf_migration_geometry.svg',
            'scope': 'Floating-point illustration; boundary offsets computed from certified rational intervals before conversion.'}


def main():
    degree = 192
    pi_low, pi_high, digits = layers.transport.certify_digits(degree)
    records = []
    for base in list(range(4, 101))+[1000]:
        value = psi_bounds(F(base), digits)
        partials, status = certified_cf(*value)
        assert len(partials) == 12
        assert partials[0:2] == [3, base*base-base]
        if base >= 5:
            assert partials[0:4] == [3, base*base-base, 1, (base-2)//3]
        records.append({'base': base, 'Psi_interval': layers.interval_record(*value),
                        'continued_fraction_prefix': partials, 'status': status})
    pi_cf, pi_cf_status = certified_cf(pi_low, pi_high)
    pi_convergents = convergents(pi_cf)
    targets = [F(22, 7), F(333, 106), F(355, 113), F(103993, 33102), F(688, 219)]
    assert all(target in pi_convergents for target in targets[:4])
    assert targets[-1] == F(355+333, 113+106)
    root_records = []
    root_intervals = {}
    for target in targets:
        root, residual, algorithm = inverse_psi(target, digits)
        offset = (root[0]-pi_high, root[1]-pi_low)
        absolute_offset = absolute_interval(*offset)
        error = absolute_interval(pi_low-target, pi_high-target)
        assert error[0] > 0 and absolute_offset[0] > 0
        scale = target.denominator**2
        ratio = (absolute_offset[0]/error[1], absolute_offset[1]/error[0])
        record = {
            'target_p_over_q': str(target), 'is_pi_convergent': target in pi_convergents,
            'root_exact_rational_interval': [str(root[0]), str(root[1])],
            'root': layers.interval_record(*root),
            'root_minus_pi': layers.interval_record(*offset),
            'q_squared_times_signed_offset': layers.interval_record(scale*offset[0], scale*offset[1]),
            'q_squared_times_absolute_offset': layers.interval_record(scale*absolute_offset[0], scale*absolute_offset[1]),
            'kappa_q_squared_abs_pi_minus_p_over_q': layers.interval_record(scale*error[0], scale*error[1]),
            'absolute_offset_over_rational_error': layers.interval_record(*ratio),
            'Psi_root_minus_target': layers.interval_record(*residual),
            'algorithm': algorithm}
        root_records.append(record)
        root_intervals[str(target)] = root
    m_low = moment_bounds(pi_high, digits)[0]
    m_high = moment_bounds(pi_low, digits)[1]
    assert 1 < m_low < m_high
    convergence_constant = (1/(m_high-1), 1/(m_low-1))
    # Exercise the ambiguity-safe derivative fallback with a deliberately
    # chosen rational target inside the enclosure of Psi(7/2).
    midpoint_value = psi_bounds(F(7, 2), digits)
    fallback_root, _, fallback_algorithm = inverse_psi(sum(midpoint_value)/2, digits)
    assert fallback_algorithm['derivative_fallback_steps'] >= 1
    result = {
        'native_digit_certificate': {'count': degree, 'method': 'Machin rational bounds and 600-bit outward intervals'},
        'pi_continued_fraction_prefix': pi_cf, 'pi_prefix_status': pi_cf_status,
        'integer_base_checks': records,
        'universal_prefix_finitely_checked': {
            'formula_for_b_at_least_5': '[3; b^2-b, 1, floor((b-2)/3)]',
            'bases': '5 through 100, and 1000',
            'exception_b4_prefix': records[0]['continued_fraction_prefix']},
        'selected_bases': [4, 5, 6, 10, 20, 100, 1000],
        'rational_boundary_preimages': root_records,
        'm_sum_n_d_n_pi_to_minus_n': layers.interval_record(m_low, m_high),
        'convergence_constant_1_over_m_minus_1': layers.interval_record(*convergence_constant),
        'ambiguity_fallback_check': {'algorithm': fallback_algorithm,
                                     'root': layers.interval_record(*fallback_root)},
        'figure': figure(digits, root_intervals, (pi_low, pi_high)),
        'scope': 'Finite exact certificates; the universal prefix and inverse asymptotic claims require the accompanying mathematical proofs. No P2848 conclusion.'
    }
    (HERE/'cf_migration_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'native_digits': degree, 'integer_bases_checked': len(records),
                      'continued_fraction_terms_per_base': 12, 'boundary_preimages': len(root_records),
                      'files': ['cf_migration_results.json', 'cf_migration_geometry.svg']}, indent=2))


if __name__ == '__main__':
    main()
