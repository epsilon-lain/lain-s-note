"""Exact finite checks and a display figure for frozen-word base transport.

Run: python3 verify_transport.py
All certificate decisions use integers and Fraction. The optional SVG uses
floating point only for display; it is not a proof of the infinite theorems.
"""

import importlib.util
import itertools
import json
import math
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
HELPER = HERE.parent / 'pi-native-space-2026-10-31' / 'verify_native.py'
spec = importlib.util.spec_from_file_location('native_certificates', HELPER)
native = importlib.util.module_from_spec(spec)
spec.loader.exec_module(native)


def outward(value, rounding):
    with localcontext() as ctx:
        ctx.prec = 30
        ctx.rounding = rounding
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def interval_record(low, high):
    assert low <= high
    return {'lower': outward(low, ROUND_FLOOR),
            'upper': outward(high, ROUND_CEILING)}


def certify_digits(count):
    low0, high0 = native.pi_bounds(160)
    grid = 2**600
    low = F((low0 * grid).__floor__(), grid)
    high = F((high0 * grid).__ceil__(), grid)
    assert low <= low0 < high0 <= high
    remainder_low = remainder_high = F(1)
    digits = []
    for _ in range(count):
        raw_low, raw_high = low * remainder_low, high * remainder_high
        digit = raw_low.__floor__()
        assert digit == raw_high.__floor__(), 'Increase precision'
        assert 0 <= digit <= 3
        digits.append(digit)
        remainder_low = F(((raw_low-digit)*grid).__floor__(), grid)
        remainder_high = F(((raw_high-digit)*grid).__ceil__(), grid)
        assert 0 <= remainder_low <= remainder_high < 1
    return low, high, digits


def evaluate_word(word, base):
    result = F(0)
    for digit in reversed(word):
        result = (digit + result) / base
    return result


def bounded_unit_and_moment(digits, base):
    count = len(digits)
    ratio = F(1, base)
    unit_low = sum((digit * ratio**n
                    for n, digit in enumerate(digits, 1)), F(0))
    moment_low = sum((n * digit * ratio**n
                      for n, digit in enumerate(digits, 1)), F(0))
    unit_tail = 3 * ratio**(count+1) / (1-ratio)
    moment_tail = (3 * ratio**(count+1)
                   * ((count+1)-count*ratio) / (1-ratio)**2)
    return (unit_low, unit_low+unit_tail), (moment_low, moment_low+moment_tail)


def finite_words(digits, pi_low, pi_high, counts, max_length=7):
    legal_by_length = {0: [()]}
    checks = []
    for length in range(1, max_length+1):
        legal = []
        for word in itertools.product(range(4), repeat=length):
            # Every finite suffix is less than the infinite boundary word.
            # Equality through the end is allowed since the certified boundary
            # has a nonzero digit later and is not eventually zero.
            lex_legal = all(word[start:] <= tuple(digits[:length-start])
                            for start in range(length))
            # An independent interval check of every suffix's numeric value.
            numeric_legal = True
            for start in range(length):
                suffix = word[start:]
                lower = evaluate_word(suffix, pi_high)
                upper = evaluate_word(suffix, pi_low)
                assert upper < 1 or lower > 1, 'Unresolved endpoint'
                if lower > 1:
                    numeric_legal = False
                    break
            assert lex_legal == numeric_legal
            if lex_legal:
                legal.append(word)
        assert len(legal) == sum(counts[:length+1])
        terminal = sum(word[-1] > 0 for word in legal)
        assert terminal == counts[length]
        legal_by_length[length] = legal
        checks.append({'length': length, 'legal_words': len(legal),
                       'positive_last_digit_words': terminal})
    return legal_by_length, checks


def gap_checks(base, unit, legal_by_length, counts, depth=5):
    low_unit, high_unit = unit
    gaps = []
    per_depth = []
    for length in range(1, depth+1):
        count = 0
        for word in legal_by_length[length]:
            if word[-1] == 0:
                continue
            right = evaluate_word(word, F(base))
            scale = F(1, base**length)
            left_low = right-scale+scale*low_unit
            left_high = right-scale+scale*high_unit
            assert left_high < right
            assert right-left_high == scale*(1-high_unit)
            assert right-left_low == scale*(1-low_unit)
            gaps.append((left_low, left_high, right))
            count += 1
        assert count == counts[length]
        per_depth.append(count)
    gaps.sort(key=lambda gap: gap[2])
    for previous, following in zip(gaps, gaps[1:]):
        assert previous[2] < following[0]
    assert gaps[0][0] > 0 and gaps[-1][2] < low_unit
    return {'depth': depth, 'gap_counts_by_depth': per_depth,
            'strictly_disjoint_gaps_checked': len(gaps)}


def render_figure(legal_by_length, values):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
    except ImportError:
        return {'generated': False, 'reason': 'matplotlib not installed'}

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'svg.fonttype': 'none'})
    figure = plt.figure(figsize=(12, 6.3), layout='constrained')
    grid = figure.add_gridspec(1, 2, width_ratios=(1.65, 1))
    covers = figure.add_subplot(grid[0, 0])
    dimension = figure.add_subplot(grid[0, 1])
    rows = [(math.pi, 3), (4, 2), (10, 1), (16, 0)]
    covers.plot([0, 1], [3, 3], color='#385f91', lw=7, solid_capstyle='butt')
    for base, row in rows[1:]:
        unit = float(sum(values[base]['unit']) / 2)
        segments = []
        for word in legal_by_length[2]:
            start = float(evaluate_word(word, F(base))) / unit
            end = min(1.0, start + base**-2)
            if start < end:
                segments.append([(start, row), (end, row)])
        covers.add_collection(LineCollection(segments, colors='#268273', lw=7))
    covers.set(xlim=(-0.015, 1.015), ylim=(-0.55, 3.55),
               yticks=[3, 2, 1, 0],
               yticklabels=['t = pi: interval', 't = 4', 't = 10', 't = 16'],
               xlabel='Position divided by U(t) (each row has hull [0, 1])',
               title='One word space, several target bases')
    covers.text(0.02, 3.32, 'Exact interval at t = pi', color='#385f91', fontsize=10)
    covers.text(0.02, -0.40, 'Other rows: depth-2 outer covers, not a finite point cloud',
                color='#46515d', fontsize=9)
    covers.spines[['top', 'right', 'left']].set_visible(False)
    covers.tick_params(axis='y', length=0)

    bases = [math.pi + (16-math.pi)*j/500 for j in range(501)]
    dimensions = [math.log(math.pi)/math.log(base) for base in bases]
    dimension.plot(bases, dimensions, color='#803e75', lw=2)
    for base in [4, 10, 16]:
        value = math.log(math.pi)/math.log(base)
        dimension.scatter([base], [value], color='#803e75', s=22)
        dimension.annotate(f'{base}: {value:.4f}', (base, value),
                           xytext=(8, 8) if base == 4 else (-4, 10),
                           textcoords='offset points', ha='left' if base == 4 else 'right',
                           fontsize=10)
    dimension.set(xlim=(math.pi, 16.5), ylim=(0.35, 1.05), xlabel='Target base t',
                  ylabel='Hausdorff dimension',
                  title='Dimension shrinks;\nsymbolic growth stays fixed')
    dimension.text(0.48, 0.88, r'$\dim_H K_t=\log\pi/\log t$',
                   transform=dimension.transAxes, ha='center', fontsize=13)
    dimension.grid(alpha=0.2)
    dimension.spines[['top', 'right']].set_visible(False)
    figure.suptitle('Frozen digits preserve the language and change its geometric scale',
                    fontsize=15)
    figure.savefig(HERE / 'transport_geometry.svg', bbox_inches='tight', pad_inches=0.2)
    # Temporary local raster for visual inspection; only the SVG is published.
    figure.savefig('/tmp/base_transport_geometry_preview.png', dpi=140,
                   bbox_inches='tight', pad_inches=0.2)
    plt.close(figure)
    return {'generated': True, 'file': 'transport_geometry.svg',
            'scope': 'Floating-point display of proved formulas and finite outer covers.'}


def main():
    degree = 192
    pi_low, pi_high, digits = certify_digits(degree)
    counts = [1]
    for length in range(1, degree+1):
        counts.append(sum(digits[k-1]*counts[length-k]
                          for k in range(1, length+1)))
    assert counts[:8] == [1, 3, 9, 28, 88, 276, 867, 2724]
    words, word_checks = finite_words(digits, pi_low, pi_high, counts)
    numerical = []
    values = {}
    for base in (4, 10, 16):
        unit, moment = bounded_unit_and_moment(digits, base)
        delta = (1-unit[1], 1-unit[0])
        assert 0 < unit[0] < unit[1] < 1
        mean = (moment[0]/(unit[1]*delta[1]),
                moment[1]/(unit[0]*delta[0]))
        weighted_partial = sum((F(counts[n], base**n)
                                for n in range(1, degree+1)), F(0))
        # c_n*pi^-n <= 1 by induction from the probability renewal law;
        # pi < 22/7 therefore gives this simple geometric tail bound.
        count_ratio = F(22, 7*base)
        count_tail = count_ratio**(degree+1)/(1-count_ratio)
        covered_gap_length = (delta[0]*weighted_partial,
                              delta[1]*weighted_partial)
        missing = (unit[0]-covered_gap_length[1],
                   unit[1]-covered_gap_length[0])
        assert 0 < missing[0] < missing[1] < delta[1]*count_tail
        numeric_record = {
            'base': base, 'U': interval_record(*unit),
            'delta_1_minus_U': interval_record(*delta),
            'M1': interval_record(*moment),
            'length_weighted_mean_gap_depth': interval_record(*mean),
            'finite_gaps': gap_checks(base, unit, words, counts),
            'gap_mass_partial_depth': degree,
            'uncovered_gap_mass': interval_record(*missing),
            'uncovered_gap_mass_upper_bound': outward(delta[1]*count_tail, ROUND_CEILING),
            'dimension_exact': f'log(pi)/log({base})',
            'dimension_display_only': math.log(math.pi)/math.log(base)
        }
        numerical.append(numeric_record)
        values[base] = {'unit': unit}
    result = {
        'native_digit_certificate': {
            'count': degree, 'pi_bounds': 'Machin identity, 160 alternating terms',
            'working_intervals': '600-bit dyadic outward enclosures at each step',
            'prefix': ''.join(map(str, digits[:64]))},
        'return_counts_first_16': counts[:16],
        'finite_word_checks': word_checks,
        'target_base_certificates': numerical,
        'figure': render_figure(words, values),
        'scope': 'Finite exact interval and integer checks. Infinite arithmetic classifications, Cantor geometry, dimensions, and transport asymptotics require the accompanying proofs. No P2848 conclusion.'
    }
    (HERE / 'transport_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'digits_certified': degree,
                      'word_lengths_checked': len(word_checks),
                      'gap_checks_per_base': numerical[0]['finite_gaps']['strictly_disjoint_gaps_checked'],
                      'files': ['transport_results.json', 'transport_geometry.svg']}, indent=2))


if __name__ == '__main__':
    main()
