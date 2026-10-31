"""Exact smallest-denominator certificates for native survival windows.

Run: python3 verify_survival_windows.py
Uses integer polynomial comparisons and Farey-neighbor certificates. The figure
is only a floating-point rendering of rigorously computed intervals.
"""

import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('cf_migration', HERE/'verify_continued_fraction_migration.py')
cf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cf)


def carry_numerator(p, q, digits):
    result, q_power = 1, 1
    for digit in digits:
        q_power *= q
        result = p*result-digit*q_power
    return result, q_power


def carry_class(p, q, digits):
    numerator, denominator = carry_numerator(p, q, digits)
    return -1 if numerator <= -2*denominator else 1 if numerator >= 2*denominator else 0


def simplest(digits):
    left, right = (3, 1), (4, 1)
    cache = {}
    runs = 0

    def classify(pair):
        if pair not in cache:
            cache[pair] = carry_class(*pair, digits)
        return cache[pair]

    assert classify(left) == -1 and classify(right) == 1
    while True:
        assert left[1]*right[0]-left[0]*right[1] == 1
        middle = (left[0]+right[0], left[1]+right[1])
        direction = classify(middle)
        if direction == 0:
            assert classify(left) == -1 and classify(right) == 1
            return F(*middle), F(*left), F(*right), len(cache), runs
        origin, step = (left, right) if direction == -1 else (right, left)

        def candidate(k):
            return origin[0]+k*step[0], origin[1]+k*step[1]

        low, high = 1, 2
        while classify(candidate(high)) == direction:
            low, high = high, 2*high
        while high-low > 1:
            midpoint = (low+high)//2
            if classify(candidate(midpoint)) == direction:
                low = midpoint
            else:
                high = midpoint
        if direction == -1:
            left = candidate(low)
        else:
            right = candidate(low)
        runs += 1


def root_bracket(digits, level):
    left, right = F(3), F(4)
    width = F(1, 2**(2*len(digits)+50))

    def compare(point):
        numerator, denominator = carry_numerator(point.numerator, point.denominator, digits)
        return numerator-level*denominator

    assert compare(left) < 0 < compare(right)
    while right-left > width:
        middle = (left+right)/2
        value = compare(middle)
        if value < 0:
            left = middle
        elif value > 0:
            right = middle
        else:
            return middle, middle
    assert compare(left) <= 0 <= compare(right)
    return left, right


def approximation_type(value, partials):
    convergents = cf.convergents(partials)
    if value in convergents:
        return {'type': 'convergent', 'index': convergents.index(value)}
    pp, p, qp, q = 0, 1, 1, 0
    for index, a in enumerate(partials):
        if index >= 1:
            numerator_diff = value.numerator-pp
            denominator_diff = value.denominator-qp
            if q and denominator_diff % q == 0:
                k = denominator_diff//q
                if 1 <= k <= a and numerator_diff == k*p:
                    return {'type': 'intermediate convergent', 'index': index,
                            'intermediate_coefficient': k, 'full_coefficient': a}
        pp, p = p, a*p+pp
        qp, q = q, a*q+qp
    return {'type': 'not identified in certified CF prefix'}


def absolute_error(value, pi_interval):
    return cf.absolute_interval(value-pi_interval[1], value-pi_interval[0])


def figure(records, plateaus):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return {'generated': False}
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'none', 'font.size': 11})
    fig, axes = plt.subplots(2, 1, figsize=(10.8, 7.5), sharex=True, layout='constrained')
    depths = [r['N'] for r in records]
    denominators = [r['minimum_denominator'] for r in records]
    metric = [float(r['Q_squared_pi_to_minus_N']['lower']) for r in records]
    axes[0].step(depths, denominators, where='post', color='#356b83', lw=1.6)
    axes[0].set(yscale='log', ylabel='Smallest denominator Q_N',
                title='Exact denominator staircase in shrinking survival windows')
    axes[1].plot(depths, metric, color='#87599e', lw=1.4)
    axes[1].set(yscale='log', xlabel='Number N of native digits read',
                ylabel=r'$Q_N^2\,\pi^{-N}$', title='Low points mean a small denominator survives a narrow window')
    for plateau in plateaus:
        if plateau['rational'] == '355/113':
            for axis in axes:
                axis.axvspan(plateau['first_N']-.5, plateau['last_N']+.5,
                             color='#e7bfd2', alpha=.45)
            axes[1].annotate('355/113 plateau',
                             (plateau['last_N'], metric[plateau['last_N']-4]),
                             xytext=(18, 0.03), textcoords='data',
                             arrowprops={'arrowstyle': '->', 'color': '#9b3e70'}, color='#9b3e70')
    for axis in axes:
        axis.grid(alpha=.2)
        axis.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Finite integer certificates; no conclusion about infinite liminf', fontsize=13)
    fig.savefig(HERE/'window_denominator_geometry.svg', bbox_inches='tight', pad_inches=.18)
    fig.savefig('/tmp/window_denominator_geometry_preview.png', dpi=140,
                bbox_inches='tight', pad_inches=.18)
    plt.close(fig)
    return {'generated': True, 'file': 'window_denominator_geometry.svg',
            'scope': 'Floating-point rendering of certified exact comparisons.'}


def main():
    pi_low, pi_high, digits = cf.layers.transport.certify_digits(192)
    partials, cf_status = cf.certified_cf(pi_low, pi_high, 100)
    assert len(partials) == 100
    records = []
    plateaus = []
    previous_q = 0
    total_queries = 0
    selected = {4, 7, 8, 14, 15, 30, 60, 100, 160}
    endpoints = []
    for n in range(4, 161):
        word = digits[:n]
        value, left, right, queries, runs = simplest(word)
        p, q = value.numerator, value.denominator
        assert q >= previous_q
        previous_q = q
        total_queries += queries
        metric = (q*q/pi_high**n, q*q/pi_low**n)
        kappa = tuple(q*q*x for x in absolute_error(value, (pi_low, pi_high)))
        classification = approximation_type(value, partials)
        record = {'N': n, 'rational': str(value), 'minimum_denominator': q,
                  'farey_left_parent': str(left), 'farey_right_parent': str(right),
                  'parent_determinant': left.denominator*right.numerator-left.numerator*right.denominator,
                  'parent_classes': [-1, 1], 'candidate_class': 0,
                  'Q_squared_pi_to_minus_N': cf.layers.interval_record(*metric),
                  'kappa_q_squared_abs_pi_minus_rational': cf.layers.interval_record(*kappa),
                  'relation_to_pi_cf': classification, 'integer_polynomial_queries': queries,
                  'accelerated_runs': runs}
        assert record['parent_determinant'] == 1
        records.append(record)
        if plateaus and plateaus[-1]['rational'] == str(value):
            plateaus[-1]['last_N'] = n
            plateaus[-1]['last_metric'] = record['Q_squared_pi_to_minus_N']
        else:
            plateaus.append({'rational': str(value), 'denominator': q, 'first_N': n, 'last_N': n,
                             'relation_to_pi_cf': classification,
                             'first_metric': record['Q_squared_pi_to_minus_N'],
                             'last_metric': record['Q_squared_pi_to_minus_N']})
        if n in selected:
            low = root_bracket(word, -2)
            high = root_bracket(word, 2)
            assert low[1] < pi_low < pi_high < high[0]
            width_bounds = (high[0]-low[1], high[1]-low[0])
            normalized_width = (width_bounds[0]*pi_low**n, width_bounds[1]*pi_high**n)
            endpoints.append({'N': n,
                              'left_root_exact_interval': [str(x) for x in low],
                              'right_root_exact_interval': [str(x) for x in high],
                              'left_root_minus_pi': cf.layers.interval_record(low[0]-pi_high, low[1]-pi_low),
                              'right_root_minus_pi': cf.layers.interval_record(high[0]-pi_high, high[1]-pi_low),
                              'width_pi_to_N': cf.layers.interval_record(*normalized_width),
                              'root_interval_width_bound': str(F(1, 2**(2*n+50)))})
    assert all(r['relation_to_pi_cf']['type'] != 'not identified in certified CF prefix' for r in records)
    # A small denominator exhaustive check, independent from Stern-Brocot.
    exhaustive = []
    for row in records:
        if row['minimum_denominator'] > 150:
            continue
        n, minimum = row['N'], row['minimum_denominator']
        for q in range(1, minimum):
            for p in range(3*q, 4*q+1):
                assert carry_class(p, q, digits[:n]) != 0
        exhaustive.append(n)
    result = {'native_digit_certificate': {'count': 192,
                'method': 'Machin rational bounds and 600-bit dyadic outward intervals'},
              'pi_cf_certificate': {'partial_quotients': partials, 'status': cf_status},
              'definitions': {'window': 'I_N={t in [3,4]: |K_N(t)|<2}, N>=4',
                              'K_N': 't^N-sum(d_j*t^(N-j), j=1..N)',
                              'Q_N': 'smallest reduced denominator in I_N'},
              'range': {'first_N': 4, 'last_N': 160, 'depths_certified': len(records),
                        'integer_polynomial_queries': total_queries,
                        'exhaustive_denominator_checks_at_depths': exhaustive},
              'records': records, 'plateaus': plateaus,
              'selected_endpoint_certificates': endpoints,
              'figure': figure(records, plateaus),
              'scope': 'Finite exact rational polynomial decisions and Farey minimality certificates. No infinite-liminf or P2848 conclusion.'}
    (HERE/'survival_window_results.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'range': result['range'], 'plateaus': plateaus,
                      'selected_widths': [{'N': x['N'], 'width_pi_to_N': x['width_pi_to_N']} for x in endpoints]}, indent=2))


if __name__ == '__main__':
    main()
