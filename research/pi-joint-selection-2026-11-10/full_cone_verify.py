"""Independent exact checks for finite full-cone certificate optimization.

Imports only the independent verifier's mathematical helpers, never either
constructor. The completeness theorem is proved in THEORY.md; this program
checks its finite hypotheses, every complete candidate list and each witness.
"""
from fractions import Fraction as F
from math import gcd, isqrt, lcm
from pathlib import Path
import hashlib
import json
import sys

from verify import sources, root_convergents, pi_interval, quality_interval

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parent


def score(p, q, lower, upper):
    """The exact minimum triangle certificate over all interval sources."""
    squared = F(p * p, q * q)
    if squared < lower:
        chosen = lower
    elif squared > upper:
        chosen = upper
    else:
        chosen = squared
    distance = abs(squared - chosen)
    return q * q * (upper - lower + distance) / 6, chosen


def candidates(lower, upper, qlimit):
    choices = set()
    checked = 0
    for end in [lower, upper]:
        for p, q in root_convergents(end.numerator, end.denominator, qlimit):
            if 3 * q <= p <= 4 * q:
                choices.add((p, q))
                checked += 1
    return choices, checked


def exhaustive_check(lower, upper, qlimit, eligible, expected):
    """Literal integer enumeration for small q, without the CF criterion."""
    den = lcm(lower.denominator, upper.denominator)
    ll = lower.numerator * (den // lower.denominator)
    uu = upper.numerator * (den // upper.denominator)
    width = uu - ll
    incumbent = uu - 9 * den
    best = None
    scanned = 0
    competitive = 0
    for q in range(1, qlimit + 1):
        q2 = q * q
        left, right = ll * q2, uu * q2
        for p in range(3 * q, 4 * q + 1):
            scanned += 1
            value = p * p * den
            displacement = max(left - value, value - right, 0)
            scaled = width * q2 + displacement
            if scaled > incumbent or gcd(p, q) != 1:
                continue
            competitive += 1
            assert (p, q) in eligible
            row = (scaled, q, p)
            if best is None or row < best:
                best = row
    assert best is not None
    assert (F(best[0], 6 * den), best[1], best[2]) == expected
    return {'q_limit': qlimit, 'integer_pairs_scanned': scanned,
            'competitive_primitive_pairs': competitive}


def rational_root_controls():
    controls = [(126025, 12769, 113, [(3, 1), (22, 7), (355, 113)]),
                (126025, 12769, 112, [(3, 1), (22, 7)]),
                (961, 100, 1000, [(3, 1), (31, 10)]),
                (9, 1, 1000, [(3, 1)])]
    for A, B, cap, expected in controls:
        assert list(root_convergents(A, B, cap)) == expected
    return len(controls)


def expand_compact(data):
    """Reconstruct and hash omitted exact records without constructor imports."""
    assert data['schema'] == 2
    nmax = data['parameters']['n_max']
    _, rs = sources(nmax)
    table = data['distinct_fractions']
    assert [row['id'] for row in table] == list(range(len(table)))
    assert len({(f['p'], f['q']) for f in table}) == len(table)
    expanded, exact_records = [], {}
    seen_ids = set()
    transition_indices = []
    previous = None
    for n, compact in enumerate(data['cases'], 1):
        assert compact['n'] == n
        r0, r1 = rs[n], rs[n + 1]
        lo, hi = sorted([r0, r1])
        delta = hi - lo
        C0 = (hi - 9) / 6
        cap = 6 * C0 / delta
        limit = isqrt(cap.numerator // cap.denominator)
        choices, checked = candidates(lo, hi, limit)
        value, q, p, R = min((score(p, q, lo, hi)[0], q, p, score(p, q, lo, hi)[1])
                             for p, q in choices)
        # Solve y/x from barycentric weights, independently of the
        # constructor's determinant/adjugate formula.
        if R == r0:
            x, y = 1, 0
        elif R == r1:
            x, y = 0, 1
        else:
            ratio = r0.denominator * (R - r0) / (r1.denominator * (r1 - R))
            assert ratio > 0
            x, y = ratio.denominator, ratio.numerator
        rawA = x * r0.numerator + y * r1.numerator
        rawB = x * r0.denominator + y * r1.denominator
        g = gcd(rawA, rawB)
        A, B = rawA // g, rawB // g
        h = B * p * p - A * q * q
        position = 'below' if F(p * p, q * q) < lo else 'above' if F(p * p, q * q) > hi else 'inside'
        w = {k: str(v) for k, v in zip(['x', 'y', 'g', 'A', 'B', 'h'], [x, y, g, A, B, h])}
        record = {'n': n, 'source0': str(r0), 'source1': str(r1),
                  'L': str(lo), 'U': str(hi), 'delta': str(delta), 'C0': str(C0),
                  'q_limit': str(limit), 'endpoint_convergents_checked': checked,
                  'unique_candidates': len(choices), 'p': str(p), 'q': str(q),
                  'certificate': str(value), 'R': str(R), 'witness': w,
                  'position': position}
        digest = hashlib.sha256(json.dumps(record, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        assert digest == compact['exact_record_sha256']
        assert compact['q_limit_bits'] == limit.bit_length()
        assert compact['endpoint_convergents_checked'] == checked
        assert compact['unique_candidates'] == len(choices)
        assert compact['position'] == position
        fid = compact['fraction_id']
        assert (str(p), str(q)) == (table[fid]['p'], table[fid]['q'])
        seen_ids.add(fid)
        if fid != previous:
            transition_indices.append(n)
        previous = fid
        exact_records[n] = record
        expanded.append({**record,
                         'certificate_decimal_enclosure': compact['certificate_decimal_enclosure'],
                         'actual_quality_interval': table[fid]['actual_quality_interval'],
                         'actual_quality_decimal_enclosure': table[fid]['actual_quality_decimal_enclosure']})
    assert seen_ids == set(range(len(table)))
    retained = sorted(set(transition_indices + [nmax]))
    assert [r['n'] for r in data['witness_records']] == retained
    for record in data['witness_records']:
        assert record == exact_records[record['n']]
    return {**data, 'schema': 1, 'cases': expanded}


def run():
    path = ROOT / 'full_cone_results.json'
    data = json.loads(path.read_text())
    original_schema = data['schema']
    if original_schema == 2:
        data = expand_compact(data)
    nmax = data['parameters']['n_max']
    assert data['schema'] == 1 and data['parameters']['n_min'] == 1
    assert data['parameters']['fixed_norm_budgets'] is False
    assert data['parameters']['fraction_domain'] == 'coprime p,q; q>0; 3<=p/q<=4'
    assert data['source_commit'] == '88c331c410290c5eb915776ecd2415bffff20dc0'
    assert len(data['cases']) == nmax
    u, rs = sources(nmax)
    widths = [F(30, (n + 1) ** 2 * u[n] * u[n + 1]) for n in range(1, nmax + 1)]
    bits = max(d.denominator.bit_length() - d.numerator.bit_length() for d in widths) + 512
    plo, phi, terms = pi_interval(bits)
    total_checked = total_unique = 0
    exhaustive = []
    for n, c in enumerate(data['cases'], 1):
        assert c['n'] == n
        r0, r1 = rs[n], rs[n + 1]
        lo, hi = sorted([r0, r1])
        delta = widths[n - 1]
        assert 9 < lo < hi <= 10
        assert delta == hi - lo
        assert lo < plo * plo < phi * phi < hi
        assert [r0, r1, lo, hi, delta] == [F(c[k]) for k in ['source0', 'source1', 'L', 'U', 'delta']]
        C0 = (hi - 9) / 6
        assert 0 < C0 <= F(1, 6) and C0 == F(c['C0'])
        cap = 6 * C0 / delta
        limit = isqrt(cap.numerator // cap.denominator)
        assert limit == int(c['q_limit'])
        assert limit * limit <= cap < (limit + 1) ** 2
        choices, checked = candidates(lo, hi, limit)
        assert (3, 1) in choices
        assert len(choices) == c['unique_candidates']
        assert checked == c['endpoint_convergents_checked']
        value, q, p, R = min((score(p, q, lo, hi)[0], q, p, score(p, q, lo, hi)[1])
                             for p, q in choices)
        assert value <= C0
        assert (value, p, q, R) == (F(c['certificate']), int(c['p']), int(c['q']), F(c['R']))
        dl, du = map(F, c['certificate_decimal_enclosure'])
        assert dl <= value <= du
        assert gcd(p, q) == 1 and 3 * q <= p <= 4 * q
        w = {key: int(val) for key, val in c['witness'].items()}
        x, y, g, A, B, h = (w[k] for k in ['x', 'y', 'g', 'A', 'B', 'h'])
        assert x >= 0 and y >= 0 and x + y > 0 and gcd(x, y) == 1
        rawA = x * r0.numerator + y * r1.numerator
        rawB = x * r0.denominator + y * r1.denominator
        assert g == gcd(rawA, rawB) and g > 0
        assert (rawA, rawB) == (g * A, g * B)
        assert gcd(A, B) == 1 and F(A, B) == R
        assert h == B * p * p - A * q * q
        assert value == (F(abs(h), B) + q * q * delta) / 6
        sq = F(p * p, q * q)
        position = 'below' if sq < lo else 'above' if sq > hi else 'inside'
        assert c['position'] == position
        assert (position == 'inside') == (h == 0)
        ql, qh = map(F, c['actual_quality_interval'])
        dl, du = map(F, c['actual_quality_decimal_enclosure'])
        tl, tu = quality_interval(p, q, plo, phi)
        assert 0 < dl <= ql <= tl < tu <= qh <= du and qh < value
        if n <= 4:
            e = exhaustive_check(lo, hi, limit, choices, (value, q, p))
            exhaustive.append({'n': n, **e})
        total_checked += checked
        total_unique += len(choices)
        if n % 64 == 0:
            print(json.dumps({'verified_full_cone_n': n, 'endpoint_convergents': total_checked}), flush=True)
    cases = data['cases']
    best = min(cases, key=lambda c: (F(c['certificate']), int(c['q']), int(c['p'])))
    actual = min(cases, key=lambda c: F(c['actual_quality_interval'][1]))
    for c in cases:
        if (c['p'], c['q']) != (actual['p'], actual['q']):
            assert F(actual['actual_quality_interval'][1]) < F(c['actual_quality_interval'][0])
    expected_summary = {
        'source_intervals': nmax, 'endpoint_convergents_checked': total_checked,
        'unique_candidates_per_interval_sum': total_unique,
        'distinct_selected_fractions': len({(c['p'], c['q']) for c in cases}),
        'inside_selected': sum(c['position'] == 'inside' for c in cases),
        'best_certificate_n': best['n'], 'best_certificate': best['certificate'],
        'best_certificate_decimal_enclosure': best['certificate_decimal_enclosure'],
        'best_certificate_fraction': best['p'] + '/' + best['q'],
        'best_actual_fraction': actual['p'] + '/' + actual['q'],
        'best_actual_n': actual['n'],
        'best_actual_quality_decimal_enclosure': actual['actual_quality_decimal_enclosure'],
        'scoring_bits': max(2 * int(c['q']).bit_length() + 256 for c in cases)}
    assert expected_summary == data['summary']
    # Stored score bounds are themselves reproduced at their original scale;
    # the strictly narrower, independently generated bounds above check them.
    score_lo, score_hi, _ = pi_interval(expected_summary['scoring_bits'])
    for c in cases:
        assert quality_interval(int(c['p']), int(c['q']), score_lo, score_hi) == tuple(map(F, c['actual_quality_interval']))
    report = {
        'checks': 'PASS', 'file': path.name,
        'input_schema': original_schema,
        'compact_exact_record_hashes_recomputed': nmax if original_schema == 2 else 0,
        'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'source_intervals_checked': nmax,
        'independent_endpoint_convergents_checked': total_checked,
        'unique_candidate_sum': total_unique,
        'optimal_source_and_primitive_weight_witnesses': nmax,
        'summary_recomputed_exactly': True,
        'independent_pi_grid_bits': bits, 'independent_taylor_terms': terms,
        'exhaustive_integer_cross_checks': exhaustive,
        'rational_root_terminal_controls': rational_root_controls(),
        'best_actual_quality_decimal_enclosure': actual['actual_quality_decimal_enclosure'],
        'method': 'Independent binomial/telescoping source sums, dyadic endpoint CFs, exact envelope minimum, direct integer witness checks, outward Taylor intervals.',
        'scope': 'Finite n only. Global minimum refers only to the stated uniform-tail triangle certificate over the full positive source cone, not actual pi quality or a new global record.'}
    (ROOT / 'full_cone_verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)
    return report


if __name__ == '__main__':
    run()
