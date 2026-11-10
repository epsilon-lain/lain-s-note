"""Independent verification of the finite Euler-filter construction.

Reconstructs each rational by polynomial division and quotient integration,
never by the generator's weighted Taylor sections. Beta moments use monomial
integration rather than the generator's product recurrence. No constructor is
imported. A separate dyadic Machin implementation checks quality intervals.
"""
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from functools import cache
from math import comb, gcd, lcm
from pathlib import Path
import hashlib
import json
import sys

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / 'pi-joint-selection-2026-11-10'))
from verify import pi_interval, quality_interval


@cache
def reconstruct_angle(s, t, n, k):
    """Divide z^n(1-z)^k by 1+c^2*z, then integrate z=x^2."""
    c = F(s, t)
    csq = c * c
    polynomial = [F(0)] * (n + k + 1)
    for j in range(k + 1):
        polynomial[n + j] = F((-1) ** j * comb(k, j))
    work = polynomial[:]
    quotient = [F(0)] * (n + k)
    for degree in range(n + k, 0, -1):
        value = work[degree] / csq
        quotient[degree - 1] = value
        work[degree] -= csq * value
        work[degree - 1] -= value
    assert all(a == 0 for a in work[1:])
    remainder = work[0]
    prefactor = (-1) ** n * c ** (2 * n + 1) * (csq / (1 + csq)) ** k
    assert prefactor * remainder / c == 1
    # atan(c) - prefactor*(integral Q(x^2)dx + remainder*atan(c)/c).
    rational = -prefactor * sum((a / (2 * i + 1) for i, a in enumerate(quotient)), F(0))
    beta = sum((F((-1) ** j * comb(k, j), 2 * n + 2 * j + 1)
                for j in range(k + 1)), F(0))
    assert beta > 0
    magnitude = abs(prefactor) * beta
    low, high = magnitude / (1 + csq), magnitude
    interval = (low, high) if n % 2 == 0 else (-high, -low)
    odd_lcm = lcm(*(range(1, 2 * (n + k), 2)))
    raw = t ** (2 * n - 1) * (s * s + t * t) ** k * odd_lcm
    assert raw % rational.denominator == 0
    return rational, interval, raw


def display(lo, hi):
    with localcontext() as context:
        context.prec = 18
        context.rounding = ROUND_FLOOR
        dl = str(Decimal(lo.numerator) / Decimal(lo.denominator))
        context.rounding = ROUND_CEILING
        du = str(Decimal(hi.numerator) / Decimal(hi.denominator))
    return [dl, du]


def digest(row):
    return hashlib.sha256(json.dumps(row, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def run():
    path = ROOT / 'fresh_probe_results.json'
    data = json.loads(path.read_text())
    assert data['schema'] == 2
    params = data['parameters']
    assert params['nmax'] == 40 and params['n_min'] == 1
    families = {'gregory': [[4, 1, 1]],
                'two_angle': [[4, 1, 2], [4, 1, 3]],
                'machin': [[16, 1, 5], [-4, 1, 239]]}
    assert params['families'] == families
    bits = params['pi_score_bits']
    lo, hi, _ = pi_interval(bits)
    tight_lo, tight_hi, terms = pi_interval(bits + 128)
    assert lo < tight_lo < tight_hi < hi
    rows, full = [], []
    for family, angles in families.items():
        for total in range(1, params['nmax'] + 1):
            for k in range(total):
                n = total - k
                r, el, eh, raw = F(0), F(0), F(0), 1
                for coefficient, s, t in angles:
                    component, bounds, cost = reconstruct_angle(s, t, n, k)
                    r += coefficient * component
                    low, high = (coefficient * bounds[0], coefficient * bounds[1])
                    if coefficient < 0:
                        low, high = high, low
                    el += low
                    eh += high
                    raw = lcm(raw, cost)
                p, q = r.numerator, r.denominator
                assert gcd(p, q) == 1 and raw % q == 0
                ql, qh = quality_interval(p, q, lo, hi)
                tl, th = quality_interval(p, q, tight_lo, tight_hi)
                assert 0 < ql < tl < th < qh
                assert el < tight_lo - r < tight_hi - r < eh
                exact = {'family': family, 'total_order': total, 'n': n, 'k': k,
                         'p': str(p), 'q': str(q), 'raw_denominator': str(raw),
                         'raw_gcd': str(raw // q),
                         'error_interval': list(map(str, (el, eh))),
                         'quality_interval': list(map(str, (ql, qh)))}
                row = {'family': family, 'total_order': total, 'n': n, 'k': k,
                       'q_bits': q.bit_length(), 'raw_denominator_bits': raw.bit_length(),
                       'raw_gcd_bits': (raw // q).bit_length(),
                       'actual_quality_decimal': display(ql, qh), 'record_sha256': digest(exact)}
                full.append(exact)
                rows.append(row)
    assert data['family_names'] == list(families)
    assert data['candidate_columns'] == ['family_id', 'total_order', 'n', 'k', 'record_sha256']
    compact = [[list(families).index(r['family']), r['total_order'], r['n'], r['k'],
                r['record_sha256']] for r in rows]
    assert compact == data['candidates']
    selected, summaries = {}, {}
    all_improvements = 0
    for family in families:
        own = [i for i, r in enumerate(rows) if r['family'] == family]
        best = min(own, key=lambda i: F(full[i]['quality_interval'][1]))
        improvements, wins = [], []
        for total in range(1, params['nmax'] + 1):
            ids = [i for i in own if rows[i]['total_order'] == total]
            base = next(i for i in ids if rows[i]['k'] == 0)
            win = min(ids, key=lambda i: F(full[i]['quality_interval'][1]))
            for i in ids:
                if (full[i]['p'], full[i]['q']) != (full[win]['p'], full[win]['q']):
                    assert F(full[win]['quality_interval'][1]) < F(full[i]['quality_interval'][0])
            if F(full[win]['quality_interval'][1]) < F(full[base]['quality_interval'][0]):
                improvements.append(total)
            wins.append({'total_order': total, 'winner_k': rows[win]['k'],
                         'winner_record_sha256': rows[win]['record_sha256'],
                         'baseline_record_sha256': rows[base]['record_sha256']})
            if total in (1, 3, 40):
                selected[str((family, total, 'winner'))] = full[win]
                selected[str((family, total, 'baseline'))] = full[base]
        for i in own:
            if (full[i]['p'], full[i]['q']) != (full[best]['p'], full[best]['q']):
                assert F(full[best]['quality_interval'][1]) < F(full[i]['quality_interval'][0])
        selected[str((family, 'best'))] = full[best]
        summaries[family] = {
            'candidate_count': len(own),
            'best_parameters': {key: rows[best][key] for key in ('total_order', 'n', 'k')},
            'best_fraction': full[best]['p'] + '/' + full[best]['q'],
            'best_actual_quality_decimal': rows[best]['actual_quality_decimal'],
            'quality_improved_total_orders': improvements, 'winners': wins}
        win = next(i for i in own if rows[i]['total_order'] == 40 and
                   rows[i]['record_sha256'] == wins[-1]['winner_record_sha256'])
        base = next(i for i in own if rows[i]['total_order'] == 40 and rows[i]['k'] == 0)
        wl, wu = map(F, full[win]['quality_interval'])
        bl, bu = map(F, full[base]['quality_interval'])
        summaries[family]['matched_order_40'] = {
            'winner_n': rows[win]['n'], 'winner_k': rows[win]['k'],
            'winner_q_bits': rows[win]['q_bits'], 'baseline_q_bits': rows[base]['q_bits'],
            'quality_ratio_decimal': display(wl / bu, wu / bl),
            'winner_quality_decimal': rows[win]['actual_quality_decimal']}
        assert improvements == list(range(2, 41))
        all_improvements += len(improvements)
    assert summaries == data['summary']
    assert selected == data['selected_witnesses']
    assert len(rows) == 2460 and all_improvements == 117
    report = {'checks': 'PASS', 'file': path.name, 'input_schema': data['schema'],
              'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
              'candidate_rationals_reconstructed': len(rows),
              'polynomial_division_angle_reconstructions': reconstruct_angle.cache_info().currsize,
              'exact_record_hashes_checked': len(rows),
              'signed_beta_error_intervals_checked': len(rows),
              'complete_family_summaries_recomputed': len(summaries),
              'strict_matched_order_quality_improvements': all_improvements,
              'independent_pi_grid_bits': bits + 128,
              'independent_taylor_terms': terms,
              'best_fractions': {key: val['best_fraction'] for key, val in summaries.items()},
              'method': 'Polynomial long division and quotient integration; beta monomial sums; exact reduced denominator and gcd; separate outward dyadic Machin checks.',
              'scope': 'Finite Euler/binomial filter experiment. Matched Taylor-baseline gains do not imply new pi records or infinite normalized-quality decay.'}
    (ROOT / 'fresh_verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    run()
