"""Exact finite ablation: keep the old weight pool, change only selection.

The old algorithm chooses the primitive pi^2 source with smallest denominator
before its small-norm search. Here every source in precisely that finite pool
is searched with precisely the old two norm/tail budgets. Pi is used only after
candidate generation, to score accepted fractions with rational enclosures.

Run: python3 experiment.py
Standard library only. Outputs are deterministic (no timing fields).
"""
from fractions import Fraction as F
from math import gcd, isqrt
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / 'pi-small-norm-construction-2026-11-04'
sys.path.insert(0, str(OLD))
from conic_attempt import source, root_convergents
from content_lift import bezout
from construct import display_interval
from verify import atan_taylor

sys.set_int_max_str_digits(0)
NMAX = 128


def weight_pool(A0, B0, A1, B1):
    """Exactly the original content_lift.weights pool, without its min()."""
    G = abs(A0 * B1 - A1 * B0)
    common, r, s = bezout(A1, B1)
    assert common == 1 and G > 0
    k = (-(r * A0 + s * B0)) % G
    assert (A0 + k * A1) % G == (B0 + k * B1) % G == 0
    choices = {(1, 0), (0, 1)}
    if k:
        num, den = k, G
        pm2, pm1, qm2, qm1 = 0, 1, 1, 0
        while den:
            a, rem = divmod(num, den)
            for j in {1, a, max(1, a - 1)}:
                if j > a:
                    continue
                p, q = j * pm1 + pm2, j * qm1 + qm2
                x, y = q, k * q - G * p
                if x > 0 and y >= 0:
                    common = gcd(x, y)
                    choices.add((x // common, y // common))
            p, q = a * pm1 + pm2, a * qm1 + qm2
            pm2, pm1, qm2, qm1 = pm1, p, qm1, q
            num, den = den, rem
    pool = []
    for x, y in choices:
        aa, bb = x * A0 + y * A1, x * B0 + y * B1
        g = gcd(aa, bb)
        A, B = aa // g, bb // g
        assert gcd(x, y) == gcd(A, B) == 1 and G % g == 0
        pool.append((B, A, x, y, g))
    pool.sort()
    assert len({(A, B) for B, A, x, y, g in pool}) == len(pool)
    return pool


def pool_hash(pool):
    raw = json.dumps([list(map(str, t)) for t in pool], separators=(',', ':'))
    return hashlib.sha256(raw.encode()).hexdigest()


def best_bound(indices, hits):
    return min(indices, key=lambda i: (F(hits[i]['quality_upper_from_conic']), i)) if indices else None


def best_quality(indices, hits):
    if not indices:
        return None
    chosen = min(indices, key=lambda i: (sum(map(F, hits[i]['quality_interval'])),
                                        F(hits[i]['quality_upper_from_conic']), i))
    same = (hits[chosen]['p'], hits[chosen]['q'])
    upper = F(hits[chosen]['quality_interval'][1])
    for i in indices:
        if (hits[i]['p'], hits[i]['q']) != same:
            assert upper < F(hits[i]['quality_interval'][0]), 'Increase scoring precision; do not guess a rank.'
    return chosen


def run():
    u, v = source(NMAX + 1)
    cases, hits = [], []
    for n in range(1, NMAX + 1):
        r0, r1 = 6 * v[n] / u[n], 6 * v[n + 1] / u[n + 1]
        delta = abs(r0 - r1)
        pool = weight_pool(r0.numerator, r0.denominator, r1.numerator, r1.denominator)
        baseline = pool[0]
        case = {'n': n, 'source_count': len(pool), 'pool_sha256': pool_hash(pool),
                'baseline_source': list(map(str, baseline)),
                'source_width': str(delta), 'root_convergents_checked': 0,
                'baseline_root_convergents_checked': 0,
                'hits_with_budgets': 0, 'baseline_hits_with_budgets': 0}
        for record in pool:
            B, A, x, y, g = record
            budgets = {'fixed25': min(25, B), 'square_root': isqrt(B)}
            limits = {name: isqrt((H * delta.denominator) // (B * delta.numerator))
                      for name, H in budgets.items()}
            for p, q in root_convergents(A, B, max(limits.values())):
                if not 3 * q <= p <= 4 * q:
                    continue
                case['root_convergents_checked'] += 1
                if record == baseline:
                    case['baseline_root_convergents_checked'] += 1
                norm = B * p * p - A * q * q
                for name, H in budgets.items():
                    if q > limits[name] or abs(norm) > H:
                        continue
                    assert gcd(p, q) == 1 and B * q * q * delta <= H
                    bound = (F(abs(norm), B) + q * q * delta) / 6
                    hits.append({'n': n, 'source': list(map(str, record)), 'budget': name,
                                 'p': str(p), 'q': str(q), 'norm': str(norm),
                                 'H': str(H), 'q_limit': str(limits[name]),
                                 'quality_upper_from_conic': str(bound)})
                    case['hits_with_budgets'] += 1
                    if record == baseline:
                        case['baseline_hits_with_budgets'] += 1
        cases.append(case)
        if n % 16 == 0:
            print(json.dumps({'completed_n': n, 'sources': sum(c['source_count'] for c in cases),
                              'checked': sum(c['root_convergents_checked'] for c in cases),
                              'hits': len(hits)}), flush=True)

    # This is the first access to a pi enclosure; all candidates already exist.
    scoring_bits = 2 * max((int(h['q']).bit_length() for h in hits), default=1) + 256
    l5, h5, _ = atan_taylor(5, scoring_bits)
    l239, h239, _ = atan_taylor(239, scoring_bits)
    lo, hi = 16 * l5 - 4 * h239, 16 * h5 - 4 * l239
    assert 3 < lo < hi < 4
    for h in hits:
        p, q = int(h['p']), int(h['q'])
        if F(p, q) < lo:
            ql, qh = q * (q * lo - p), q * (q * hi - p)
        else:
            assert F(p, q) > hi
            ql, qh = q * (p - q * hi), q * (p - q * lo)
        assert 0 < ql < qh < F(h['quality_upper_from_conic'])
        h['quality_interval'] = [str(ql), str(qh)]
        h['quality_decimal_enclosure'] = display_interval(ql, qh)

    base_indices = []
    for c in cases:
        joint = [i for i, h in enumerate(hits) if h['n'] == c['n']]
        base = [i for i in joint if hits[i]['source'] == c['baseline_source']]
        base_indices.extend(base)
        c['joint_best_bound'] = best_bound(joint, hits)
        c['baseline_best_bound'] = best_bound(base, hits)
        c['joint_best_quality'] = best_quality(joint, hits)
        c['baseline_best_quality'] = best_quality(base, hits)
    all_indices = list(range(len(hits)))
    distinct = sorted({(int(h['p']), int(h['q'])) for h in hits}, key=lambda t: (t[1], t[0]))
    base_distinct = sorted({(int(hits[i]['p']), int(hits[i]['q'])) for i in base_indices}, key=lambda t: (t[1], t[0]))
    improved_bound, improved_quality, newly_successful = [], [], []
    for c in cases:
        j, b = c['joint_best_bound'], c['baseline_best_bound']
        if j is not None and b is None:
            newly_successful.append(c['n'])
        elif b is not None and F(hits[j]['quality_upper_from_conic']) < F(hits[b]['quality_upper_from_conic']):
            improved_bound.append(c['n'])
        j, b = c['joint_best_quality'], c['baseline_best_quality']
        if b is not None and (hits[j]['p'], hits[j]['q']) != (hits[b]['p'], hits[b]['q']):
            assert F(hits[j]['quality_interval'][1]) < F(hits[b]['quality_interval'][0])
            improved_quality.append(c['n'])
    summary = {
        'source_count': sum(c['source_count'] for c in cases),
        'root_convergents_checked': sum(c['root_convergents_checked'] for c in cases),
        'baseline_source_count': len(cases),
        'baseline_root_convergents_checked': sum(c['baseline_root_convergents_checked'] for c in cases),
        'hits_with_budgets': len(hits), 'baseline_hits_with_budgets': len(base_indices),
        'distinct_fractions': [f'{p}/{q}' for p, q in distinct],
        'baseline_distinct_fractions': [f'{p}/{q}' for p, q in base_distinct],
        'new_to_baseline_fractions': [f'{p}/{q}' for p, q in distinct if (p, q) not in base_distinct],
        'newly_successful_n': newly_successful, 'improved_bound_n': improved_bound,
        'improved_actual_quality_n': improved_quality,
        'max_hit_n': max((h['n'] for h in hits), default=None),
        'baseline_max_hit_n': max((hits[i]['n'] for i in base_indices), default=None),
        'joint_best_bound': best_bound(all_indices, hits),
        'baseline_best_bound': best_bound(base_indices, hits),
        'joint_best_quality': best_quality(all_indices, hits),
        'baseline_best_quality': best_quality(base_indices, hits)}
    out = {'schema': 1, 'base_commit': '88c331c410290c5eb915776ecd2415bffff20dc0',
           'parameters': {'n_min': 1, 'n_max': NMAX, 'weight_pool': 'original',
                          'budget_names': ['fixed25', 'square_root'],
                          'score_pi_bits': scoring_bits},
           'pi_scoring_interval': [str(lo), str(hi)], 'cases': cases, 'hits': hits,
           'summary': summary,
           'scope': 'Finite same-pool, same-budget selection ablation. No global weight optimum, new pi record, or infinite-success claim.'}
    (ROOT / 'results.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(summary, indent=2), flush=True)
    return out


if __name__ == '__main__':
    run()
