"""Independent finite verification of the joint-selection experiment.

No import of experiment.py or either original constructor. Apéry sources are
rebuilt from binomial and telescoping sums. Quadratic-root convergents are
extracted from strict dyadic intervals, not the constructor's surd recurrence.
"""
from fractions import Fraction as F
from math import comb, gcd, isqrt
from pathlib import Path
import argparse
import hashlib
import json
import sys

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parent


def sources(nmax):
    u = [sum(comb(n, j) ** 2 * comb(n + j, j) for j in range(n + 1))
         for n in range(nmax + 2)]
    ratio = F(0)
    r = [F(0)]
    for n in range(1, nmax + 2):
        ratio += F(5 * (-1) ** (n - 1), n * n * u[n] * u[n - 1])
        r.append(6 * ratio)
    assert r[1] == 10
    return u, r


def rational_cf(value):
    """Euclidean quotients of an exact rational number."""
    num, den = value.numerator, value.denominator
    while den:
        a, num = divmod(num, den)
        yield a
        num, den = den, num


def bezout(a, b):
    # Backward substitution in the Euclidean algorithm.
    steps = []
    while b:
        q, rem = divmod(a, b)
        steps.append(q)
        a, b = b, rem
    x, y = 1, 0
    for q in reversed(steps):
        x, y = y, x - q * y
    return a, x, y


def source_pool(r0, r1, pool='original'):
    A0, B0, A1, B1 = r0.numerator, r0.denominator, r1.numerator, r1.denominator
    G = abs(A0 * B1 - A1 * B0)
    g, r, s = bezout(A1, B1)
    assert g == 1 and r * A1 + s * B1 == 1
    k = -(r * A0 + s * B0) % G
    assert (A0 + k * A1) % G == (B0 + k * B1) % G == 0
    weights = {(1, 0), (0, 1)}
    columns = [(0, 1), (1, 0)]
    if k:
        for a in rational_cf(F(k, G)):
            js = {1, a, max(1, a - 1)}
            if pool != 'original':
                raise ValueError('Unknown weight pool: ' + pool)
            old, new = columns[-2], columns[-1]
            for j in js:
                if j > a:
                    continue
                p, q = old[0] + j * new[0], old[1] + j * new[1]
                x, y = q, k * q - G * p
                if x > 0 and y >= 0:
                    d = gcd(x, y)
                    weights.add((x // d, y // d))
            columns.append((old[0] + a * new[0], old[1] + a * new[1]))
    candidates = []
    for x, y in weights:
        aa, bb = x * A0 + y * A1, x * B0 + y * B1
        content = gcd(aa, bb)
        assert G % content == 0 and gcd(x, y) == 1
        A, B = aa // content, bb // content
        assert min(r0, r1) <= F(A, B) <= max(r0, r1)
        candidates.append((B, A, x, y, content))
    result = sorted(candidates)
    assert len({(c[0], c[1]) for c in result}) == len(result)
    return result


def root_convergents(A, B, qmax):
    """Exact Euclidean CF from dyadic brackets, including rational-square case."""
    root = isqrt(A * B)
    pm2, pm1, qm2, qm1 = 0, 1, 1, 0
    if root * root == A * B:
        for a in rational_cf(F(root, B)):
            p, q = a * pm1 + pm2, a * qm1 + qm2
            if q > qmax:
                return
            yield p, q
            pm2, pm1, qm2, qm1 = pm1, p, qm1, q
        return
    bits = 4 * qmax.bit_length() + 2 * B.bit_length() + 128
    scale = 1 << bits
    lower = isqrt(A * scale * scale // B)
    ln, ld, hn, hd = lower, scale, lower + 1, scale
    assert ln * ln * B < A * ld * ld
    assert hn * hn * B > A * hd * hd
    while True:
        a = ln // ld
        upper_floor = (hn - 1) // hd
        assert a == upper_floor, 'Insufficient dyadic precision; no guessed branch.'
        p, q = a * pm1 + pm2, a * qm1 + qm2
        if q > qmax:
            return
        yield p, q
        pm2, pm1, qm2, qm1 = pm1, p, qm1, q
        assert ln > a * ld
        ln, ld, hn, hd = hd, hn - a * hd, ld, ln - a * ld


def atan_interval(t, bits):
    """Taylor terms rounded outward on an integer dyadic grid."""
    scale = 1 << bits
    lo = hi = 0
    power = t
    j = 0
    while True:
        low, rem = divmod(scale, (2 * j + 1) * power)
        if low == 0:
            if j % 2:
                lo -= 1
            else:
                hi += 1
            return F(lo, scale), F(hi, scale), j
        high = low + bool(rem)
        if j % 2:
            lo -= high
            hi -= low
        else:
            lo += low
            hi += high
        j += 1
        power *= t * t


def pi_interval(bits):
    # Exact Gaussian identity fixes the Machin angle relation; positive angle
    # 4 atan(1/5)-atan(1/239) is less than 4/5 < pi/2.
    z = (1, 0)
    for _ in range(4):
        z = (5 * z[0] - z[1], z[0] + 5 * z[1])
    assert z == (476, 480) == (2 * 239 - 2, 2 * 239 + 2)
    assert 4 * (F(1, 5) - F(1, 375)) - F(1, 239) > 0
    l5, h5, t5 = atan_interval(5, bits)
    l239, h239, t239 = atan_interval(239, bits)
    lo, hi = 16 * l5 - 4 * h239, 16 * h5 - 4 * l239
    assert 3 < lo < hi < 4
    return lo, hi, [t5, t239]


def quality_interval(p, q, lo, hi):
    r = F(p, q)
    assert r < lo or r > hi
    if r < lo:
        return q * (q * lo - p), q * (q * hi - p)
    return q * (p - q * hi), q * (p - q * lo)


def pool_hash(pool):
    raw = [[str(v) for v in row] for row in pool]
    return hashlib.sha256(json.dumps(raw, separators=(',', ':')).encode()).hexdigest()


def best_bound(indices, hits):
    if not indices:
        return None
    return min(indices, key=lambda i: (F(hits[i]['quality_upper_from_conic']), i))


def best_quality(indices, hits):
    if not indices:
        return None
    winner = min(indices, key=lambda i: (sum(map(F, hits[i]['quality_interval'])),
                                         F(hits[i]['quality_upper_from_conic']), i))
    pq = (hits[winner]['p'], hits[winner]['q'])
    for i in indices:
        if (hits[i]['p'], hits[i]['q']) != pq:
            assert F(hits[winner]['quality_interval'][1]) < F(hits[i]['quality_interval'][0])
    return winner


def run(filename='results.json', output='verification.json'):
    path = ROOT / filename
    data = json.loads(path.read_text())
    par = data['parameters']
    assert data['schema'] == 1
    assert data['base_commit'] == '88c331c410290c5eb915776ecd2415bffff20dc0'
    assert par['n_min'] == 1 and par['weight_pool'] == 'original'
    assert par['budget_names'] == ['fixed25', 'square_root']
    nmax = par['n_max']
    u, rs = sources(nmax)
    assert len(data['cases']) == nmax
    widths = [F(30, (n + 1) ** 2 * u[n] * u[n + 1]) for n in range(1, nmax + 1)]
    bits = max(d.denominator.bit_length() - d.numerator.bit_length() for d in widths) + 512
    plo, phi, terms = pi_interval(bits)
    slo, shi = map(F, data['pi_scoring_interval'])
    assert slo <= plo < phi <= shi
    total_sources = total_roots = baseline_roots = 0
    expected_hits = []
    for n, case in enumerate(data['cases'], 1):
        assert case['n'] == n
        r0, r1 = rs[n], rs[n + 1]
        delta = widths[n - 1]
        assert abs(r0 - r1) == delta == F(case['source_width'])
        assert min(r0, r1) < plo * plo < phi * phi < max(r0, r1)
        pool = source_pool(r0, r1, par['weight_pool'])
        assert len(pool) == case['source_count']
        assert pool_hash(pool) == case['pool_sha256']
        assert pool[0] == tuple(map(int, case['baseline_source']))
        baseline = pool[0]
        checked = base_checked = hit_count = base_hit_count = 0
        for row in pool:
            B, A, x, y, g = row
            budgets = {'fixed25': min(25, B), 'square_root': isqrt(B)}
            limits = {name: isqrt(H * delta.denominator // (B * delta.numerator))
                      for name, H in budgets.items()}
            for p, q in root_convergents(A, B, max(limits.values())):
                if not 3 * q <= p <= 4 * q:
                    continue
                checked += 1
                base_checked += row == baseline
                norm = B * p * p - A * q * q
                for name, H in budgets.items():
                    if q > limits[name] or abs(norm) > H:
                        continue
                    assert gcd(p, q) == 1 and B * q * q * delta <= H
                    bound = (F(abs(norm), B) + q * q * delta) / 6
                    expected_hits.append((n, row, name, p, q, norm, H, limits[name], bound))
                    hit_count += 1
                    base_hit_count += row == baseline
        assert checked == case['root_convergents_checked']
        assert base_checked == case['baseline_root_convergents_checked']
        assert hit_count == case['hits_with_budgets']
        assert base_hit_count == case['baseline_hits_with_budgets']
        total_sources += len(pool)
        total_roots += checked
        baseline_roots += base_checked
        if n % 32 == 0:
            print(json.dumps({'verified_n': n, 'sources': total_sources, 'root_convergents': total_roots}), flush=True)
    hits = data['hits']
    assert len(expected_hits) == len(hits)
    for expected, h in zip(expected_hits, hits):
        n, row, name, p, q, norm, H, limit, bound = expected
        assert (n, row, name) == (h['n'], tuple(map(int, h['source'])), h['budget'])
        assert [p, q, norm, H, limit] == [int(h[k]) for k in ['p', 'q', 'norm', 'H', 'q_limit']]
        assert F(h['quality_upper_from_conic']) == bound
        ql, qh = map(F, h['quality_interval'])
        dl, du = map(F, h['quality_decimal_enclosure'])
        tl, tu = quality_interval(p, q, plo, phi)
        assert 0 < dl <= ql <= tl < tu <= qh <= du
        assert qh < bound
        assert (ql, qh) == quality_interval(p, q, slo, shi)
    base_indices = []
    improved_bound = []
    improved_quality = []
    newly_successful = []
    for case in data['cases']:
        inds = [i for i, h in enumerate(hits) if h['n'] == case['n']]
        base = [i for i in inds if hits[i]['source'] == case['baseline_source']]
        base_indices.extend(base)
        for key, value in [('joint_best_bound', best_bound(inds, hits)),
                           ('baseline_best_bound', best_bound(base, hits)),
                           ('joint_best_quality', best_quality(inds, hits)),
                           ('baseline_best_quality', best_quality(base, hits))]:
            assert case[key] == value
        j, b = case['joint_best_bound'], case['baseline_best_bound']
        if j is not None and b is None:
            newly_successful.append(case['n'])
        elif b is not None and F(hits[j]['quality_upper_from_conic']) < F(hits[b]['quality_upper_from_conic']):
            improved_bound.append(case['n'])
        j, b = case['joint_best_quality'], case['baseline_best_quality']
        if b is not None and (hits[j]['p'], hits[j]['q']) != (hits[b]['p'], hits[b]['q']):
            improved_quality.append(case['n'])
    distinct = sorted({(int(h['p']), int(h['q'])) for h in hits}, key=lambda t: (t[1], t[0]))
    bd = sorted({(int(hits[i]['p']), int(hits[i]['q'])) for i in base_indices}, key=lambda t: (t[1], t[0]))
    all_indices = list(range(len(hits)))
    expected_summary = {
        'source_count': total_sources, 'root_convergents_checked': total_roots,
        'baseline_source_count': nmax, 'baseline_root_convergents_checked': baseline_roots,
        'hits_with_budgets': len(hits), 'baseline_hits_with_budgets': len(base_indices),
        'distinct_fractions': [f'{p}/{q}' for p, q in distinct],
        'baseline_distinct_fractions': [f'{p}/{q}' for p, q in bd],
        'new_to_baseline_fractions': [f'{p}/{q}' for p, q in distinct if (p, q) not in bd],
        'newly_successful_n': newly_successful, 'improved_bound_n': improved_bound,
        'improved_actual_quality_n': improved_quality,
        'max_hit_n': max((h['n'] for h in hits), default=None),
        'baseline_max_hit_n': max((hits[i]['n'] for i in base_indices), default=None),
        'joint_best_bound': best_bound(all_indices, hits),
        'baseline_best_bound': best_bound(base_indices, hits),
        'joint_best_quality': best_quality(all_indices, hits),
        'baseline_best_quality': best_quality(base_indices, hits)}
    assert data['summary'] == expected_summary
    report = {
        'checks': 'PASS', 'file': path.name,
        'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'source_indices': nmax, 'source_pools_rebuilt': nmax,
        'primitive_sources_rebuilt': total_sources,
        'independent_root_convergents_checked': total_roots,
        'baseline_root_convergents_checked': baseline_roots,
        'certified_hits': len(hits), 'summary_recomputed_exactly': True,
        'independent_pi_grid_bits': bits, 'independent_taylor_terms': terms,
        'distinct_fractions': expected_summary['distinct_fractions'],
        'improved_actual_quality_n': improved_quality,
        'method': 'Separate implementation: binomial/telescoping sources, exact dyadic root brackets, complete weight pools, integer norm budgets, outward Taylor scoring.',
        'scope': 'All listed finite cases and comparisons verified. No global weight optimum or infinite-success theorem.'}
    (ROOT / output).write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', default='results.json')
    parser.add_argument('--output', default='verification.json')
    args = parser.parse_args()
    run(args.input, args.output)
