"""Exact normal/time ideal capacity checks. Python standard library only."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, factorial, prod
from pathlib import Path
import json


def rank(matrix):
    if not matrix:
        return 0
    a = [[F(x) for x in row] for row in matrix]
    r = 0
    for c in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        pivot = a[r][c]
        a[r] = [x / pivot for x in a[r]]
        for i in range(r+1, len(a)):
            if a[i][c]:
                factor = a[i][c]
                a[i] = [x-factor*y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def determinant(matrix):
    a = [row[:] for row in matrix]
    if not a:
        return 1
    old, sign = 1, 1
    for k in range(len(a)-1):
        p = next((i for i in range(k, len(a)) if a[i][k]), None)
        if p is None:
            return 0
        if p != k:
            a[p], a[k] = a[k], a[p]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                numerator = a[i][j]*pivot-a[i][k]*a[k][j]
                assert numerator % old == 0
                a[i][j] = numerator//old
            a[i][k] = 0
        old = pivot
    return sign*a[-1][-1]


def columns(H, w0, w):
    return [(h, a, b) for a in range(H//w[0]+1)
            for b in range((H-w[0]*a)//w[1]+1)
            for h in range((H-w[0]*a-w[1]*b)//w0+1)]


def period_entry(row, col, Z):
    j, t, a, b = row
    h, x, y = col
    if a > x or b > y:
        return 0
    d = x+y-a-b
    return comb(x, a)*comb(y, b)*sum(
        comb(d, k)*(j*Z)**(d-k)*h**(t-k)*factorial(t)//factorial(t-k)
        for k in range(min(d, t)+1))


def coordinate_coefficient(adapted, raw):
    # y=Y-1, T=X1, U=X2-X1; adapted exponent is (time,T,U).
    p, q, nu = adapted
    h, a, b = raw
    return comb(h, p)*comb(b, nu) if p <= h and nu <= b and q == a+b-nu else 0


def small_checks():
    checks = 0
    records = []
    for H, w0, w in ((3, 1, (1, 2)), (4, 2, (1, 2)), (5, 2, (2, 3))):
        ambient = columns(H, w0, w)
        transform = [[coordinate_coefficient(ad, raw) for raw in ambient] for ad in ambient]
        assert rank(transform) == len(ambient)
        chosen = ambient[::2]
        targets = [(j, t, a, b) for j in range(2) for t in range(2)
                   for a, b in ((0, 0), (1, 0), (0, 1))]
        for Z in (0, 1, 2):
            actual_rank = rank([[period_entry(row, col, Z) for col in chosen] for row in targets])
            for k in range(4):
                for r in range(4):
                    outside = [ad for ad in ambient if not (ad[0] >= k and ad[2] >= r)]
                    quotient = [[coordinate_coefficient(ad, raw) for raw in chosen] for ad in outside]
                    intersection_dimension = len(chosen)-rank(quotient)
                    survivors = sum(t >= k and a+b >= r for j, t, a, b in targets)
                    full_ideal = sum(h >= k and nu >= r for h, q, nu in ambient)
                    assert intersection_dimension >= max(0, len(chosen)-len(ambient)+full_ideal)
                    assert actual_rank <= len(chosen)-intersection_dimension+survivors
                    checks += 1
        records.append({'H': H, 'w0': w0, 'w': w, 'ambient_columns': len(ambient),
                        'selected_columns': len(chosen), 'target_rows': len(targets)})
    return {'exact_selected_ideal_rank_checks': checks, 'models': records}


def suffix_grid(counts, xmax, ymax):
    grid = [[0]*(ymax+2) for _ in range(xmax+2)]
    for (x, y), value in counts.items():
        grid[x][y] += value
    for x in range(xmax, -1, -1):
        for y in range(ymax, -1, -1):
            grid[x][y] += grid[x+1][y]+grid[x][y+1]-grid[x+1][y+1]
    return grid


def phase_counterexample():
    H, K, w0, w, v0, v, L, Q, N0 = 2000, 2, 24, (40, 40), 1, (400, 400), 3, 2, 12
    source, target = [], []
    buckets = defaultdict(list)
    for j in range(K):
        for a in range(H//v[0]+1):
            for b in range(H//v[1]+1):
                for s in range(max(0, H-v[0]*a-v[1]*b)):
                    n = s+a+b
                    aa, bb = (a, b) if n < N0 else (L*a//Q+n%L, L*b//Q+(n//L)%L)
                    t = n-aa-bb
                    assert t >= 0 and aa >= a and bb >= b
                    assert w[0]*aa+w[1]*bb <= H
                    source.append((j, s, a, b)); target.append((j, t, aa, bb))
                    if j == 0:
                        buckets[n].append((a, b))
    assert len(set(source)) == len(set(target)) == len(source)
    M = len(source)
    ambient = columns(H, w0, w)
    chosen = sorted(ambient, key=lambda c: (w[0]*c[1]+w[1]*c[2], c[1], c[2], c[0]))[:M]
    assert len(chosen) == M
    multiplicities = Counter((a, b) for h, a, b in chosen)
    source_times = Counter(row[1] for row in source)
    target_times = Counter(row[1] for row in target)
    T = H//w0+1
    source_type_count = len({(a, b) for j, s, a, b in source})
    assert K*source_type_count*6*T <= M
    slow = tlow = 0
    cuts = []
    for r in range(1, T+1):
        slow += source_times[r-1]; tlow += target_times[r-1]
        dr = sum(max(n-r, 0) for n in multiplicities.values())
        assert 0 <= tlow-slow <= M-slow-dr
        cuts.append(M-dr-tlow)
    gamma = [sum(row[2+i] for row in target)-sum(row[2+i] for row in source) for i in range(2)]
    residual = [sum(c[1+i] for c in chosen)-sum(row[2+i] for row in target) for i in range(2)]
    assert min(gamma+residual) >= 0
    seen = set()
    for n, types in buckets.items():
        if n < N0:
            continue
        types = tuple(sorted(types)); phase = (n%L, (n//L)%L)
        if (types, phase) in seen:
            continue
        seen.add((types, phase))
        nodes = [(L*a//Q+phase[0], L*b//Q+phase[1]) for a, b in types]
        pascal = [[comb(x, a)*comb(y, b) if x >= a and y >= b else 0
                   for x, y in nodes] for a, b in types]
        formula = F(1)
        for beta in types:
            for b in beta:
                formula *= F(prod(L*b//Q-L*k//Q for k in range(b)), factorial(b))
        assert determinant(pascal) == formula > 0
    amax = max(a+b for j, t, a, b in target)
    visible = sum(nu <= amax for h, q, nu in ambient)
    assert visible < M
    col_counts = Counter((h, nu) for h, q, nu in ambient)
    row_counts = Counter((t, a+b) for j, t, a, b in target)
    cs = suffix_grid(col_counts, H//w0, H//w[1])
    ts = suffix_grid(row_counts, max(t for j, t, a, b in target), amax)
    best = (len(ambient), None)
    mixed_checks = 0
    for k in range(H//w0+2):
        for r in range(H//w[1]+2):
            surviving = ts[k][r] if r <= amax else 0
            bound = len(ambient)-cs[k][r]+surviving
            if bound < best[0]:
                best = (bound, (k, r))
            mixed_checks += 1
    return {'parameters': {'H': H, 'K': K, 'm': 2, 'v0': v0, 'v': v,
                            'w0': w0, 'w': w, 'theta': '1/10', 'dilation': '3/2', 'low_threshold': N0},
            'source_rows': M, 'target_rows': len(target), 'ambient_columns': len(ambient),
            'selected_columns': len(chosen), 'gamma': gamma, 'residual': residual,
            'time_cuts': len(cuts), 'minimum_time_cut_slack': min(cuts),
            'pascal_blocks_checked': len(seen), 'max_target_transverse_order': amax,
            'visible_column_bound': visible, 'rank_deficiency_at_least': M-visible,
            'mixed_ideal_cuts_checked': mixed_checks, 'best_mixed_rank_bound': best[0],
            'best_mixed_cut': best[1], 'best_rank_deficiency_at_least': M-best[0],
            'normal_ninth_power_quotient_dimension': sum(nu < 9 for h, q, nu in ambient),
            'normal_ninth_power_surviving_rows': sum(a+b >= 9 for j, t, a, b in target),
            'scope': 'All time cuts and transport pass, but every target period minor is zero by exact normal-ideal capacity. The equal-weight source is itself rank deficient; this is not a counterexample in the original separated-weight rank-eligible regime.'}


if __name__ == '__main__':
    result = {'status': 'passed', 'small_checks': small_checks(),
              'phase_counterexample': phase_counterexample()}
    Path(__file__).with_name('mixed_normal_cut_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
