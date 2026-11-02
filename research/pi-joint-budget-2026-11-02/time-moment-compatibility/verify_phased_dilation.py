"""Exact finite matching, filtration, and Pascal checks; standard library only."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import comb, factorial, prod
from pathlib import Path
import json


def ceil_fraction(x):
    return -(-x.numerator//x.denominator)


def determinant(matrix):
    a = [row[:] for row in matrix]
    n = len(a)
    if not n:
        return 1
    old, sign = 1, 1
    for k in range(n-1):
        p = next((i for i in range(k, n) if a[i][k]), None)
        if p is None:
            return 0
        if p != k:
            a[p], a[k] = a[k], a[p]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = a[i][j]*pivot-a[i][k]*a[k][j]
                assert numerator % old == 0
                a[i][j] = numerator//old
            a[i][k] = 0
        old = pivot
    return sign*a[-1][-1]


def phase(N, m, L):
    out = []
    for _ in range(m):
        out.append(N % L)
        N //= L
    return tuple(out)


def pascal_formula(types, L, Q):
    result = Fraction(1)
    for beta in types:
        for b in beta:
            result *= Fraction(prod((L*b)//Q-(L*k)//Q for k in range(b)),
                               factorial(b))
    assert result.denominator == 1 and result >= 1
    return result.numerator


def main():
    H, K, m = 4000, 2, 2
    v0, w0, w = 1, 40, (50, 60)
    theta, L, Q = Fraction(1, 12), 3, 2
    dilation = Fraction(L, Q)
    cutoff = dilation*theta
    v = tuple(int(x/theta) for x in w)
    S = m*(L-1)
    N0 = ceil_fraction(cutoff*H/min(w))+S
    assert w0 > (S+2)*(m+1)*v0
    assert min(v) > v0 and cutoff < 1
    assert {phase(n, m, L) for n in range(L**m)} == set(product(range(L), repeat=m))

    types = [(a, b) for a in range((H-1)//v[0]+1)
             for b in range((H-1-v[0]*a)//v[1]+1)]
    source, targets = [], set()
    buckets = defaultdict(list)
    time_source, time_target, multiplicity = Counter(), Counter(), Counter()
    sum_beta, sum_a = [0]*m, [0]*m
    low_count, outside = 0, 0
    for j in range(K):
        for beta in types:
            for s in range(H-sum(x*y for x, y in zip(v, beta))):
                N = s+sum(beta)
                if j == 0:
                    buckets[N].append(beta)
                if N < N0:
                    a = beta
                    low_count += 1
                else:
                    shift = phase(N, m, L)
                    a = tuple((L*b)//Q+r for b, r in zip(beta, shift))
                t = N-sum(a)
                assert t >= 0 and t+sum(a) == s+sum(beta)
                assert all(aa >= b for aa, b in zip(a, beta))
                assert sum(x*y for x, y in zip(w, a)) <= H
                target = (j, t, a)
                assert target not in targets
                targets.add(target)
                source.append((j, s, beta))
                time_source[s] += 1
                time_target[t] += 1
                multiplicity[a] += 1
                outside += v0*t+sum(x*y for x, y in zip(v, a)) >= H
                for i in range(m):
                    sum_beta[i] += beta[i]
                    sum_a[i] += a[i]
    M = len(source)
    assert len(targets) == M

    # A deterministic cheapest-transverse-degree raw-column choice.
    all_columns = [(w[0]*a+w[1]*b, a, b, h)
                   for a in range(H//w[0]+1)
                   for b in range((H-w[0]*a)//w[1]+1)
                   for h in range((H-w[0]*a-w[1]*b)//w0+1)]
    all_columns.sort()
    assert len(all_columns) >= M
    columns = all_columns[:M]
    counts = Counter((a, b) for _, a, b, _ in columns)
    T = H//w0+1
    count_certificate = K*len(types)*(S+2)*T
    assert count_certificate <= M
    target_low = source_low = 0
    cut_slacks = []
    for r in range(1, T+1):
        target_low += time_target[r-1]
        source_low += time_source[r-1]
        d_r = sum(max(k-r, 0) for k in counts.values())
        C_r = M-source_low-d_r
        E_r = target_low-source_low
        assert 0 <= E_r <= C_r
        assert M-d_r >= Fraction(M*r, T)
        assert target_low <= K*len(types)*(S+2)*r
        cut_slacks.append(C_r-E_r)

    # Verify every distinct changed-bucket Pascal matrix, including its phase.
    seen = set()
    min_pascal = None
    for N, bs in buckets.items():
        bs = tuple(sorted(bs))
        bs_set = set(bs)
        for beta in bs:
            assert all(tuple(g) in bs_set for g in product(*(range(b+1) for b in beta)))
        if N < N0:
            continue
        shift = phase(N, m, L)
        key = (bs, shift)
        if key in seen:
            continue
        seen.add(key)
        nodes = [tuple((L*b)//Q+r for b, r in zip(beta, shift)) for beta in bs]
        matrix = [[prod(comb(a, b) if a >= b else 0 for a, b in zip(node, beta))
                   for node in nodes] for beta in bs]
        value = determinant(matrix)
        assert value == pascal_formula(bs, L, Q) and value > 0
        min_pascal = value if min_pascal is None else min(min_pascal, value)

    column_sums = [sum(col[1+i] for col in columns) for i in range(m)]
    gamma = [a-b for a, b in zip(sum_a, sum_beta)]
    residual = [a-b for a, b in zip(column_sums, sum_a)]
    assert all(g >= 0 for g in gamma) and all(e >= 0 for e in residual)
    mean = Fraction(sum(x*y for x, y in zip(w, sum_a)), M*H)
    collision = Fraction(sum(n*n for n in multiplicity.values()), M*H)
    result = {
        'status': 'passed',
        'parameters': {'H': H, 'm': m, 'K': K, 'v0': v0, 'w0': w0,
                       'w': w, 'theta': str(theta), 'dilation': str(dilation),
                       'cone_cutoff': str(cutoff), 'phase_period': L**m},
        'source_rows': M, 'source_type_count': len(types),
        'ambient_raw_columns': len(all_columns), 'selected_raw_columns': len(columns),
        'target_rows_distinct': len(targets), 'outside_rows': outside,
        'outside_fraction': str(Fraction(outside, M)),
        'low_bucket_threshold': N0, 'low_bucket_rows': low_count,
        'all_cut_counting_left_side': count_certificate,
        'all_cut_counting_right_side': M,
        'nontrivial_cuts_checked': T, 'minimum_target_tail_minus_column_filtration': min(cut_slacks),
        'distinct_pascal_matrices_checked': len(seen), 'minimum_pascal_determinant': min_pascal,
        'gamma': gamma, 'column_minus_target_coordinate_surplus': residual,
        'target_mean': str(mean), 'normalized_collision_cost': str(collision),
        'scope': 'Exact coordinatewise bucket matching, distinctness, all intrinsic raw-column '
                 'time cuts, nonzero Pascal transport minors, and coordinate degree margins. '
                 'No actual-period minor or aggregate size conclusion.'}
    Path(__file__).with_name('phased_dilation_verification.json').write_text(
        json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
