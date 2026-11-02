#!/usr/bin/env python3
"""Exact finite checks for the shifted-column time-grade theorem.

Only the Python standard library is used. Z=1 is a formal specialization of
the period-normalized jet polynomials, not an assertion that exp(1)=1.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, factorial
from pathlib import Path
import json


def matmul(a, b):
    bt = list(zip(*b))
    return [[sum((x*y for x, y in zip(row, col)), F(0)) for col in bt]
            for row in a]


def transpose(a):
    return [list(c) for c in zip(*a)]


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def det(a):
    a = [list(map(F, r)) for r in a]
    n = len(a)
    ans = F(1)
    for i in range(n):
        k = next((k for k in range(i, n) if a[k][i]), None)
        if k is None:
            return F(0)
        if k != i:
            a[i], a[k] = a[k], a[i]
            ans = -ans
        v = a[i][i]
        ans *= v
        for k in range(i+1, n):
            if not a[k][i]:
                continue
            c = a[k][i]/v
            for j in range(i+1, n):
                a[k][j] -= c*a[i][j]
            a[k][i] = F(0)
    return ans


def inverse(a):
    n = len(a)
    out = [list(map(F, r))+e for r, e in zip(a, eye(n))]
    for i in range(n):
        k = next(k for k in range(i, n) if out[k][i])
        out[i], out[k] = out[k], out[i]
        v = out[i][i]
        out[i] = [x/v for x in out[i]]
        for k in range(n):
            if k == i or not out[k][i]:
                continue
            c = out[k][i]
            out[k] = [x-c*y for x, y in zip(out[k], out[i])]
    return [r[n:] for r in out]


def rank(a):
    a = [list(map(F, row)) for row in a]
    k = 0
    for j in range(len(a[0])):
        p = next((p for p in range(k, len(a)) if a[p][j]), None)
        if p is None:
            continue
        a[k], a[p] = a[p], a[k]
        for p in range(k+1, len(a)):
            c = a[p][j]/a[k][j]
            for cidx in range(j, len(a[0])):
                a[p][cidx] -= c*a[k][cidx]
        k += 1
        if k == len(a):
            break
    return k


@lru_cache(None)
def stirling(n, k):
    if n == k == 0:
        return 1
    if n == 0 or k == 0 or k > n:
        return 0
    return k*stirling(n-1, k)+stirling(n-1, k-1)


def jet(row, col, zcenter=F(1)):
    j, t, a = row
    h, alpha = col
    if h > t or a > alpha:
        return F(0)
    return comb(alpha, a)*sum(
        (F(comb(alpha-a, d))*((j*zcenter)**(alpha-a-d))
         *F(factorial(h)*stirling(t-d, h), factorial(t-d))
         for d in range(min(t, alpha-a)+1)), F(0))


def transport(row, target, zeta=F(1, 10)):
    j, s, beta = row
    k, t, a = target
    if j != k or s+beta != t+a or a < beta:
        return F(0)
    return F(comb(a, beta))*(-zeta)**(a-beta)


def case(height):
    # w0=1/4, w1=1, v0=1, v1=4, theta=1/4, K=2.
    ds = [(height-s+3)//4 for s in range(height)]
    source = [(j, s, b) for s in range(height)
              for j in range(2) for b in range(ds[s])]
    columns = [(s, a) for s in range(height) for a in range(2*ds[s])]
    outside = [(j, n-a, a) for j in range(2) for n in range(height)
               for a in range(n+1) if n+3*a >= height]
    n = len(source)
    assert len(columns) == n
    assert all(F(h, 4)+a <= height for h, a in columns)
    delta = sum(s for j, s, b in source)-sum(h for h, a in columns)
    assert delta == 0
    degree = sum(a for h, a in columns)-sum(b for j, s, b in source)
    assert degree == sum(d*d for d in ds)
    b = [[jet(r, c) for c in columns] for r in source]
    bo = [[jet(r, c) for c in columns] for r in outside]
    assert det(b) == 1
    assert det([[jet(r, c, F(2)) for c in columns] for r in source]) == 2**degree
    u = [[transport(r, q) for q in source] for r in source]
    tout = [[transport(r, q) for q in outside] for r in source]
    assert det(u) == 1
    rr = matmul(inverse(u), tout)
    # Independently check the alternating boundary formula.
    boundary_checks = 0
    for i, (j, s, beta) in enumerate(source):
        for k, (jj, t, a) in enumerate(outside):
            value = F(0)
            if jj == j and t+a == s+beta and a >= beta:
                gamma = a-beta
                gap = height-s-4*beta
                alternating = sum((-1)**kk*comb(gamma, kk)
                                  for kk in range(gamma+1) if 3*kk < gap)
                value = F(comb(a, beta))*(-F(1, 10))**gamma*alternating
            assert rr[i][k] == value
            assert not rr[i][k] or t < s
            boundary_checks += 1
    bi = inverse(b)
    x = matmul(bo, bi)
    for i, (j, t, a) in enumerate(outside):
        for k, (jj, ss, bb) in enumerate(source):
            assert not x[i][k] or ss <= t
    e = matmul(rr, bo)
    nn = matmul(e, bi)
    assert any(any(row) for row in nn)
    for i, (j, s, beta) in enumerate(source):
        for k, (jj, ss, bb) in enumerate(source):
            assert not nn[i][k] or ss < s
    power = eye(n)
    nilpotence_index = None
    for k in range(1, height+1):
        power = matmul(power, nn)
        if not any(any(row) for row in power):
            nilpotence_index = k
            break
    assert nilpotence_index is not None
    for eta in (F(-1), F(0), F(1), F(2), F(3, 7)):
        assert det([[bv+eta*ev for bv, ev in zip(br, er)]
                    for br, er in zip(b, e)]) == 1
    # Every coefficient of the marker polynomial follows from strict
    # triangularity above; the five point evaluations are additional checks.
    result = dict(H=height, source_rows=n, outside_rows=len(outside),
                  time_surplus=delta, spatial_degree=degree,
                  boundary_entry_checks=boundary_checks,
                  boundary_rank=rank(e), nilpotence_index=nilpotence_index,
                  marker_checks=5)
    if height == 4:
        buckets = [(j, nn0) for j in range(2) for nn0 in range(height)]
        counts = {'all': 0, 'identity': 0, 'nonidentity_zero': 0}
        for aa in product(*(range(nn0+1) for j, nn0 in buckets)):
            target = [(j, nn0-a, a) for (j, nn0), a in zip(buckets, aa)]
            tt = [[transport(r, q) for q in target] for r in source]
            assert det(tt) != 0
            pb = [[jet(r, c) for c in columns] for r in target]
            counts['all'] += 1
            if set(target) == set(source):
                assert abs(det(pb)) == 1
                counts['identity'] += 1
            else:
                assert det(pb) == 0
                counts['nonidentity_zero'] += 1
        assert counts == {'all': 576, 'identity': 1, 'nonidentity_zero': 575}
        result['exhaustive_targets'] = counts
        grams = []
        for eps in (F(1, 2), F(1, 4), F(1, 8)):
            rt = [[v*eps**(2*(source[i][1]-outside[k][1])-1)
                   if v else F(0) for k, v in enumerate(row)]
                  for i, row in enumerate(rr)]
            xt = [[v*eps**(2*(outside[i][1]-source[k][1])+1)
                   if v else F(0) for k, v in enumerate(row)]
                  for i, row in enumerate(x)]
            rg = matmul(rt, transpose(rt))
            xg = matmul(transpose(xt), xt)
            bound = det([[v+F(i == k) for k, v in enumerate(row)]
                         for i, row in enumerate(rg)])*det(
                            [[v+F(i == k) for k, v in enumerate(row)]
                             for i, row in enumerate(xg)])
            assert bound >= 1
            grams.append(dict(epsilon=str(eps), bound_exact=str(bound),
                              bound_decimal=float(bound)))
        result['midpoint_gram_bounds'] = grams
    return result


def main():
    results = [case(h) for h in (4, 8)]
    output = dict(status='passed', arithmetic='exact rational',
                  formal_center_specializations=[1, 2], cases=results,
                  scope='Finite checks support the written proof. No claim about pi bad approximability.')
    Path(__file__).with_name('time_grade_verification.json').write_text(
        json.dumps(output, indent=2)+'\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
