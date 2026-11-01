"""Exact standard-library checks for the small-escape Hardy Gram bound.

The kernels use real nodes 0,6,...: unitary rotation has the same Gram
determinants as imaginary nodes. These are rational surrogate nodes, not
an assertion that exp(6i)=1. The written proof handles the actual period.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, prod
from pathlib import Path
import json


def determinant(a):
    if not a:
        return F(1)
    a = [[F(x) for x in r] for r in a]
    out = F(1)
    for k in range(len(a)):
        p = next((i for i in range(k, len(a)) if a[i][k]), None)
        if p is None:
            return F(0)
        if p != k:
            a[k], a[p] = a[p], a[k]
            out = -out
        diagonal = a[k][k]
        out *= diagonal
        for i in range(k+1, len(a)):
            ratio = a[i][k]/diagonal
            for j in range(k+1, len(a)):
                a[i][j] -= ratio*a[k][j]
    return out


def inverse(a):
    n = len(a)
    a = [[F(x) for x in r]+[F(i == j) for j in range(n)]
         for i, r in enumerate(a)]
    for k in range(n):
        p = next(i for i in range(k, n) if a[i][k])
        a[k], a[p] = a[p], a[k]
        diagonal = a[k][k]
        a[k] = [x/diagonal for x in a[k]]
        for i in range(n):
            if i == k:
                continue
            factor = a[i][k]
            if factor:
                a[i] = [x-factor*y for x, y in zip(a[i], a[k])]
    return [r[n:] for r in a]


def is_psd(a):
    """Exact diagonal-pivot Schur algorithm for real symmetric matrices."""
    a = [[F(x) for x in r] for r in a]
    while a:
        if any(a[i][i] < 0 for i in range(len(a))):
            return False
        p = next((i for i in range(len(a)) if a[i][i] > 0), None)
        if p is None:
            return all(x == 0 for row in a for x in row)
        ids = [i for i in range(len(a)) if i != p]
        a = [[a[i][j]-a[i][p]*a[p][j]/a[p][p] for j in ids] for i in ids]
    return True


def kernel(x, y, s, t, radius):
    radius2 = F(radius)**2
    base = 1-F(x*y)/radius2
    return sum(F(factorial(s+t-c), factorial(s-c)*factorial(t-c)*factorial(c))
               * (F(y)/radius2)**(s-c)*(F(x)/radius2)**(t-c)
               * (1/radius2)**c/base**(s+t-c+1)
               for c in range(min(s, t)+1))


def gram(rows, radius):
    return [[kernel(6*j, 6*l, s, t, radius) if a == b else F(0)
             for l, t, b in rows] for j, s, a in rows]


def mul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def power(a, n):
    out = [F(1)]
    for _ in range(n):
        out = mul(out, a)
    return out


def shifted(p, center):
    return [sum(p[d]*comb(d, s)*center**(d-s) for d in range(s, len(p)))
            for s in range(len(p))]


def cardinal(K, n, j, s):
    x = 6*j
    lam = [F(1)]
    for l in range(K):
        if l != j:
            lam = mul(lam, power([F(-6*l, x-6*l), F(1, x-6*l)], n))
    local = shifted(lam, x)
    assert local[0] == 1
    reciprocal = [F(1)]
    for d in range(1, n-s):
        reciprocal.append(-sum(local[e]*reciprocal[d-e]
                               for e in range(1, min(d, len(local)-1)+1)))
    correction = [F(0)]*(n-s)
    for l, value in enumerate(reciprocal):
        term = power([-x, 1], l)
        for d, coefficient in enumerate(term):
            correction[d] += value*coefficient
    return mul(mul(lam, power([-x, 1], s)), correction)


def q_bound(K):
    L = F(142*K-44, 7)
    return L*F((K+1)**2, 4)*(L/3)**(K-1)/factorial(K-1)


def check_cardinals():
    records = []
    for K in (1, 2, 3):
        R = 14*K
        for n in (1, 2, 3, 4):
            rows = [(j, s, 0) for j in range(K) for s in range(n)]
            polynomials = [cardinal(K, n, j, s) for j, s, _ in rows]
            for i, p in enumerate(polynomials):
                for r, (j, s, _) in enumerate(rows):
                    value = shifted(p, 6*j)[s] if s < len(p) else F(0)
                    assert value == (i == r)
            weighted = [[(p[d] if d < len(p) else 0)*R**d
                         for d in range(K*n)] for p in polynomials]
            pg = [[sum(x*y for x, y in zip(a, b)) for b in weighted]
                  for a in weighted]
            g = gram(rows, R)
            gi = inverse(g)
            diff = [[pg[i][j]-gi[i][j] for j in range(K*n)] for i in range(K*n)]
            assert is_psd(diff)
            assert all(pg[i][i] <= n*n*q_bound(K)**(2*n) for i in range(K*n))
            trace_inverse = sum(gi[i][i] for i in range(K*n))
            assert trace_inverse <= K*n**3*q_bound(K)**(2*n)
            records.append({'K': K, 'n': n, 'cardinal_polynomials': K*n,
                            'inverse_operator_majorization': 'exact PSD check passed',
                            'trace_inverse': str(trace_inverse),
                            'Q_K': str(q_bound(K))})
    return records


def check_exchanges():
    # Complete source packet H=3, v0=1, v1=2, theta=1/2, w1=1.
    K, H, R = 2, 3, 28
    source = [(j, s, b) for j in range(K) for b in range(H)
              for s in range(H) if s+2*b < H]
    ambient = [(j, s, a) for j in range(K) for a in range(H)
               for s in range(H) if s+a < H and a <= H]
    outside = [r for r in ambient if r not in source]
    assert (len(source), len(outside)) == (8, 4)
    g = gram(source, R)
    gi = inverse(g)
    gd = determinant(g)
    assert gd > 0
    bound = F(9604, 4563)*K*H**3*q_bound(K)**(2*H)
    count = 0
    rank_counts = {}
    for k in (1, 2, 3, 4):
        removals = list(combinations(range(len(source)), k))
        # All one- and two-exchange choices; representative higher exchanges.
        if k >= 3:
            removals = removals[:4]
        for removed in removals:
            kept = [r for i, r in enumerate(source) if i not in removed]
            invminor = determinant([[gi[i][j] for j in removed] for i in removed])
            assert determinant(gram(kept, R)) == gd*invminor
            for additions in combinations(outside, k):
                target = kept+list(additions)
                target_determinant = determinant(gram(target, R))
                product_diagonals = prod(kernel(6*j, 6*j, s, s, R)
                                         for j, s, a in additions)
                assert all(kernel(6*j, 6*j, s, s, R) <= F(9604, 4563)
                           for j, s, a in additions)
                assert 0 < target_determinant <= gd*invminor*product_diagonals
                assert target_determinant <= gd*bound**k
                count += 1
                rank_counts[k] = rank_counts.get(k, 0)+1
    return {'complete_source_rows': len(source), 'ambient_rows': len(ambient),
            'source_gram_determinant': str(gd), 'exchange_cases': count,
            'counts_by_number_of_exchanges': rank_counts,
            'source_parameters': {'K': K, 'H': H, 'v0': 1, 'v1': 2,
                                  'theta': '1/2', 'w1': 1},
            'kernel_nodes': [0, 6], 'radius': R}


def main():
    records = check_cardinals()
    exchanges = check_exchanges()
    result = {'status': 'passed', 'hermite_inverse_cases': records,
              'minor_exchange_checks': exchanges,
              'scope': 'Exact rational Hardy-kernel checks at surrogate nodes, '
                       'not finite polynomial interpolation rank or an evaluation '
                       'of pi. The written Hermite proof gives the uniform '
                       'actual-period upper bound, with parameters fixed before H.'}
    Path(__file__).with_name('small_escape_verification.json').write_text(
        json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'hermite_inverse_cases': len(records),
                      'cardinal_polynomials': sum(r['cardinal_polynomials'] for r in records),
                      'exchange_cases': exchanges['exchange_cases']}, indent=2))


if __name__ == '__main__':
    main()
