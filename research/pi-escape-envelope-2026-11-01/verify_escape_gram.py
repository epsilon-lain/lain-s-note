"""Exact rational checks for escape-count Gram and spectral envelopes.

All computations use the Python standard library. Generic matrix obstructions
are explicitly not asserted to be actual-period interpolation packets.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import comb, prod
from pathlib import Path
import json


def transpose(a):
    return [list(row) for row in zip(*a)]


def mul(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def add_identity(a):
    return [[x+F(i == j) for j, x in enumerate(row)] for i, row in enumerate(a)]


def det(a):
    if not a:
        return F(1)
    a = [[F(x) for x in row] for row in a]
    answer = F(1)
    for k in range(len(a)):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return F(0)
        if pivot != k:
            a[pivot], a[k] = a[k], a[pivot]
            answer = -answer
        value = a[k][k]
        answer *= value
        for i in range(k+1, len(a)):
            ratio = a[i][k]/value
            for j in range(k+1, len(a)):
                a[i][j] -= ratio*a[k][j]
    return answer


def inverse(a):
    n = len(a)
    a = [[F(x) for x in row]+identity(n)[i] for i, row in enumerate(a)]
    for k in range(n):
        pivot = next(i for i in range(k, n) if a[i][k])
        a[k], a[pivot] = a[pivot], a[k]
        value = a[k][k]
        a[k] = [x/value for x in a[k]]
        for i in range(n):
            if i != k:
                ratio = a[i][k]
                a[i] = [x-ratio*y for x, y in zip(a[i], a[k])]
    return [row[n:] for row in a]


def submatrix(a, rows, cols):
    return [[a[i][j] for j in cols] for i in rows]


def elementary_gram(g):
    n = len(g)
    return [sum(det(submatrix(g, s, s)) for s in combinations(range(n), k))
            for k in range(n+1)]


def determinant_marker_coefficients(a):
    """Coefficients of det(I+t A), by direct signed permutation products."""
    n = len(a)
    out = [F(0)]*(n+1)
    for p in permutations(range(n)):
        sign = (-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1, n))
        poly = [F(1)]
        for i, j in enumerate(p):
            next_poly = [F(0)]*(len(poly)+1)
            for k, coefficient in enumerate(poly):
                next_poly[k] += coefficient*F(i == j)
                next_poly[k+1] += coefficient*a[i][j]
            poly = next_poly
        out = [x+sign*y for x, y in zip(out, poly)]
    return out


def check_packet(b, c, r):
    m, n = len(b), len(c)
    x = mul(c, inverse(b))
    db = det(b)
    ex = elementary_gram(mul(transpose(x), x))
    er = elementary_gram(mul(r, transpose(r)))
    periods = [F(0)]*(m+1)
    for chosen in combinations(range(m+n), m):
        k = sum(j >= m for j in chosen)
        minor = det([row for j, row in enumerate(b+c) if j in chosen])
        periods[k] += minor*minor
    assert periods == [db*db*v for v in ex]
    marker = determinant_marker_coefficients(mul(r, x))
    products = []
    minor_checks = 0
    for k in range(m+1):
        coefficient = F(0)
        rsum = xsum = F(0)
        for inside in combinations(range(m), k):
            for outside in combinations(range(n), k):
                rm = det(submatrix(r, inside, outside))
                xm = det(submatrix(x, outside, inside))
                coefficient += rm*xm
                rsum += rm*rm
                xsum += xm*xm
                minor_checks += 1
        assert coefficient == marker[k]
        assert rsum == er[k] and xsum == ex[k]
        assert coefficient*coefficient <= rsum*xsum
        products.append(str(coefficient))
    full = det(add_identity(mul(r, x)))
    for u in (F(1, 7), F(1), F(5), F(13, 2)):
        q = sum(v*u**(-k) for k, v in enumerate(er))
        p = sum(v*u**k for k, v in enumerate(ex))
        assert full*full <= q*p
        weighted = [[sum(b[t][i]*b[t][j] for t in range(m))
                     + u*sum(c[t][i]*c[t][j] for t in range(n))
                     for j in range(m)] for i in range(m)]
        assert det(weighted) == sum(v*u**k for k, v in enumerate(periods))
    # The support-aware weights are rational; Sylvester's identity removes roots.
    e = [F(i+1, i+2) for i in range(m)]
    d = [F(2*i+3, i+1) for i in range(n)]
    left = [[e[i]*sum(r[i][k]*r[j][k]/d[k] for k in range(n))
             for j in range(m)] for i in range(m)]
    right = [[sum(x[k][i]*d[k]*x[k][j] for k in range(n))/e[i]
              for j in range(m)] for i in range(m)]
    support_bound = det(add_identity(left))*det(add_identity(right))
    assert full*full <= support_bound
    return {'M': m, 'outside_rows': n, 'matching_minor_checks': minor_checks,
            'marker_coefficients': products,
            'relative_determinant': str(full),
            'support_weighted_squared_upper_bound': str(support_bound)}


def check_known_spectrum():
    # U has orthonormal columns and V is orthogonal, all over Q.
    u = [[F(3,5), 0], [F(4,5), 0], [0, 1]]
    v = [[F(5,13), F(12,13)], [-F(12,13), F(5,13)]]
    assert mul(transpose(u), u) == identity(2)
    assert mul(transpose(v), v) == identity(2)
    x = mul(mul(u, [[3, 0], [0, 2]]), v)
    records = []
    for k in (0, 1, 2):
        minors = [abs(det(submatrix(x, ii, jj)))
                  for ii in combinations(range(3), k)
                  for jj in combinations(range(2), k)]
        top = prod((3, 2)[:k])
        f = max(minors)
        count = comb(3, k)*comb(2, k)
        assert f <= top and count*f*f >= top*top
        records.append({'k': k, 'maximum_coordinate_minor': str(f),
                        'top_singular_product': top, 'minor_count': count})
    return records


def check_singular_source():
    b = [[F(1),F(0)],[F(0),F(0)]]
    c = [[F(0),F(1)],[F(1),F(1)]]
    # G=[[2,1],[1,2]], detG=3. The whitened eigenvalues are 0 and 2/3:
    # G^{-1} B^T B has the same eigenvalues as the Hermitian whitened matrix.
    g = [[sum(b[t][i]*b[t][j]+c[t][i]*c[t][j] for t in range(2))
          for j in range(2)] for i in range(2)]
    q_similar = mul(inverse(g),mul(transpose(b),b))
    assert det(g) == 3 and det(q_similar) == 0
    assert sum(q_similar[i][i] for i in range(2)) == F(2,3)
    counts = [F(0)]*3
    for chosen in combinations(range(4),2):
        k = sum(j>=2 for j in chosen)
        value = det([row for j,row in enumerate(b+c) if j in chosen])
        counts[k] += value*value
    assert counts == [0,2,1]
    for u in (F(1,3),F(1),F(5)):
        assert sum(counts[k]*u**k for k in range(3)) == 3*u*(F(2,3)+u/F(3))
    return {'source_rank':1,'full_column_rank':2,'whitened_eigenvalues':['0','2/3'],
            'escape_square_sum_coefficients':[str(x) for x in counts]}


def check_obstructions():
    records = []
    for m in (2, 3, 5):
        for h in (1, 3, 8):
            r, s = F(1, 2**h), F(2**(2*h))
            rr = [[r*F(i == j) for j in range(m)]+[F(0)]*m for i in range(m)]
            aligned = [[s*F(i == j) for j in range(m)] for i in range(m)]+[[F(0)]*m for _ in range(m)]
            disjoint = aligned[m:]+aligned[:m]
            assert mul(transpose(aligned), aligned) == mul(transpose(disjoint), disjoint)
            da = det(add_identity(mul(rr, aligned)))
            dd = det(add_identity(mul(rr, disjoint)))
            assert da == (1+r*s)**m and dd == 1
            bound = sum(F(comb(m,k))*(r*s)**k for k in range(m+1))
            assert bound == da
            # u=r/s attains the optimum of the scalar weighted bound.
            weight = r/s
            squared = (1+r*r/weight)**m*(1+weight*s*s)**m
            assert squared == da*da
            # A source-relative single-exit spike with bounded matrix entries.
            b = identity(m);b[0][0] = F(1, 2**(m*h))
            c = [[F(j == 0) for j in range(m)]]
            tiny = [[F(i == 0, 2**h)] for i in range(m)]
            actual = det([[b[i][j]+mul(tiny,c)[i][j] for j in range(m)] for i in range(m)])
            assert actual == F(1,2**(m*h))+F(1,2**h)
            assert actual/det(b) == 1+2**((m-1)*h)
            records.append({'M':m,'H':h,'aligned_relative_determinant':str(da),
                            'disjoint_relative_determinant':str(dd),
                            'single_exit_relative_determinant':str(actual/det(b))})
    return records


def main():
    packets = [
        check_packet([[2,1],[1,1]], [[1,2],[3,-1],[2,4]],
                     [[F(1,3),0,F(1,5)],[0,F(1,7),-F(1,2)]]),
        check_packet([[1,1,0],[0,1,1],[1,0,2]],
                     [[1,0,2],[2,1,0],[-1,2,1],[0,1,3]],
                     [[F(1,2),0,F(1,5),0],[0,F(1,3),0,F(1,7)],
                      [F(1,11),0,-F(1,13),F(1,17)]])]
    result = {'status':'passed','packets':packets,
              'known_spectrum_envelope_checks':check_known_spectrum(),
              'singular_source_check':check_singular_source(),
              'generic_obstructions':check_obstructions(),
              'scope':'Exact finite linear-algebra identities and generic countermodels. '
                      'No assertion that these countermodels occur in actual-period packets; '
                      'no period spectral growth or Diophantine conclusion is certified.'}
    Path(__file__).with_name('escape_gram_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','packets':len(packets),
                      'matching_minor_checks':sum(p['matching_minor_checks'] for p in packets),
                      'known_spectrum_checks':len(result['known_spectrum_envelope_checks']),
                      'singular_source_checks':1,
                      'generic_obstruction_cases':len(result['generic_obstructions'])},indent=2))


if __name__ == '__main__':
    main()
