"""Exact finite checks for the multi-factor interpolation construction.

This is not a test of the separated-weight theorem, an asymptotic proof,
or any assertion that pi is or is not badly approximable.

Weights in these tests: (w0,v0,theta,w1,w2,F0)=(1,2,1/2,2,4,5).
The two truncated logarithms really differ. All arithmetic is rational.
"""

from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from math import comb

W0, V0, W1, W2 = 1, 2, 2, 4
THETA = F(1, 2)
T1, T2 = 5, 10


def G(n):
    return {k: F((-1) ** (k + 1), k) for k in range(1, n)}


def multiply(p, q, H):
    """Sparse product in t,u1,u2,delta, truncated by the row budget."""
    result = defaultdict(F)
    for (s, b1, b2, e), c in p.items():
        for (t, a1, a2, f), d in q.items():
            idx = (s + t, b1 + a1, b2 + a2, e + f)
            if V0 * idx[0] + 4 * idx[1] + 8 * idx[2] < H:
                result[idx] += c * d
    return {idx: c for idx, c in result.items() if c}


def special_column(H, h, a, k, j):
    # r1=1, r2=1+delta, so jr2-jr1=j*delta.
    x1 = {(0, 0, 0, 0): F(j), (0, 1, 0, 0): F(1)}
    x1.update({(t, 0, 0, 0): c for t, c in G(T1).items()})
    difference = {(0, 0, 0, 1): F(j), (0, 1, 0, 0): F(-1),
                  (0, 0, 1, 0): F(1)}
    g1, g2 = G(T1), G(T2)
    for t in range(1, T2):
        c = g2.get(t, 0) - g1.get(t, 0)
        if c:
            difference[(t, 0, 0, 0)] = c
    out = {(h, 0, 0, 0): F(1)}
    for _ in range(a):
        out = multiply(out, x1, H)
    for _ in range(k):
        out = multiply(out, difference, H)
    return out


def ceil_fraction(x):
    return -(-x.numerator // x.denominator)


def params(H, h):
    residual = F(H - V0 * h)
    U = ceil_fraction(THETA * residual / W1) - 1
    A = int(THETA * residual / (4 * W2))
    return U, A, A


def det_fraction(matrix):
    matrix = [[F(x) for x in row] for row in matrix]
    result = F(1)
    for j in range(len(matrix)):
        pivot = next((i for i in range(j, len(matrix)) if matrix[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            matrix[j], matrix[pivot] = matrix[pivot], matrix[j]
            result = -result
        d = matrix[j][j]
        result *= d
        for i in range(j + 1, len(matrix)):
            factor = matrix[i][j] / d
            for k in range(j + 1, len(matrix)):
                matrix[i][k] -= factor * matrix[j][k]
            matrix[i][j] = 0
    return result


def block_symbolic():
    # A nontrivial 2-by-2 grid of derivatives and exponents.
    H, h = 52, 10
    U, A, B = params(H, h)
    assert (U, A, B) == (7, 1, 1)
    assert h + W1 * A + W2 * (U + 1 + B) <= H
    aa = list(range(A + 1))
    bb = list(range(B + 1))
    # Directional rows are partial_z^a partial_v^b / (a! b!),
    # where z=u1 and v=u2-u1. They are combinations of admissible rows.
    assert V0 * h + (W2 / THETA) * (A + B) < H
    rows = list(product(aa, bb))
    cols = list(product(aa, bb))
    matrix = [[comb(l, i) * comb(U + 1 + b, j) if l >= i else 0
               for l, b in cols] for i, j in rows]
    # Each nonzero determinant term has the same delta exponent:
    # sum(k over columns)-sum(j over rows). Evaluate its coefficient at d=1.
    block_exponent = sum(U + 1 + b for _, b in cols) - sum(j for _, j in rows)
    coefficient = det_fraction(matrix)
    assert coefficient == 1
    assert block_exponent == (A + 1) * (B + 1) * (U + 1) == 32
    total = 0
    checked = 0
    for a, b in cols:
        k = U + 1 + b
        e = k - U
        total += e
        for center in (0, 1):
            column = special_column(H, h, a, k, center)
            if center == 1:
                assert column
            for (s, b1, b2, exponent), coefficient in column.items():
                assert s >= h and V0 * s + 4 * b1 + 8 * b2 < H
                assert exponent >= e, (a, b, center, s, b1, b2, exponent, e)
                checked += 1
    assert total == 6
    print('H=52: U=7, A=B=1, 4 columns, guaranteed full-column factor delta^6')
    print('Independent 4x4 directional-jet subdeterminant: delta^', block_exponent)
    print('Exact nonzero sparse terms checked against their column factor:', checked)


def series_multiply(p, q, cut):
    out = [F(0)] * (cut + 1)
    for i, c in enumerate(p):
        for j, d in enumerate(q):
            if i + j <= cut:
                out[i + j] += c * d
    return out


def monomial_entry(col, row, r1, r2):
    h, a, b = col
    j, s, da, db = row
    if a < da or b < db:
        return F(0)
    p = [F(comb(h, t)) if t <= h else F(0) for t in range(s + 1)]
    for r, exponent, trunc in ((r1, a - da, T1), (r2, b - db, T2)):
        g = G(trunc)
        series = [F(j * r)] + [g.get(t, F(0)) for t in range(1, s + 1)]
        for _ in range(exponent):
            p = series_multiply(p, series, s)
    return comb(a, da) * comb(b, db) * p[s]


def add_basis(basis, vector):
    vector = vector[:]
    for pivot, item in basis:
        if vector[pivot]:
            factor = vector[pivot]
            vector = [a - factor * b for a, b in zip(vector, item)]
    pivot = next((i for i, c in enumerate(vector) if c), None)
    if pivot is None:
        return False
    factor = vector[pivot]
    vector = [x / factor for x in vector]
    basis.append((pivot, vector))
    basis.sort(key=lambda pair: pair[0])
    return True


def full_row_example():
    H = 18
    rows = [(j, s, b1, b2) for j in range(2) for s in range(H)
            for b1 in range(H) for b2 in range(H)
            if 2 * s + 4 * b1 + 8 * b2 < H]
    monomials = [(h, a, b) for h in range(H + 1)
                 for a in range(H // 2 + 1) for b in range(H // 4 + 1)
                 if h + 2 * a + 4 * b <= H]
    special = []
    for h in (5, 6, 7, 8):
        U, A, B = params(H, h)
        assert A == B == 0
        k = U + 1
        assert h + W2 * k <= H
        special.append((h, 0, k))

    def vector(h, a, k):
        # Expand (Y-1)^h X1^a (X2-X1)^k in original monomial columns.
        out = [F(0)] * len(rows)
        for ypower in range(h + 1):
            cy = (-1) ** (h - ypower) * comb(h, ypower)
            for x2power in range(k + 1):
                cx = (-1) ** (k - x2power) * comb(k, x2power)
                col = (ypower, a + k - x2power, x2power)
                for pos, row in enumerate(rows):
                    out[pos] += cy * cx * monomial_entry(col, row, F(1), F(2))
        return out

    basis = []
    for spec in special:
        assert add_basis(basis, vector(*spec)), spec
    used = 0
    for col in monomials:
        if add_basis(basis, [monomial_entry(col, row, F(1), F(2)) for row in rows]):
            used += 1
        if len(basis) == len(rows):
            break
    assert len(basis) == len(rows)
    assert used == len(rows) - len(special)
    print('H=18: full row rank', len(rows), 'with', len(special),
          'engineered columns and', used, 'original monomial columns')


def asymptotic_count():
    H = 1_000_000
    tau = F(1, 20)
    assert 2 * tau < 1
    assert F(W0, V0) + THETA * tau * (1 + F(2 * W2, W1)) < 1
    hmin = ceil_fraction((1 - 2 * tau) * H / V0)
    hmax = int((1 - tau) * H / V0)
    R = 0
    for h in range(hmin, hmax + 1):
        U, A, B = params(H, h)
        assert W0 * h + W1 * A + W2 * (U + 1 + B) <= H
        R += (A + 1) * (B + 1) * (B + 2) // 2
    limit = F(15, 512) * THETA ** 3 * tau ** 4 / (V0 * W2 ** 3)
    ratio = F(R, H ** 4) / limit
    assert F(99, 100) < ratio < F(101, 100), float(ratio)
    print('H=1000000: R/H^4 divided by predicted limit =', float(ratio))


if __name__ == '__main__':
    block_symbolic()
    full_row_example()
    asymptotic_count()
