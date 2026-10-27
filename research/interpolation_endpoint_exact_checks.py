"""Exact finite checks for the endpoint-interpolation research note.

These checks verify small examples, not the asymptotic theorem or pi's
approximation properties. Run with Python 3; no third-party packages.
"""
from fractions import Fraction
from math import comb


def poly_mul(left, right, max_degree):
    out = [Fraction(0)] * (max_degree + 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if i + j <= max_degree:
                out[i + j] += a * b
    return out


def determinant(matrix):
    matrix = [row[:] for row in matrix]
    value = Fraction(1)
    for col in range(len(matrix)):
        pivot = next((r for r in range(col, len(matrix)) if matrix[r][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
            value = -value
        diagonal = matrix[col][col]
        value *= diagonal
        for row in range(col + 1, len(matrix)):
            if matrix[row][col] == 0:
                continue
            scale = matrix[row][col] / diagonal
            for k in range(col + 1, len(matrix)):
                matrix[row][k] -= scale * matrix[col][k]
            matrix[row][col] = 0
    return value


def entry(h, alpha, j, order, beta, center_step):
    if alpha < beta:
        return Fraction(0)
    p = [Fraction(comb(h, k)) if k <= h else Fraction(0) for k in range(order + 1)]
    # G(t) = t - t^2/2 + t^3/3 - t^4/4.
    g = [Fraction(j) * center_step]
    g += [Fraction((-1) ** (k + 1), k) if k <= 4 else Fraction(0)
          for k in range(1, order + 1)]
    for _ in range(alpha - beta):
        p = poly_mul(p, g, order)
    return comb(alpha, beta) * p[order]


def check_explicit_minor(H, center_step):
    # m=1, K=2, theta=2/5, w0=w1=1, v0=8/5.
    groups = []
    for s in range(H + 1):
        count = sum(1 for beta in range(H + 1) if 16 * s + 25 * beta < 10 * H)
        if count:
            groups.append((s, count))
    rows = [(j, s, beta) for s, count in groups for j in range(2) for beta in range(count)]
    cols = [(s, alpha) for s, count in groups for alpha in range(2 * count)]
    assert len(rows) == len(cols)
    assert all(h + alpha <= H for h, alpha in cols)
    matrix = [[entry(h, alpha, j, s, beta, center_step) for h, alpha in cols]
              for j, s, beta in rows]
    exponent = sum(count * count for _, count in groups)
    assert abs(determinant(matrix)) == abs(center_step) ** exponent
    return len(rows), exponent


def check_two_coordinate_minor(r1, r2):
    # Rows (j,s)=(0,0),(0,1),(1,0),(1,1), columns 1,X2,X2^2,X1.
    matrix = [[Fraction(1), Fraction(0), Fraction(0), Fraction(0)],
              [Fraction(0), Fraction(1), Fraction(0), Fraction(1)],
              [Fraction(1), r2, r2 * r2, r1],
              [Fraction(0), Fraction(1), 2 * r2, Fraction(1)]]
    assert determinant(matrix) == 2 * r2 * (r2 - r1)


def check_separation(x, a1, q1, a2, q2):
    assert 0 < q1 <= q2
    assert a1 * q2 != a2 * q1
    k1 = q1 * q1 * abs(x - Fraction(a1, q1))
    k2 = q2 * q2 * abs(x - Fraction(a2, q2))
    ratio = Fraction(q2, q1)
    assert 1 <= k1 * ratio + k2 / ratio


def main():
    for H in range(5, 10):
        size, exponent = check_explicit_minor(H, Fraction(2))
        check_explicit_minor(H, Fraction(-1, 2))
        print(f'one-coordinate H={H}, size={size}, exponent={exponent}: OK')
    for r1 in range(-3, 4):
        for r2 in range(-3, 4):
            check_two_coordinate_minor(Fraction(r1), Fraction(r2))
    print('two-coordinate 4-by-4 determinant: 49 exact checks OK')
    x = Fraction(3141592653589793, 10**15)  # finite rational test input, not exact pi
    check_separation(x, 22, 7, 355, 113)
    print('two-rational separation inequality: exact arithmetic check OK')


if __name__ == '__main__':
    main()
