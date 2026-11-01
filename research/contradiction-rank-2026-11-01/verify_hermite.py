"""Exact finite checks of the multivariate monomial Hermite minor.

Uses the actual Taylor coefficient matrix, including unequal G_i(t).
Finite checks do not certify the external interpolation theorem's thresholds.
"""
from fractions import Fraction as F
from itertools import product
from math import comb
import json


def multiply(p, q, degree):
    out = [F(0)] * (degree + 1)
    for a, x in enumerate(p):
        for b, y in enumerate(q):
            if a + b <= degree:
                out[a + b] += x * y
    return out


def power(p, n, degree):
    result = [F(1)]
    for _ in range(n):
        result = multiply(result, p, degree)
    return result


def determinant(matrix):
    a = [list(map(F, row)) for row in matrix]
    out = F(1)
    for k in range(len(a)):
        pivot = next((j for j in range(k, len(a)) if a[j][k]), None)
        if pivot is None:
            return F(0)
        if pivot != k:
            a[pivot], a[k] = a[k], a[pivot]
            out *= -1
        value = a[k][k]
        out *= value
        for j in range(k + 1, len(a)):
            factor = a[j][k] / value
            for ell in range(k + 1, len(a)):
                a[j][ell] -= factor * a[k][ell]
    return out


def build(H, theta, v0, w0, weights, centers, G, K):
    groups = []
    for s in range(int(F(H) / v0) + 1):
        bounds = [range(int(theta * H / w) + 1) for w in weights[1:]]
        for gamma in product(*bounds):
            remaining = theta * (H - v0 * s) - sum(w * b for w, b in zip(weights[1:], gamma))
            J = sum(weights[0] * b < remaining for b in range(int(theta * H / weights[0]) + 1))
            if J:
                groups.append((s, gamma, J))
    groups.sort(key=lambda g: (g[0], sum(g[1]), g[1]))
    rows, columns = [], []
    for s, gamma, J in groups:
        rows.extend((j, s, (b,) + gamma) for j in range(K) for b in range(J))
        columns.extend((s, (a,) + gamma) for a in range(K * J))

    def entry(row, column):
        j, s, beta = row
        h, alpha = column
        if any(b > a for b, a in zip(beta, alpha)):
            return F(0)
        poly = [F(comb(h, ell)) for ell in range(min(h, s) + 1)]
        factor = 1
        for b, a, r, g in zip(beta, alpha, centers, G):
            factor *= comb(a, b)
            poly = multiply(poly, power([j * r] + list(g), a - b, s), s)
        return factor * (poly[s] if len(poly) > s else F(0))

    matrix = [[entry(row, col) for col in columns] for row in rows]
    predicted = F(1)
    exponent = 0
    for _, _, J in groups:
        for j in range(K):
            for k in range(j + 1, K):
                predicted *= ((k - j) * centers[0]) ** (J * J)
                exponent += J * J
    column_cost = sum(sum(w * a for w, a in zip(weights, alpha)) for _, alpha in columns)
    row_cost = sum(sum(w * b for w, b in zip(weights, beta)) for _, _, beta in rows)
    assert column_cost - row_cost == weights[0] * exponent
    actual = determinant(matrix)
    assert abs(actual) == abs(predicted)
    allowed = all(w0 * s + sum(w * a for w, a in zip(weights, alpha)) <= H
                  for s, alpha in columns)
    return dict(rows=len(rows), determinant=str(actual), predicted_abs=str(abs(predicted)),
                exponent=exponent, column_cost=str(column_cost), row_cost=str(row_cost),
                exact_cost_gap=str(column_cost - row_cost), all_columns_in_budget=allowed,
                transverse_derivative_rows=sum(any(beta[1:]) for _, _, beta in rows))


def main():
    cases = []
    inputs = [
        (2, 4, (F(1), F(1))),
        (2, 7, (F(1), F(2))),
        (3, 4, (F(1), F(1), F(1))),
    ]
    for m, H, weights in inputs:
        for diagonal in (False, True):
            centers = tuple(F(3) if diagonal else F(i + 2, i + 1) for i in range(m))
            G = tuple(tuple(F((-1) ** (ell + 1), ell) for ell in range(1, 3 + i))
                      for i in range(m))
            result = build(H, F(1, 3), F(1), F(1, 2), weights, centers, G, 2)
            assert result['transverse_derivative_rows'] > 0
            assert result['all_columns_in_budget']
            cases.append(dict(dimension=m, H=H, weights=list(map(str, weights)),
                              diagonal_centers=diagonal, centers=list(map(str, centers)), **result))
    print(json.dumps(dict(status='passed', exact_determinant_checks=len(cases), cases=cases,
                          scope='Actual finite coefficient matrices; theorem proof is in hermite.md; '
                                'no pi approximation conclusion, no external weight-threshold certification.'), indent=2))


if __name__ == '__main__':
    main()
