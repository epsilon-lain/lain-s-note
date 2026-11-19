"""Exchange one column to force a nonzero diagonal vanishing order.

Usage: python verify_forced.py [path/to/verify_newton.py]. Standard library.
"""
from math import comb
from pathlib import Path
from runpy import run_path
from sys import argv


def forced_entry(row, degree, modulus, imaginary_unit):
    j, s, b1, b2 = row
    if s or b1 + b2 > degree:
        return 0
    # q1^(degree-b1) q2^(degree-b2) times the (b1,b2) jet of (X1-X2)^degree.
    # The logarithm cancels exactly, so every positive-time row vanishes.
    cross = 22 * 113 - 355 * 7
    return (comb(degree, b1) * comb(degree - b1, b2) * (-1)**b2
            * 7**b2 * 113**b1
            * pow(2 * imaginary_unit * j * cross, degree - b1 - b2, modulus)) % modulus


def expansion_coefficients(matrix, target, prime):
    # Solve M x = target in F_prime, with columns in the given order.
    m = [row + [c] for row, c in zip(matrix, target)]
    for k in range(len(m)):
        i = next(i for i in range(k, len(m)) if m[i][k])
        m[k], m[i] = m[i], m[k]
        u = pow(m[k][k], -1, prime)
        m[k] = [(x * u) % prime for x in m[k]]
        for i in range(len(m)):
            if i != k and m[i][k]:
                t = m[i][k]
                m[i] = [(x - t * y) % prime for x, y in zip(m[i], m[k])]
    return [row[-1] for row in m]


def experiment(height, verifier):
    rows = verifier['row_packet'](height)
    cols = verifier['greedy_basis'](rows, verifier['columns'](height))
    ent = verifier['entry']
    p, unit, q = 65537, 256, 113
    matrix = [[ent(row, col, p, unit) for col in cols] for row in rows]
    degree = height // 3
    max_order = max(b1 + b2 for _, _, b1, b2 in rows)
    assert degree > max_order
    coords = expansion_coefficients(matrix, [forced_entry(row, degree, p, unit) for row in rows], p)
    # Among columns replaceable over F_p, keep as many high-degree Newton factors as possible.
    replace = min((i for i, c in enumerate(coords) if c), key=lambda i: (cols[i][2], i))
    remaining = cols[:replace] + cols[replace+1:]
    nominal = sum(c[2] for c in remaining) + degree - sum(row[3] for row in rows)
    guarantee = sum(max(0, sum(c[2] >= 2*k for c in remaining) - sum(r[3] >= k for r in rows))
                    for k in range(1, height // 6 + 1))
    precision = nominal + 100
    modulus = q**precision
    values = []
    for initial_root in (15, -15):
        root = verifier['lifted_i'](q, initial_root, precision)
        square = [[forced_entry(row, degree, modulus, root) if i == replace
                   else ent(row, col, modulus, root) for i, col in enumerate(cols)] for row in rows]
        values.append(verifier['determinant_valuation'](square, q, precision))
    assert min(values) >= guarantee
    out = (cols[replace], coords[replace], degree, degree-max_order, nominal, guarantee, values, max(0, nominal-min(values)))
    print(f'H={height}; replaced {out[0]} (exchange coefficient {out[1]} mod {p}) ' 
          f'by (X1-X2)^{degree}; diagonal order >= {out[3]}; '
          f'nominal={nominal}; guaranteed removal={guarantee}; valuations={values}; '
          f'actual 113 denominator exponent={out[-1]}')
    return out


if __name__ == '__main__':
    verifier_path = Path(argv[1]) if len(argv) > 1 else Path(__file__).with_name('verify_newton.py')
    verifier = run_path(str(verifier_path))
    assert experiment(14, verifier) == ((0, 0, 1), 63811, 4, 1, 81, 28, [31, 31], 50)
    assert experiment(18, verifier) == ((0, 0, 1), 28663, 6, 2, 177, 71, [74, 74], 103)
