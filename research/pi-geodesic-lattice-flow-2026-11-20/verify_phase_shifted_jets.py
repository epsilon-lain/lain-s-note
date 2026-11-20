#!/usr/bin/env python3
"""Exact checks for the phase-shifted 28-row periodic jet minors."""
from fractions import Fraction as F
from math import comb

from verify_period_jet_blocks import determinant
from verify_split_period_jets import rank


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


ROWS = (
    [(s, 0, 0) for s in range(8)]
    + [(s, 1, 0) for s in range(4)]
    + [(s, 0, 1) for s in range(2)]
)
ANCHORS = ROWS
CORE = (
    [(h, 1, 1) for h in range(2)]
    + [(h, 2, 0) for h in range(4)]
    + [(h, 1, 0) for h in range(4, 8)]
)
PHASE = [(2, 2, 1), (3, 2, 1), (4, 0, 1), (5, 0, 1)]
PURE_PHASE = [(h, 0, 0) for h in range(8, 12)]


def terms(column):
    h, a, b = column
    yield F(1), h, a, b
    if b == 1:
        # Use W^h U^a (Z - L_(7-h)(W)).
        for n in range(1, 8 - h):
            yield F((-1)**n, n), h + n, a, 0


def tangent(h, s, phase):
    return choose(h, s) * phase**(h-s) if h >= s else F(0)


def jet(row, column, center, delta, step, phase):
    s, eta, k = row
    value = F(0)
    for c, h, a, b in terms(column):
        if a < eta or b < k:
            continue
        normal = choose(a, eta) * (center * delta)**(a-eta)
        const = (center * step)**(b-k) * tangent(h, s, phase)
        log = F(0)
        if b - k == 1:
            log = sum(
                F((-1)**(l+1), l) * tangent(h, s-l, phase) / (1+phase)**l
                for l in range(1, s+1)
            )
        value += c * normal * (const + log)
    return value


def full_matrix(last_columns, delta, step, phase):
    columns = ANCHORS + CORE + last_columns
    return [
        [jet(row, col, j, delta, step, j*phase) for col in columns]
        for j in (0, 1) for row in ROWS
    ]


def leading_coefficient(columns):
    return determinant([
        [F(choose(h, s) * choose(a, eta) * choose(b, k)) for h, a, b in columns]
        for s, eta, k in ROWS
    ])


def main():
    for new_columns, expected, u_order, v_order in (
        (PHASE, F(20), 14, 8),
        (PURE_PHASE, F(160000), 10, 32),
    ):
        columns = CORE + new_columns
        assert sum(a for _, a, _ in columns) - sum(eta for _, eta, _ in ROWS) == u_order
        assert sum(h for h, _, _ in columns) - sum(s for s, _, _ in ROWS) == v_order
        assert leading_coefficient(columns) == expected
        assert max(h + 3*(a+b) for col in ANCHORS + columns for _, h, a, b in terms(col)) <= 16

        for delta, step in ((F(2, 7), F(3)), (F(-3, 5), F(-2))):
            matrix = full_matrix(new_columns, delta, step, F(0))
            assert rank(matrix) == 24
            assert all(matrix[i][j] == (i == j) for i in range(14) for j in range(14))
            assert all(matrix[i][j] == 0 for i in range(14) for j in range(14, 28))
            for phase in (F(1, 9), F(-1, 13)):
                value = determinant(full_matrix(new_columns, delta, step, phase))
                assert value != 0
                if new_columns == PURE_PHASE:
                    assert value == expected * delta**u_order * phase**v_order
            if new_columns == PHASE:
                phase = F(1, 10**6)
                ratio = determinant(full_matrix(new_columns, delta, step, phase))
                ratio /= delta**u_order * phase**v_order
                assert abs(ratio - 20 * step**4) < F(1, 100)
        print(f"rank 24 on V=0, leading coefficient {expected}, U order {u_order}, V order {v_order}")


if __name__ == "__main__":
    main()
