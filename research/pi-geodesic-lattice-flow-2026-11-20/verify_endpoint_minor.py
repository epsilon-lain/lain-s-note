#!/usr/bin/env python3
"""Check the joint column/row denominator clearing on a small rational jet matrix."""
from fractions import Fraction
from itertools import product
from math import comb, lcm


def multiply(left, right, cutoff):
    result = [Fraction(0)] * (cutoff + 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right[:cutoff + 1 - i]):
            result[i + j] += a * b
    return result


def coefficient(row, column, centers, truncations):
    j, s, beta = row
    h, alpha = column
    if any(a < b for a, b in zip(alpha, beta)):
        return Fraction(0)
    series = [Fraction(comb(h, k)) if k <= h else Fraction(0) for k in range(s + 1)]
    for a, b, center, truncation in zip(alpha, beta, centers, truncations):
        series = [comb(a, b) * term for term in series]
        coordinate = [j * center] + [
            Fraction((-1)**(k+1), k) if k < truncation else Fraction(0)
            for k in range(1, s + 1)
        ]
        for _ in range(a-b):
            series = multiply(series, coordinate, s)
    return series[s]


def determinant(matrix):
    grid = [row[:] for row in matrix]
    result = Fraction(1)
    for j in range(len(grid)):
        pivot = next((i for i in range(j, len(grid)) if grid[i][j]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != j:
            grid[j], grid[pivot] = grid[pivot], grid[j]
            result = -result
        value = grid[j][j]
        result *= value
        for i in range(j + 1, len(grid)):
            if grid[i][j]:
                factor = grid[i][j] / value
                for k in range(j + 1, len(grid)):
                    grid[i][k] -= factor * grid[j][k]
    return result


def greedy_columns(columns, vectors):
    pivots = {}
    chosen = []
    for column, vector in zip(columns, vectors):
        reduced = vector[:]
        for pivot, known in sorted(pivots.items()):
            coeff = reduced[pivot]
            if coeff:
                reduced = [a-coeff*b for a,b in zip(reduced, known)]
        first = next((i for i, term in enumerate(reduced) if term), None)
        if first is not None:
            pivots[first] = [term / reduced[first] for term in reduced]
            chosen.append(column)
    return chosen


def main():
    # Real rational centers suffice to test the clearing identity, which also
    # holds over Gaussian rationals. The weights here describe a toy matrix.
    centers = (Fraction(355, 113), Fraction(104348, 33215))
    denominators = tuple(center.denominator for center in centers)
    truncations = (7, 10)
    weights = (2, 3)
    H, w0, v0, K = 18, 2, 6, 3

    def weight(column):
        h, alpha = column
        return w0*h + sum(w*a for w,a in zip(weights,alpha))

    columns = [(h, (a,b)) for h,a,b in product(range(H//w0+1),
               range(H//weights[0]+1), range(H//weights[1]+1))
               if weight((h, (a,b))) <= H]
    columns.sort(key=weight)
    rows = [(j,s,(a,b)) for j,s,a,b in product(range(K), range(4), range(3), range(2))
            if v0*s + 2*(weights[0]*a+weights[1]*b) < H]
    vectors = [[coefficient(r, c, centers, truncations) for r in rows] for c in columns]
    lower_lcms = tuple(lcm(*range(1,T)) for T in truncations)
    joint_lcms = tuple(lcm(q,L) for q,L in zip(denominators, lower_lcms))
    chosen = greedy_columns(columns, vectors)
    assert len(chosen) == len(rows)
    exact = [[coefficient(row,col,centers,truncations) for col in chosen] for row in rows]
    cleared = []
    for row, values in zip(rows, exact):
        cleared_row = []
        for (_, alpha), value in zip(chosen, values):
            for b, a, d in zip(row[2], alpha, joint_lcms):
                value *= Fraction(d)**(a-b)
            assert value.denominator == 1
            cleared_row.append(value)
        cleared.append(cleared_row)
    direct_det = determinant(exact)
    integral_det = determinant(cleared)
    assert direct_det and integral_det and integral_det.denominator == 1
    exponents = [sum(col[1][i] for col in chosen)-sum(row[2][i] for row in rows)
                 for i in range(len(weights))]
    old_multiplier = 1
    new_multiplier = 1
    for i, power in enumerate(exponents):
        old_multiplier *= denominators[i]**power
        old_multiplier *= lower_lcms[i]**(len(rows)*(H//weights[i]))
        new_multiplier *= joint_lcms[i]**power
    assert direct_det*new_multiplier == integral_det
    print(f"Rank {len(chosen)} from {len(columns)} columns.")
    print(f"Old per-entry clearing multiplier: {old_multiplier.bit_length()} bits.")
    print(f"New column-and-row clearing multiplier: {new_multiplier.bit_length()} bits.")
    print(f"Exact determinant denominator: {direct_det.denominator.bit_length()} bits.")
    print(f"Sum of chosen column weights: {sum(map(weight,chosen))}; ceiling: {H*len(rows)}.")


if __name__ == "__main__":
    main()
