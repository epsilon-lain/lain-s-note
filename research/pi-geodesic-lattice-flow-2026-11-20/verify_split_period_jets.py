#!/usr/bin/env python3
"""Exact rank and minor-order checks for periodic and split jet blocks."""
from fractions import Fraction
from itertools import combinations
from math import comb

from verify_period_jet_blocks import (
    determinant,
    pencil_value,
    verify,
    verify_two_coordinates,
)


def split_matrix(centers, multiplicity, delta):
    return [
        [
            comb(degree, order) * (center * delta) ** (degree - order)
            if degree >= order else Fraction(0)
            for degree in range(centers * multiplicity)
        ]
        for center in range(centers)
        for order in range(multiplicity)
    ]


def rank(matrix):
    work = [row[:] for row in matrix]
    next_row = 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(next_row, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[next_row], work[pivot] = work[pivot], work[next_row]
        for row in range(next_row + 1, len(work)):
            factor = work[row][column] / work[next_row][column]
            if factor:
                for j in range(column, len(work[0])):
                    work[row][j] -= factor * work[next_row][j]
        next_row += 1
    return next_row


def minor_orders(centers, multiplicity):
    size = centers * multiplicity
    # The first multiplicity rows are [I, 0]. Clearing the first
    # multiplicity columns leaves the following lower-right block.
    block = split_matrix(centers, multiplicity, Fraction(1))
    minimum = [0]
    for order in range(1, (centers - 1) * multiplicity + 1):
        nonzero_orders = (
            sum(columns) - sum(row % multiplicity for row in rows)
            for rows in combinations(range(multiplicity, size), order)
            for columns in combinations(range(multiplicity, size), order)
            if determinant([[block[row][col] for col in columns] for row in rows])
        )
        minimum.append(min(nonzero_orders))
    return [right - left for left, right in zip(minimum, minimum[1:])]


def verify_split_two_coordinates():
    centers, cutoff, time_weight, light, heavy = 2, 16, 2, 2, 3
    groups = []
    for time in range((cutoff + time_weight - 1) // time_weight):
        for heavy_order in (1, 0):
            budget = cutoff - time_weight * time - 2 * centers * heavy * heavy_order
            count = (budget + 2 * centers * light - 1) // (2 * centers * light)
            if count > 0:
                groups.append((time, heavy_order, count))
    rows = [
        (time, heavy_order, center, light_order)
        for time, heavy_order, count in groups
        for center in range(centers)
        for light_order in range(count)
    ]
    columns = [
        (time, heavy_degree, light_degree)
        for time, heavy_degree, count in groups
        for light_degree in range(centers * count)
    ]
    max_weight = max(time + heavy * (g + d) for time, g, d in columns)
    assert max_weight == 12 < cutoff

    def coefficient(step, power):
        matrix = []
        for time, heavy_order, center, light_order in rows:
            row = []
            for h, g, d in columns:
                if h > time or g < heavy_order or d - light_order != power:
                    row.append(Fraction(0))
                    continue
                light_jet = comb(d, light_order) * center**power
                heavy_jet = comb(g, heavy_order) * pencil_value(
                    center * step, g - heavy_order, time - h
                )[time - h]
                row.append(light_jet * heavy_jet)
            matrix.append(row)
        return matrix

    size = len(rows)
    coefficients = [coefficient(Fraction(2), power) for power in range(4)]
    truncation_ranks = []
    for length in range(1, 5):
        truncated = [
            [
                coefficients[t - u][i][j] if t >= u else Fraction(0)
                for u in range(length)
                for j in range(size)
            ]
            for t in range(length)
            for i in range(size)
        ]
        truncation_ranks.append(rank(truncated))
    assert truncation_ranks == [16, 38, 62, 90]

    order = sum(count**2 for _, _, count in groups)
    for step in (Fraction(2), Fraction(-5, 3)):
        coefficients = [coefficient(step, power) for power in range(4)]
        assert rank(coefficients[0]) == 16
        for delta in (Fraction(1, 7), Fraction(-2, 11)):
            matrix = [
                [
                    sum(matrix[i][j] * delta**power
                        for power, matrix in enumerate(coefficients))
                    for j in range(size)
                ]
                for i in range(size)
            ]
            assert rank(matrix) == 28
            assert determinant(matrix) == delta**order
    assert order == 22
    print("Two normal coordinates, new columns: rank 16 -> 28, determinant order 22")
    print(f"Truncated jet ranks modulo delta^1 through delta^4: {truncation_ranks}")


def main():
    for centers, cutoff in ((2, 12), (3, 18), (4, 22)):
        verify(centers, cutoff)
    verify_two_coordinates()
    verify_split_two_coordinates()

    delta = Fraction(1, 7)
    for centers, multiplicity in ((2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (4, 2)):
        rank_at_zero = rank(split_matrix(centers, multiplicity, Fraction(0)))
        split = split_matrix(centers, multiplicity, delta)
        rank_after_split = rank(split)
        assert rank_at_zero == multiplicity
        assert rank_after_split == centers * multiplicity

        observed = minor_orders(centers, multiplicity)
        expected = [
            ell + (ell - 1) // (centers - 1)
            for ell in range(1, (centers - 1) * multiplicity + 1)
        ]
        assert observed == expected
        vanishing_order = comb(centers, 2) * multiplicity**2
        assert sum(observed) == vanishing_order
        if (centers, multiplicity) == (3, 2):
            assert determinant(split) == 16 * delta**12

        print(
            f"K={centers}, b={multiplicity}: rank {rank_at_zero}"
            f" -> {rank_after_split}, positive orders {observed},"
            f" determinant order {vanishing_order}"
        )

    left, right = Fraction(355, 113), Fraction(104348, 33215)
    assert left.denominator * right.denominator * (right - left) == -1


if __name__ == "__main__":
    main()
