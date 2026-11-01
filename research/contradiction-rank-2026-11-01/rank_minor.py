"""Minimum-cost full-row minor over Q, with exact prefix-rank certificate.

This is the standard greedy algorithm for a represented linear matroid.
It certifies finite rational matrices; it does not estimate a pi determinant.
"""
from fractions import Fraction as F


def determinant(matrix):
    n = len(matrix)
    if any(len(row) != n for row in matrix):
        raise ValueError('determinant requires a square matrix')
    a = [list(map(F, row)) for row in matrix]
    value = F(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            value = -value
        value *= a[j][j]
        for i in range(j + 1, n):
            ratio = a[i][j] / a[j][j]
            for k in range(j + 1, n):
                a[i][k] -= ratio * a[j][k]
    return value


def cheapest_minor(matrix, costs):
    """Return rank, cheapest basis, prefix ranks and nonzero determinant.

    If rank < row count, `determinant` and `minimum_cost` are None.
    Equal-cost columns are processed in input order. The prefix ranks at cost
    boundaries, and the minimum total cost, do not depend on that tie order.
    """
    costs = list(map(F, costs))
    m, n = len(matrix), len(costs)
    if not m or any(len(row) != n for row in matrix):
        raise ValueError('require a nonempty rectangular matrix matching costs')
    columns = [[F(matrix[i][j]) for i in range(m)] for j in range(n)]
    basis, chosen, profile = {}, [], []
    ordering = sorted(range(n), key=lambda j: (costs[j], j))
    for pos, j in enumerate(ordering):
        vector = columns[j][:]
        for pivot in sorted(basis):
            multiplier = vector[pivot]
            if multiplier:
                vector = [x - multiplier * y for x, y in zip(vector, basis[pivot])]
        pivot = next((i for i, x in enumerate(vector) if x), None)
        if pivot is not None:
            divisor = vector[pivot]
            basis[pivot] = [x / divisor for x in vector]
            chosen.append(j)
        if pos + 1 == n or costs[ordering[pos + 1]] != costs[j]:
            profile.append((costs[j], len(basis)))
    if len(basis) < m:
        return dict(rank=len(basis), chosen=chosen, profile=profile,
                    minimum_cost=None, determinant=None)
    minor = [[row[j] for j in chosen] for row in matrix]
    det = determinant(minor)
    value = sum((costs[j] for j in chosen), F(0))
    certificate = F(0)
    previous_rank = 0
    for cost, rank in profile:
        certificate += cost * (rank - previous_rank)
        previous_rank = rank
    assert det and certificate == value
    return dict(rank=m, chosen=chosen, profile=profile,
                minimum_cost=value, determinant=det)
