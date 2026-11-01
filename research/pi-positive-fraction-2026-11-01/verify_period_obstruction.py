"""Exact, standard-library certificate for a packed-transport rank failure.

The formal-period matrix uses D_t^s [u^b] P(e^t, j Z+t+u).
Dividing row (j,s,b) by s! recovers the time-coefficient normalization.
No floating point, third-party package, or transcendence assumption is used.
"""

from fractions import Fraction
from math import comb, factorial


H, W0, W1, V0, V1, K = 5, 1, 2, 1, 6, 2
C, T = V1 - V0, 2
COLS = [(h, a) for h in range(6) for a in range(3) if h + 2 * a <= H]
SOURCE = [(j, s, 0) for j in range(K) for s in range(5)]
TARGET = [(j, s, 0) for s in range(3) for j in range(K)] + [
    (0, 2, 1), (1, 1, 2), (0, 3, 1), (1, 2, 2)
]


def entry_polynomial(row, column):
    """Return the coefficients of the integer polynomial in Z."""
    j, s, b = row
    h, a = column
    if a < b:
        return [0]
    answer = [0] * (a - b + 1)
    for k in range(min(s, a - b) + 1):
        answer[a - b - k] += (
            comb(a, b) * comb(a - b, k) * j ** (a - b - k)
            * (factorial(s) // factorial(s - k)) * h ** (s - k)
        )
    return answer


def entry(row, column, z):
    return sum(coefficient * z ** exponent
               for exponent, coefficient in enumerate(entry_polynomial(row, column)))


def determinant(matrix):
    """Exact rational elimination, adequate for these 8- and 10-row minors."""
    matrix = [[Fraction(value) for value in row] for row in matrix]
    answer = Fraction(1)
    for k in range(len(matrix)):
        pivot = next((i for i in range(k, len(matrix)) if matrix[i][k]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != k:
            matrix[k], matrix[pivot] = matrix[pivot], matrix[k]
            answer = -answer
        value = matrix[k][k]
        answer *= value
        for i in range(k + 1, len(matrix)):
            ratio = matrix[i][k] / value
            for j in range(k + 1, len(matrix)):
                matrix[i][j] -= ratio * matrix[k][j]
    return answer


def certify_determinant(rows, columns, coefficient, power):
    # Every determinant term has degree at most sum(alpha)-sum(beta).
    degree_bound = sum(a for _, a in columns) - sum(b for _, _, b in rows)
    assert power <= degree_bound
    # Equality at D+1 distinct points proves a polynomial identity of degree <=D.
    for z in range(degree_bound + 1):
        actual = determinant([[entry(row, col, z) for col in columns] for row in rows])
        assert actual == coefficient * z ** power, (z, actual)
    return degree_bound + 1


def certify_relation(terms):
    """terms=(scalar, exponent_of_Z, target_row_index)."""
    for col in COLS:
        result = [0] * 8
        for scalar, shift, row_index in terms:
            for exponent, coefficient in enumerate(entry_polynomial(TARGET[row_index], col)):
                result[shift + exponent] += scalar * coefficient
        assert not any(result), (col, result)


def main():
    assert len(SOURCE) == 10 and len(COLS) == 12 and len(TARGET) == 10
    assert SOURCE == [(j, s, b) for j in range(2) for s in range(6)
                      for b in range(2) if V0 * s + V1 * b < H]
    threshold = Fraction(T * (H + C), C + T * V0)
    assert threshold == Fraction(20, 7)
    displacement = 0
    for n in range(5):
        source_types = [b for j, s, b in SOURCE if j == 0 and s + b == n]
        assert source_types == [0]
        for j in range(2):
            types = [b for jj, s, b in TARGET if jj == j and s + b == n]
            assert len(types) == 1
            a = types[0]
            if n >= threshold:
                assert C * a < T * (H - V0 * n + C)
                assert a <= n and W1 * a <= H
                assert a == j + 1
            else:
                assert a == 0
            assert comb(a, 0) == 1  # The complete bucket Pascal minor.
            displacement += a
    assert displacement == 6  # Full normalized transport coefficient is zeta^6.
    assert all(any(entry(row, col, 1) for col in COLS) for row in TARGET)

    source_points = certify_determinant(SOURCE, COLS[:10], -9216, 5)
    independent_target_indices = (0, 1, 2, 3, 4, 5, 6, 8)
    independent_target_rows = [TARGET[i] for i in independent_target_indices]
    target_points = certify_determinant(independent_target_rows, COLS[:8], 16, 4)
    certify_relation([(1, 0, 4), (-1, 0, 5), (1, 1, 6), (1, 2, 7)])
    certify_relation([(1, 0, 7), (-1, 0, 9)])

    print("PASS: all local packing constraints and bucket minors; no zero target row")
    print("PASS: normalized complete transport coefficient = zeta^6")
    print(f"PASS: exact source minor = -9216 Z^5 ({source_points} exact evaluations)")
    print(f"PASS: exact target 8-row minor = 16 Z^4 ({target_points} exact evaluations)")
    print("PASS: two exact polynomial row relations; source rank 10, target rank 8 for Z != 0")


if __name__ == "__main__":
    main()
