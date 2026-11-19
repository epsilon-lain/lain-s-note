from fractions import Fraction
from math import gcd


A = Fraction(2, 5)
B = Fraction(3, 7)
P, Q = 93, 86
LAYERS = (-2, -1, 0)
STARTS = (74, 37, 86)


def remainder_interval(m: int, x: Fraction, terms: int = 12):
    partial = sum(
        (Fraction((-1) ** j) * x ** (2 * (m + j) + 1) /
         (2 * (m + j) + 1) for j in range(terms)),
        Fraction(0),
    )
    next_term = x ** (2 * (m + terms) + 1) / (2 * (m + terms) + 1)
    assert terms % 2 == 0
    return partial, partial + next_term


def bareiss_det(matrix):
    a = [list(map(int, row)) for row in matrix]
    n = len(a)
    sign, previous = 1, 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k + 1, n) if a[i][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign *= -1
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                z = a[i][j] * pivot - a[i][k] * a[k][j]
                assert z % previous == 0
                a[i][j] = z // previous
        for i in range(k + 1, n):
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def main():
    assert gcd(P, Q) == 1
    points = []
    for c, x in zip(LAYERS, STARTS):
        assert (P * x + c) % Q == 0
        n = (P * x + c) // Q
        points.extend([(x, n), (x + Q, n + P)])

    # All points lie on the prescribed three integer offset layers.
    assert [Q * n - P * m for m, n in points] == [
        -2, -2, -1, -1, 0, 0
    ]
    matrix = [[1, m, n, m*m, m*n, n*n] for m, n in points]
    det = bareiss_det(matrix)
    assert det == 4

    for m, n in points:
        ra = remainder_interval(m, A)
        rb = remainder_interval(n, B)
        ratio = (ra[0] / rb[1], ra[1] / rb[0])
        print(
            f"(m,n)=({m},{n}), ratio interval=("
            f"{float(ratio[0]):.12f},{float(ratio[1]):.12f})"
        )

    print("det(A) =", det, "; det(A) mod 5 =", det % 5)
    print("PASS: exact rational tail intervals and shifted six-point determinant")


if __name__ == "__main__":
    main()
