from fractions import Fraction
from math import gcd, lcm


A = Fraction(2, 5)
B = Fraction(3, 7)
P, Q = 332, 307
POINTS = [
    (823, 890),
    (1437, 1554),
    (430, 465),
    (1044, 1129),
    (651, 704),
    (1265, 1368),
]


def truncation(m: int, x: Fraction) -> Fraction:
    return sum(
        (Fraction((-1) ** j) * x ** (2 * j + 1) / (2 * j + 1)
         for j in range(m)),
        Fraction(0),
    )


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
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
        for i in range(k + 1, n):
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def main():
    assert gcd(P, Q) == 1 and P % 2 == 0 and Q % 2 == 1
    offsets = [Q*n - P*m for m, n in POINTS]
    assert offsets == [-6, -6, -5, -5, -4, -4]
    assert all((m - n) % 2 != 0 for m, n in POINTS)

    rationals = [
        truncation(m, A) + truncation(n, B)
        for m, n in POINTS
    ]
    common_denominator = 1
    for r in rationals:
        common_denominator = lcm(common_denominator, r.denominator)

    matrix = [
        [1, m, n, m*m, m*n, common_denominator*r]
        for (m, n), r in zip(POINTS, rationals)
    ]
    assert all(isinstance(z, int) for row in matrix for z in row)
    determinant = bareiss_det(matrix)
    assert determinant != 0
    assert determinant % 5 == 4

    print("P,Q =", P, Q)
    print("offsets =", offsets)
    print("all six pairs have opposite parity")
    print("common denominator bits =", common_denominator.bit_length())
    print("integer determinant bits =", abs(determinant).bit_length())
    print("integer determinant mod 5 =", determinant % 5)
    print("PASS: exact rational construction, clearing, Bareiss and mod-5 certificate")


if __name__ == "__main__":
    main()
