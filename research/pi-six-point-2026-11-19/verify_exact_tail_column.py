from fractions import Fraction
from math import lcm


POINTS = [
    (430, 465),
    (258, 279),
    (209, 226),
    (37, 40),
    (479, 518),
    (307, 332),
]
A = Fraction(2, 5)
B = Fraction(3, 7)


def arctan_truncation(m: int, x: Fraction) -> Fraction:
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
            pivot_row = next(
                (i for i in range(k + 1, n) if a[i][k] != 0), None
            )
            if pivot_row is None:
                return 0
            a[k], a[pivot_row] = a[pivot_row], a[k]
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
    rational_tail_complements = [
        arctan_truncation(m, A) + arctan_truncation(n, B)
        for m, n in POINTS
    ]

    common_denominator = 1
    for value in rational_tail_complements:
        common_denominator = lcm(common_denominator, value.denominator)

    integer_matrix = []
    for (m, n), value in zip(POINTS, rational_tail_complements):
        integer_matrix.append([
            1, m, n, m*m, m*n, common_denominator * value
        ])

    determinant = bareiss_det(integer_matrix)
    assert determinant != 0
    assert determinant % 5 == 4

    # The real tail-column determinant equals minus the rational determinant.
    # Multiplying it by common_denominator gives this exact integer up to sign.
    print("common denominator bit length:", common_denominator.bit_length())
    print("integer determinant decimal digits:", len(str(abs(determinant))))
    print("integer determinant mod 5:", determinant % 5)
    print("rational determinant magnitude (float):",
          float(Fraction(abs(determinant), common_denominator)))
    print("PASS: exact rational clearing and nonzero modular certificate")


if __name__ == "__main__":
    main()
