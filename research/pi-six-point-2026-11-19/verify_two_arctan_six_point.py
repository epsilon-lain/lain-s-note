from fractions import Fraction


A = Fraction(2, 5)
B = Fraction(3, 7)
POINTS = [
    (430, 465),
    (258, 279),
    (209, 226),
    (37, 40),
    (479, 518),
    (307, 332),
]
RATIO_BOUNDS = [
    (Fraction(100636, 100000), Fraction(100637, 100000)),
    (Fraction(101560, 100000), Fraction(101561, 100000)),
    (Fraction(99831, 100000), Fraction(99832, 100000)),
    (Fraction(100634, 100000), Fraction(100635, 100000)),
    (Fraction(102373, 100000), Fraction(102374, 100000)),
    (Fraction(103316, 100000), Fraction(103317, 100000)),
]


def remainder_interval(m: int, x: Fraction, terms: int = 12):
    """Exact alternating-series enclosure for integral_0^x t^(2m)/(1+t^2) dt."""
    s = sum(
        (
            (-1) ** j * x ** (2 * (m + j) + 1) / (2 * (m + j) + 1)
            for j in range(terms)
        ),
        Fraction(0),
    )
    nxt = x ** (2 * (m + terms) + 1) / (2 * (m + terms) + 1)
    assert terms % 2 == 0
    return s, s + nxt


def det_bareiss(matrix):
    """Exact determinant of an integer matrix, using fraction-free elimination."""
    a = [list(map(int, row)) for row in matrix]
    n = len(a)
    sign = 1
    previous_pivot = 1
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
                assert numerator % previous_pivot == 0
                a[i][j] = numerator // previous_pivot
        for i in range(k + 1, n):
            a[i][k] = 0
        previous_pivot = pivot
    return sign * a[-1][-1]


def main():
    rows = []
    for (m, n), (lower, upper) in zip(POINTS, RATIO_BOUNDS):
        assert (m - n) % 2 != 0
        a_lower, a_upper = remainder_interval(m, A)
        b_lower, b_upper = remainder_interval(n, B)
        ratio_lower = a_lower / b_upper
        ratio_upper = a_upper / b_lower
        assert lower < ratio_lower <= ratio_upper < upper
        c = 86 * n - 93 * m
        rows.append([1, m, n, m * m, m * n, n * n])
        print(
            f"(m,n)=({m},{n}), c={c}, "
            f"ratio in ({float(lower):.5f},{float(upper):.5f})"
        )

    assert [86 * n - 93 * m for m, n in POINTS] == [0, 0, -1, -1, 1, 1]
    determinant = det_bareiss(rows)
    assert determinant == 64
    assert determinant % 5 == 4
    print("determinant =", determinant)
    print("determinant mod 5 =", determinant % 5)
    print("PASS: exact rational tail enclosures and exact integer determinant")


if __name__ == "__main__":
    main()
