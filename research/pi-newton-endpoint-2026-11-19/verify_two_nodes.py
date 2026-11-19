#!/usr/bin/env python3
"""Exact raw jets for two centers. Standard-library Python, no decimal arithmetic."""
from fractions import Fraction as F
from math import comb, factorial, gcd, prod


def multiply_series(left, right):
    order = len(left)
    return [sum(left[k] * right[s-k] for k in range(s+1)) for s in range(order)]


def determinant(rows):
    """Bareiss after independent integer clearing of each original row."""
    from math import lcm

    matrix = []
    scale = 1
    for row in rows:
        denominator = lcm(*(entry.denominator for entry in row))
        matrix.append([int(entry * denominator) for entry in row])
        scale *= denominator
    sign = previous = 1
    size = len(matrix)
    for col in range(size-1):
        pivot_row = next((r for r in range(col, size) if matrix[r][col]), None)
        if pivot_row is None:
            return F(0)
        if pivot_row != col:
            matrix[pivot_row], matrix[col] = matrix[col], matrix[pivot_row]
            sign *= -1
        pivot = matrix[col][col]
        for i in range(col+1, size):
            for j in range(col+1, size):
                value = pivot * matrix[i][j] - matrix[i][col] * matrix[col][j]
                quotient, remainder = divmod(value, previous)
                assert remainder == 0
                matrix[i][j] = quotient
            matrix[i][col] = 0
        previous = pivot
    return F(sign * matrix[-1][-1], scale)


def raw_rows(n, t, r1, r2):
    """Coefficient of x^s in the actual raw columns, Y=1+x, X_i=j*r_i+log(1+x)."""
    height = 2*n-2
    v0 = F(height, n)
    theta = F(1, height+1)
    assert v0*(n-1) < height <= v0*n and 1/theta >= height
    largest_pure = 2*n-t-1
    assert max(largest_pure, t) <= height
    log = [F(0)] + [F((-1)**(s+1), s) for s in range(1, n)]
    rows = []
    for j in (0, 1):
        x1 = [j*r1] + log[1:]
        x2 = [j*r2] + log[1:]
        powers = [[F(1)] + [F(0)]*(n-1)]
        for _ in range(largest_pure):
            powers.append(multiply_series(powers[-1], x2))
        columns = powers + [multiply_series(x1, powers[a]) for a in range(t)]
        rows.extend(zip(*columns))
    return rows


def single_node_factor(n, t):
    """Positive Hasse determinant in the order of increasing exponents."""
    exponents = list(range(t)) + list(range(n, 2*n-t))
    vandermonde = prod(b-a for i, a in enumerate(exponents) for b in exponents[i+1:])
    coefficient, remainder = divmod(vandermonde, prod(factorial(s) for s in range(n)))
    assert remainder == 0
    return coefficient


def formula(n, t, r1, r2):
    delta = r1-r2
    if t == n:
        return delta**(n-1) * (delta+n*r2)
    return ((-1)**(t*(n-t)) * single_node_factor(n,t)
            * delta**t * r2**((n-t)**2))


def largest_part(number, base):
    """All prime powers dividing number whose underlying prime divides base."""
    supported = 1
    while (shared := gcd(number, base)) != 1:
        number //= shared
        supported *= shared
    return supported


def check_raw():
    count = 0
    for n in range(2, 11):
        for t in range(1, n+1):
            for r1, r2 in ((F(22,7),F(355,113)), (F(2,3),F(-7,5)), (F(3),F(3))):
                actual = determinant(raw_rows(n, t, r1, r2))
                assert actual == formula(n, t, r1, r2), (n,t,r1,r2)
                count += 1
    return count


def check_denominators():
    # All determinants are homogeneous in the centers in their displayed formulas.
    # Replacing 2i by 2 changes the result only by a Gaussian unit.
    count = 0
    for first, second in ((F(22,7),F(355,113)), (F(333,106),F(355,113)),
                          (F(13,5),F(25,9)), (F(0),F(1,3))):
        p1, q1 = first.as_integer_ratio()
        p2, q2 = second.as_integer_ratio()
        assert q2 % 2 == 1 and gcd(q1, q2) == 1
        cross = p1*q2 - p2*q1
        for n in range(2, 51):
            for t in sorted({1, n-1, n}):
                value = formula(n, t, 2*first, 2*second)
                if t < n:
                    degree = t+(n-t)**2
                    coefficient = single_node_factor(n,t)
                    expected = q2**degree // gcd(q2**degree, coefficient)
                    if t == n-1:
                        assert coefficient == n
                        assert (n*n-n+1)-degree == (n-1)**2
                    if t == 1:
                        assert coefficient == comb(2*n-2,n-1)
                else:
                    degree = n
                    top = p1*q2+(n-1)*p2*q1
                    expected = q2**degree // gcd(q2**degree, top)
                    assert value * q1**n * q2**n == 2**n * cross**(n-1) * top
                assert largest_part(value.denominator, q2) == expected
                count += 1
    return count


if __name__ == '__main__':
    print('Raw full-log minor checks:', check_raw())
    print('Exact denominator checks:', check_denominators())
