"""Verify integer outputs for the 221-triangle arcsine model.

Build the exact solver and pass its path:
    c++ -O3 solve_arcsin_integer_block.cpp -lgmpxx -lgmp -o /tmp/arcsin-solve
    python3 verify_arcsin_integer_outputs.py /tmp/arcsin-solve
"""
from fractions import Fraction
from math import comb, gcd, lcm
from random import Random
from subprocess import run
from sys import argv

from verify_common_hypotenuse import equations, pivot_columns, triangle_data


def valuation(n, prime, limit):
    if n == 0:
        return limit
    count = 0
    while count < limit and n % prime == 0:
        n //= prime
        count += 1
    return count


def local_echelon(matrix, rhs, prime, precision):
    """Eliminate exact equations over Z_p, tracking lost digits."""
    modulus = prime**precision
    rows = [[entry % modulus for entry in row] + [value % modulus]
            for row, value in zip(matrix, rhs, strict=True)]
    columns = list(range(len(matrix[0])))
    pivots = []
    for k in range(len(rows)):
        order, i, j = min(
            (valuation(rows[i][j], prime, precision), i, j)
            for i in range(k, len(rows)) for j in columns
        )
        assert order < precision
        if order:
            divisor = prime**order
            modulus //= divisor
            precision -= order
            for r, row in enumerate(rows):
                if r < k:
                    rows[r] = [entry % modulus for entry in row]
                else:
                    assert all(entry % divisor == 0 for entry in row)
                    rows[r] = [entry // divisor for entry in row]
        rows[k], rows[i] = rows[i], rows[k]
        inverse = pow(rows[k][j], -1, modulus)
        rows[k] = pivot = [entry * inverse % modulus for entry in rows[k]]
        for r in range(k + 1, len(rows)):
            factor = rows[r][j]
            if factor:
                rows[r] = [(a - factor * b) % modulus for a, b in zip(rows[r], pivot)]
        columns.remove(j)
        pivots.append(j)
    return precision, rows, pivots, columns


def affine_image(local, row, prime, exponent):
    """Return the least valuation of a step and the affine offset modulo p^exponent."""
    precision, rows, pivots, free = local
    assert precision >= exponent
    modulus = prime**exponent
    row = [entry % modulus for entry in row]
    offset = 0
    for pivot, equation in zip(pivots, rows):
        factor = row[pivot]
        if factor:
            offset = (offset + factor * equation[-1]) % modulus
            row = [(a - factor * b) % modulus for a, b in zip(row, equation)]
    assert all(row[j] == 0 for j in pivots)
    return min(valuation(row[j], prime, exponent) for j in free), offset


def numerator_row(hypotenuse, points, order):
    odd_lcm = lcm(*range(1, 2*order + 2, 2))
    denominator = 4**order * odd_lcm * hypotenuse**(2*order + 1)
    terms = [(2*n+1, comb(2*n, n) * 4**(order-n) * (odd_lcm//(2*n+1))
              * hypotenuse**(2*(order-n))) for n in range(10, order+1)]
    coefficients = [sum(c * (x**m + y**m) for m, c in terms) for x, y in points]
    return coefficients, denominator


def solve_on_columns(solver, matrix, rhs, columns):
    data = [str(len(matrix))]
    for row, b in zip(matrix, rhs, strict=True):
        data.append(" ".join(str(row[j]) for j in columns) + " " + str(b))
    output = run([solver], input="\n".join(data) + "\n", capture_output=True, text=True, check=True)
    x = [Fraction(0)] * len(matrix[0])
    for j, a in zip(columns, output.stdout.splitlines(), strict=True):
        x[j] = Fraction(a)
    d = lcm(*(a.denominator for a in x))
    u = [int(a*d) for a in x]
    assert all(sum(a*b for a, b in zip(row, u)) == d*c for row, c in zip(matrix, rhs))
    return d, u


def integer_solution(solver, matrix, rhs):
    first = pivot_columns(matrix, 1_000_003)
    columns = list(range(len(matrix[0])))
    Random(2848).shuffle(columns)
    second = pivot_columns([[row[j] for j in columns] for row in matrix], 1_000_003)
    d1, u1 = solve_on_columns(solver, matrix, rhs, first)
    d2, u2 = solve_on_columns(solver, matrix, rhs, [columns[j] for j in second])
    assert gcd(d1, d2) == 1
    m = pow(d1, -1, d2)
    u = [m*a + ((1 - m*d1)//d2)*b for a, b in zip(u1, u2)]
    assert all(sum(a*b for a, b in zip(row, u)) == c for row, c in zip(matrix, rhs))
    return u


def main(solver):
    primes, hypotenuse, points = triangle_data()
    matrix = equations(primes, points)
    rhs = [0] * (len(matrix) - 1) + [2]
    integer_solution(solver, matrix, rhs)
    data = {n: numerator_row(hypotenuse, points, n) for n in range(10, 31)}

    numerator_powers_at_two = {}
    for p in [2] + primes:
        local = local_echelon(matrix, rhs, p, 175 if p == 2 else 110)
        for n, (row, denominator) in data.items():
            e = valuation(denominator, p, 1000)
            precision = e + 40 if p == 2 else e
            g, t = affine_image(local, row, p, precision)
            assert t % p**min(g, e) == 0
            if p == 2 and min(g, valuation(t, p, precision)) > e:
                numerator_powers_at_two[n] = min(g, valuation(t, p, precision)) - e
    print("minimum achievable denominator = 1: orders 10..30")
    print("certified positive numerator valuations at 2:", numerator_powers_at_two)
    assert numerator_powers_at_two[11] == 22

    row, denominator = data[11]
    augmented = [equation + [0] for equation in matrix] + [row + [-denominator]]
    u = integer_solution(solver, augmented, rhs + [0])
    n = 11
    brute = sum(
        (Fraction(comb(2*k, k), 4**k*(2*k+1)*hypotenuse**(2*k+1))
         * sum(w*(x**(2*k+1) + y**(2*k+1)) for w, (x, y) in zip(u[:-1], points)))
        for k in range(n+1)
    )
    assert brute == Fraction(u[-1], 1)
    print(f"direct integer witness at N=11: numerator has {abs(u[-1]).bit_length()} bits")
    print(f"witness l1 height has {sum(map(abs, u[:-1])).bit_length()} bits")


if __name__ == "__main__":
    main(argv[1])
