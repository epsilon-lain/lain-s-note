"""Finite rank certificate for simultaneous arcsine moment and residue cancellation.

Run with Python 3. The selected 161 columns form a square integer matrix
whose determinant is a unit modulo 1000003. The rational weights are given
by this submatrix's inverse applied to (0,...,0,2).
"""
from itertools import product
from math import gcd, isqrt, prod

GAUSSIAN_FACTORS = [(2, 1), (3, 2), (4, 1), (5, 2), (6, 1),
                    (5, 4), (7, 2), (6, 5), (8, 3), (8, 5)]
MODULUS = 1_000_003


def triangle_data():
    primes = [a*a + b*b for a, b in GAUSSIAN_FACTORS]
    hypotenuse = prod(primes)
    points = set()
    for signs in product((-1, 1), repeat=len(primes) - 1):
        x, y = 1, 0
        for (a, b), sign in zip(GAUSSIAN_FACTORS, (1,) + signs):
            c, d = a*a - b*b, 2*a*b*sign
            x, y = x*c - y*d, x*d + y*c
        x, y = sorted((abs(x), abs(y)))
        assert x*x + y*y == hypotenuse**2 and gcd(x, y) == 1
        if 10*y <= 9*hypotenuse:
            points.add((x, y))
    return primes, hypotenuse, sorted(points)


def equations(primes, points):
    rows = []
    for p in primes:
        for r in range(1, (p+1)//2):
            row = [sum(int(a % p == r) - int(a % p == p-r) for a in (x,y))
                   for x, y in points]
            if any(row):
                rows.append(row)
    assert len(rows) == 150
    for exponent in range(1, 20, 2):
        rows.append([x**exponent + y**exponent for x, y in points])
    rows.append([1] * len(points))
    return rows


def pivot_columns(rows, prime):
    """Eliminate with rows in equation order, returning pivot columns in that order."""
    basis = {}
    pivots = []
    for equation in rows:
        row = [u % prime for u in equation]
        for column in sorted(basis):
            factor = row[column]
            if factor:
                row = [(u-factor*v) % prime for u,v in zip(row,basis[column])]
        pivot = next((j for j,u in enumerate(row) if u), None)
        if pivot is None:
            raise ValueError("Dependent equations modulo the verification prime")
        inverse = pow(row[pivot], -1, prime)
        basis[pivot] = [(u*inverse) % prime for u in row]
        pivots.append(pivot)
    return pivots


if __name__ == "__main__":
    primes, hypotenuse, points = triangle_data()
    assert all(all(p % d for d in range(2,isqrt(p)+1)) for p in primes+[MODULUS])
    assert hypotenuse == 1021090952484265 and len(points) == 221
    rows = equations(primes, points)
    pivots = pivot_columns(rows, MODULUS)
    assert len(pivots) == 161 and len(set(pivots)) == 161
    print(f"H={hypotenuse}, triangles={len(points)}, rank={len(pivots)}, free weights=60")
    print("Base primes:", *primes)
