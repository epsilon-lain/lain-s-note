"""Compute the exact maximal packing area on a shared rational grid."""
from argparse import ArgumentParser
from fractions import Fraction
from math import gcd, isqrt


def turn(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def hull_chain(points):
    chain = []
    for point in points:
        while len(chain) >= 2 and turn(chain[-2], chain[-1], point) <= 0:
            chain.pop()
        chain.append(point)
    return chain


def disk_hull(radius):
    # Each column is a segment between its top and bottom points.
    # Only those endpoints are needed to determine the convex hull.
    extent = [isqrt(radius**2 - x**2) for x in range(-radius, radius + 1)]
    lower = [(x, -height) for x, height in zip(range(-radius, radius + 1), extent)]
    upper = [(x, height) for x, height in zip(range(-radius, radius + 1), extent)]
    return hull_chain(lower)[:-1] + hull_chain(reversed(upper))[:-1]


def twice_area(hull):
    return sum(
        x1 * y2 - y1 * x2
        for (x1, y1), (x2, y2) in zip(hull, hull[1:] + hull[:1])
    )


def atan_bracket(base, terms=40):
    total = sum(
        Fraction((-1)**k, (2 * k + 1) * base**(2 * k + 1))
        for k in range(terms)
    )
    next_sum = total + Fraction((-1)**terms, (2 * terms + 1) * base**(2 * terms + 1))
    return sorted((total, next_sum))


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("grids", nargs="*", type=int, default=[100, 101, 106, 212, 424, 1696, 3392, 6784])
    args = parser.parse_args()
    if any(grid <= 0 for grid in args.grids):
        parser.error("grid size must be positive")
    low5, high5 = atan_bracket(5)
    low239, high239 = atan_bracket(239)
    pi_lower = 16 * low5 - 4 * high239
    benchmark = Fraction(333, 106)
    for grid in args.grids:
        hull = disk_hull(grid)
        slots = twice_area(hull)
        common_denominator = 2 * grid**2
        area = Fraction(slots, common_denominator)
        assert all(x * x + y * y <= grid**2 for x, y in hull)
        assert all(turn(hull[i - 2], hull[i - 1], hull[i]) > 0 for i in range(len(hull)))
        assert len(hull)**3 <= (80 * grid)**2
        # Twice the integer area also follows from the lattice point count
        # and Pick's theorem, independently of the shoelace products.
        points = sum(2 * isqrt(grid**2 - x**2) + 1 for x in range(-grid, grid + 1))
        boundary = sum(
            gcd(abs(x2 - x1), abs(y2 - y1))
            for (x1, y1), (x2, y2) in zip(hull, hull[1:] + hull[:1])
        )
        assert slots == 2 * points - boundary - 2
        gap_lower = pi_lower - area
        assert gap_lower > 0 and (50 * gap_lower)**3 * grid**4 > 1
        fits = common_denominator % benchmark.denominator == 0 and benchmark <= area
        print(
            f"N={grid}: vertices={len(hull)}, area={area}, "
            f"cancellation={gcd(slots, common_denominator)}, 333/106 feasible={fits}"
        )


if __name__ == "__main__":
    main()
