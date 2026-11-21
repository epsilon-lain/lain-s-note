"""Certify the boundary parity of integer disk hulls with integer arithmetic."""

from fractions import Fraction
from math import gcd, isqrt


def lower_chain(points):
    chain = []
    for x, y in points:
        while len(chain) >= 2:
            ax, ay = chain[-2]
            bx, by = chain[-1]
            if (bx - ax) * (y - ay) > (by - ay) * (x - ax):
                break
            chain.pop()
        chain.append((x, y))
    return chain


def integral_hull(radius):
    # Each vertical column lies between its two lattice endpoints.
    endpoints = set()
    for x in range(-radius, radius + 1):
        y = isqrt(radius * radius - x * x)
        endpoints.update(((x, -y), (x, y)))
    points = sorted(endpoints)
    hull = lower_chain(points)[:-1] + lower_chain(reversed(points))[:-1]
    boundary = 0
    twice_area = 0
    for a, b in zip(hull, hull[1:] + hull[:1]):
        boundary += gcd(abs(b[0] - a[0]), abs(b[1] - a[1]))
        twice_area += a[0] * b[1] - b[0] * a[1]
    return boundary, twice_area


def main():
    radii = sorted(set(range(1, 513)) | {2**n for n in range(12)})
    values = {}
    for radius in radii:
        boundary, twice_area = integral_hull(radius)
        h = isqrt(2 * radius**2)
        assert boundary % 4 == 0 and twice_area % 4 == 0
        assert boundary // 4 % 2 == twice_area // 4 % 2 == h % 2
        values[radius] = twice_area
        if radius > 1 and (radius & (radius - 1)) == 0 and h % 2:
            assert Fraction(twice_area, 2 * radius**2).denominator == radius**2 // 2
    increments = [values[2**(n + 1)] - 4 * values[2**n] for n in range(11)]
    residuals = [increments[n + 3] - 4 * increments[n] for n in range(8)]
    assert residuals == [64, 68, 32, 52, 0, 56, 232, 320]
    for n, residual in enumerate(residuals):
        assert residual // 4 % 2 == isqrt(2 * (2**(n + 4))**2) % 2
    print("Passed all radii 1..512 and all dyadic radii 1..2048.")
    print("Three-step residuals:", residuals)


if __name__ == "__main__":
    main()
