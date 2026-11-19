"""Exact affine recognition of alternating equal-area waves. Standard library."""
from fractions import Fraction as Q
from itertools import product


def subtract(a, b):
    return a[0] - b[0], a[1] - b[1]


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def decode(a, b, c, d):
    u, v, f = subtract(c, a), subtract(d, b), subtract(c, b)
    cross = det(u, v)
    if not (cross and det(v, f) and det(u, f)):
        raise ValueError("The four vertices do not give a nondegenerate decoding.")
    k = Q(cross, det(v, f))
    ell = Q(cross, det(u, f))
    origin = (b[0] - u[0] / k, b[1] - u[1] / k)
    assert origin == (c[0] - v[0] / ell, c[1] - v[1] / ell)
    return origin, k, ell


def vertices(digits):
    points = [(1, 0), (2, 1)]
    for digit in digits:
        a, b = points[-2:]
        points.append((a[0] + digit * b[0], a[1] + digit * b[1]))
    return points


def affine(p, matrix, origin):
    return tuple(sum(x * y for x, y in zip(row, p)) + shift
                 for row, shift in zip(matrix, origin))


def check():
    count = 0
    cases = [(((1, 0), (0, 1)), (0, 0)),
             (((3, 5), (2, -1)), (7, -9)),
             (((Q(2, 3), Q(3, 4)), (-1, 2)), (Q(8, 7), -2))]
    for digits in product((1, 2, 4, 15), repeat=3):
        for matrix, origin in cases:
            pts = [affine(p, matrix, origin) for p in vertices(digits)]
            unit_area = abs(det(*matrix))
            for n in (0, 1):
                o, k, ell = decode(*pts[n:n+4])
                assert (o, k, ell) == (origin, digits[n], digits[n+1])
                a, b, c = pts[n:n+3]
                assert abs(det(subtract(b, a), subtract(c, a))) == k * unit_area
                count += 1
    shifted = vertices([3, 1, 2])
    last = shifted[-1]
    shifted[-1] = (last[0], last[1] + 2)
    assert decode(*shifted[:4])[0] != decode(*shifted[1:])[0]
    print(f'{count} exact four-vertex decodings and triangle areas verified.')
    print('A perturbed fifth vertex was rejected by the common-origin test.')


if __name__ == '__main__':
    check()
