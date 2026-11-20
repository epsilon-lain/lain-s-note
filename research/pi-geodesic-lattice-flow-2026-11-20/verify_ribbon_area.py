#!/usr/bin/env python3
"""Exact finite checks for the rational-slope ribbon on the unit two-torus.

The curve is (t, theta*t) for 0 <= t < T, thickened vertically by an
interval of full width lam/T. Tests take T >= 1 for lam <= 1 and T >= 2 for 1 < lam <= 2.
It uses only the Python standard library; run python3 verify_ribbon_area.py.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import combinations
from math import gcd


def ceil(x):
    return -(-x.numerator // x.denominator)


def primitive_lanes(theta, time, width):
    """Enumerate all reduced p/q with 0<q<T and |q theta-p| < width."""
    a, b = theta.as_integer_ratio()
    w_num, w_den = width.as_integer_ratio()
    bound = b * w_num
    modulus = b * w_den
    lanes = []
    for q in range(1, ceil(time)):
        center = q * a * w_den
        first = (center - bound) // modulus + 1
        last = (center + bound - 1) // modulus
        for p in range(first, last + 1):
            if gcd(p, q) == 1:
                lanes.append((p, q, abs(q * theta - p)))
    return lanes


def area_from_primitive_lanes(theta, time, lam):
    width = lam / time
    lanes = primitive_lanes(theta, time, width)
    if lam <= F(1, 2):
        assert len(lanes) <= 1
    return lam - sum((time - q) * (width - error) for _, q, error in lanes)


def circular_union(theta, n, width):
    """Total arc length of n intervals, by sorting centers on the circle."""
    if n == 0:
        return F(0)
    a, b = theta.as_integer_ratio()
    # A gap g/b contributes min(g/b,width) to the union.
    centers = sorted((j * a) % b for j in range(n))
    cap = b * width.numerator
    scale = width.denominator
    total = min((centers[0] + b - centers[-1]) * scale, cap)
    total += sum(min((right - left) * scale, cap)
                 for left, right in zip(centers, centers[1:]))
    return F(total, b * scale)


def area_from_fibers(theta, time, lam):
    n = time.numerator // time.denominator
    part = time - n
    width = lam / time
    area = circular_union(theta, n, width)
    if part:
        area = (1 - part) * area + part * circular_union(theta, n + 1, width)
    return area


def twice_triangle_area(z, u, v):
    return ((u[0] - z[0])*(v[1] - z[1])
            - (u[1] - z[1])*(v[0] - z[0]))


def wide_terms(theta, time, lam):
    """Primitive edge, noncollinear triple, and quadruple corrections.

    Every shape is translated so its lex-smallest (q,p) is (0,0).
    For T>=2 and lam<=2 the vertical width is <=1, so two distinct
    points with the same q never occur in a positive-area overlap.
    """
    assert time >= 2 and 1 < lam <= 2
    width = lam/time
    lanes = primitive_lanes(theta, time, width)
    points = [(q, p) for p, q, _ in lanes]
    edges = sum(((time-q)*(width-error) for _, q, error in lanes), F(0))
    higher = []
    for num_offsets in (2, 3):
        correction = F(0)
        for tail in combinations(points, num_offsets):
            shape = ((0, 0),) + tail
            signed = [q*theta-p for q, p in shape]
            height = width - (max(signed) - min(signed))
            if height <= 0:
                continue
            minors = [abs(twice_triangle_area(*triple))
                      for triple in combinations(shape, 3)]
            if not any(minors):  # Collinear chains were handled by primitive edges.
                continue
            # Area <= 2 and strict overlap imply each such triple is unimodular.
            assert all(m == 1 for m in minors)
            correction += (time - tail[-1][0])*height
        higher.append(correction)
    triples, quads = higher
    return edges, triples, quads


def wide_checks():
    trials = 0
    for theta in map(F, ['0', '-1/3', '1/2', '1/3', '2/11', '17/5']):
        for time in map(F, ['2', '13/6', '5/2', '3', '19/6', '4', '25/6',
                            '31/6', '13/2', '8', '23/2']):
            for lam in map(F, ['9/8', '5/4', '3/2', '7/4', '2']):
                edges, triples, quads = wide_terms(theta, time, lam)
                area = area_from_fibers(theta, time, lam)
                assert lam - edges + triples - quads == area, (theta, time, lam)
                if triples:
                    assert area == 1  # A noncollinear triple implies full coverage.
                trials += 1
    print(f'{trials} exact wide-ribbon checks passed, including lambda=2.')
    examples = [(F(1,3), F(25,6), F(3,2)), (F(1,3), F(4), F(3,2))]
    for i, (theta, time, lam) in enumerate(examples):
        edges, triples, quads = wide_terms(theta, time, lam)
        assert triples > 0
        assert (quads > 0) == (i == 0)
        assert area_from_fibers(theta, time, lam) == lam-edges+triples-quads == 1
        print(f'theta={theta} T={time} lambda={lam}: '
              f'{lam} - {edges} + {triples} - {quads} = 1')


def arctan_bounds(b, terms):
    total = sum(F((-1) ** j, (2*j + 1) * b ** (2*j + 1))
                for j in range(terms))
    neighbor = total + F((-1) ** terms, (2*terms + 1) * b ** (2*terms + 1))
    return min(total, neighbor), max(total, neighbor)


def theta_enclosure():
    # Machin: pi = 16 atan(1/5) - 4 atan(1/239).
    a, b = arctan_bounds(5, 51)
    c, d = arctan_bounds(239, 17)
    pi_lo, pi_hi = 16*a - 4*d, 16*b - 4*c
    lo, hi = 1/(1 + pi_hi), 1/(1 + pi_lo)
    mid = (lo + hi) / 2
    den = 10**65
    center = F(round(mid*den), den)
    error = max(abs(center - lo), abs(center - hi))
    assert error < F(1, 10**64)
    return lo, hi, center, error


def dec(x, places=15):
    with localcontext() as ctx:
        ctx.prec = max(65, places + 20)
        return format(Decimal(x.numerator)/Decimal(x.denominator), f'.{places}f')


def certified_lanes(center, slope_error, time, lam):
    width = lam / time
    # Every active lane for the true slope is contained in this outer list.
    outer = primitive_lanes(center, time, width + time*slope_error)
    # If all outer candidates are strictly active, the list is exact.
    assert all(residual + q*slope_error < width for _, q, residual in outer)
    return outer


def certified_area(center, slope_error, time, lam):
    n = ceil(time)
    # Each of the n intervals moves by at most j*slope_error. The symmetric
    # difference of an interval and its translate has length at most 2*j*error.
    # Averaging the two neighboring fiber counts preserves the larger bound.
    uncertainty = n*(n-1)*slope_error
    a = area_from_fibers(center, time, lam)
    assert a == area_from_primitive_lanes(center, time, lam)
    return a, uncertainty


def equal_rounded(low, high, places):
    assert dec(low, places) == dec(high, places)
    return dec(low, places)


def sqrt_dec_bounds(x, places=26):
    # Build rigorous rational square-root endpoints using integer square root.
    from math import isqrt
    scale = 10**places
    k = isqrt(x.numerator * scale**2 // x.denominator)
    return F(k, scale), F(k + 1, scale)


def finite_checks():
    trials = 0
    for theta in map(F, ['0', '1/2', '-4/9', '7/13', '17/5', '113/468']):
        for time in map(F, ['1', '7/6', '3/2', '2', '13/6', '17/4',
                            '8', '57/5', '199/7', '101/2']):
            for lam in map(F, ['1/9', '1/2', '3/4', '1']):
                assert area_from_fibers(theta, time, lam) == (
                    area_from_primitive_lanes(theta, time, lam)
                ), (theta, time, lam)
                trials += 1
    print(f'{trials} exact rational checks passed, including seam and equality cases.')


def main():
    finite_checks()
    wide_checks()
    lo, hi, center, slope_error = theta_enclosure()
    q, p = 468, 113
    dlo, dhi = q*lo - p, q*hi - p
    assert 0 < dlo < dhi
    lam = F(1, 2)
    print('theta=1/(1+pi), lane p/q=113/468, lambda=1/2')
    print('support T in (468, ' +
          equal_rounded(lam/dhi, lam/dlo, 12) + '...), endpoint rounded')
    tlo = sqrt_dec_bounds(lam*q/dhi)[0]
    thi = sqrt_dec_bounds(lam*q/dlo)[1]
    alo = 2*sqrt_dec_bounds(lam*q*dlo)[0] - q*dhi
    ahi = 2*sqrt_dec_bounds(lam*q*dhi)[1] - q*dlo
    print('minimum time and area (rounded):',
          equal_rounded(tlo, thi, 12), equal_rounded(alo, ahi, 12))
    print('This minimum is on the stated lane support.')
    checks = [(F('11337/2'), lam), (F(8019), F(1)), (F(16038), F(1)),
              (F(32076), F(1)), (F('280001/2'), F('3/4'))]
    areas = {}
    for time, weight in checks:
        lanes = certified_lanes(center, slope_error, time, weight)
        value, rad = certified_area(center, slope_error, time, weight)
        areas[time, weight] = value, rad
        print(f'T={time} lambda={weight} area={equal_rounded(value-rad, value+rad, 14)}'
              f' active={[(p, q) for p, q, _ in lanes]}')
        if weight == F(3, 4):
            assert [(p, q) for p, q, _ in lanes] == [(33102, 137095), (33215, 137563)]
            assert 33102*137563 - 33215*137095 in (-1, 1)
            print('Both strictly active, determinant:',
                  33102*137563 - 33215*137095)
            print('full width:', dec(weight/time, 15))
            for p_lane, q_lane, residual in lanes:
                radius = q_lane*slope_error
                print(f'absolute residual for {p_lane}/{q_lane}:',
                      equal_rounded(residual-radius, residual+radius, 15))
    # An independently certified finite reconstruction from three measured areas.
    n = 8019
    vals = [areas[F(i*n), F(1)] for i in (1, 2, 4)]
    (a, ea), (b, eb), (c, ec) = vals
    reconstructed = F(4*n, 3)*(2*a - 3*b + c)
    radius = F(4*n, 3)*(2*ea + 3*eb + ec)
    assert reconstructed - radius > F(935, 2)
    assert reconstructed + radius < F(937, 2)
    print('three-area formula recovers shell return time:', round(reconstructed))


if __name__ == '__main__':
    main()
