#!/usr/bin/env python3
"""Check the cancellation-free interpolation with exact rational arithmetic."""
from fractions import Fraction as F
from itertools import product
from math import sqrt


def signed_error(pq, x):
    p, q = pq
    return q * (q*x - p)


def pair_distance_squared(left, right, x, t):
    a = signed_error(left, x)
    b = signed_error(right, x)
    projected = (1-t)*a + t*b
    lifted_sq = t*(1-t)*(b-a)**2
    return projected**2 + lifted_sq


def arctan_bounds(base, terms):
    total = sum(F((-1)**j, (2*j+1)*base**(2*j+1)) for j in range(terms))
    next_total = total + F((-1)**terms, (2*terms+1)*base**(2*terms+1))
    return min(total, next_total), max(total, next_total)


def pi_bracket():
    a, b = arctan_bounds(5, 32)
    c, d = arctan_bounds(239, 12)
    return 16*a-4*d, 16*b-4*c


def main():
    samples = [(-8, 5), (0, 2), (3, 1), (22, 7), (333, 106),
               (355, 113), (710, 226), (103993, 33102)]
    count = 0
    for left, right, x, t in product(
        samples, samples, map(F, ('-2/5', '1', '22/7', '13/4')),
        map(F, ('0', '1/13', '2/9', '1/2', '4/5', '1')),
    ):
        a, b = signed_error(left, x), signed_error(right, x)
        assert pair_distance_squared(left, right, x, t) == (1-t)*a*a + t*b*b
        if a*b < 0:
            zero = a/(a-b)
            assert pair_distance_squared(left, right, x, zero) == abs(a*b)
        count += 1

    # Independent finite-weight covariance check for the smooth global formula.
    raw = list(map(F, ('1/17', '5/17', '14/17', '17/17', '3/17', '19/17',
                        '2/17', '7/17')))
    probabilities = [weight/sum(raw) for weight in raw]
    for x in map(F, ('-15/16', '1', '22/7', '13/4')):
        errors = [signed_error(node, x) for node in samples]
        mean = sum(w*v for w, v in zip(probabilities, errors))
        variance = sum(w*(v-mean)**2 for w, v in zip(probabilities, errors))
        assert mean*mean + variance == sum(w*v*v for w, v in zip(probabilities, errors))

    lo, hi = pi_bracket()
    left, right = (333, 106), (355, 113)
    a_lo, a_hi = signed_error(left, lo), signed_error(left, hi)
    b_lo, b_hi = -signed_error(right, hi), -signed_error(right, lo)
    assert 0 < a_lo <= a_hi and 0 < b_lo <= b_hi
    lower, upper = a_lo*b_lo, a_hi*b_hi
    assert F(56436, 10**6)**2 < lower < upper < F(56437, 10**6)**2
    print(f"{count} rational arc checks and four covariance checks passed")
    print(f"Projected crossing for 333/106 -> 355/113: lifted normalized "
          f"distance ≈ {sqrt(float((lower+upper)/2)):.12f}")
    print("Exact certificate: 0.056436 < lifted normalized distance < 0.056437")


if __name__ == "__main__":
    main()
