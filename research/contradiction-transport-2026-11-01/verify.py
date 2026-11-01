"""Reproduce finite exact checks and algebraic continuum identities."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import factorial, log, sqrt
from random import Random
import json
from transport import dot, energy, exhaustive, identity_terms, leq, solve


def main():
    rng = Random(20261101)
    counts = Counter()
    for dimension in (1, 2):
        destinations = list(product(range(3), repeat=dimension))
        for _ in range(40):
            supply = Counter(rng.choice(destinations) for _ in range(rng.randint(1, 4)))
            c = F(rng.randint(1, 8), rng.randint(1, 8))
            weights = tuple(F(rng.randint(1, 4), rng.randint(1, 4)) for _ in range(dimension))
            best, occupancy, flows = solve(supply, destinations, c, weights)
            assert best == exhaustive(supply, destinations, c, weights)
            identity_terms(supply, occupancy, flows, c, weights)
            counts['independent_exhaustive_optima'] += 1
            counts['quadratic_identities'] += 1

    examples = []
    for dimension in (1, 2):
        for B in (1, 3, 6, 12):
            L = 3
            destinations = list(product(range(6), repeat=dimension))
            supply = {a: B * (L - sum(a)) for a in destinations if sum(a) < L}
            c, weights = F(1, 6), (F(1),) * dimension
            best, occupancy, flows = solve(supply, destinations, c, weights)
            before = energy(supply, c, weights)
            gap, squared, drift = identity_terms(supply, occupancy, flows, c, weights)
            counts['quadratic_identities'] += 1
            monotone = all(c * supply.get(b, 0) + dot(weights, b)
                           <= c * supply.get(a, 0) + dot(weights, a)
                           for b in destinations for a in destinations if leq(b, a))
            if B <= 6:
                assert monotone and best == before and occupancy == supply
                counts['no_move_certificates'] += 1
            else:
                assert not monotone and best < before
                counts['strictly_better_upward_transport'] += 1
            examples.append(dict(dimension=dimension, B=B, mass=sum(supply.values()),
                                 before=str(before), optimum=str(best),
                                 gap=str(gap), squared=str(squared), drift=str(drift)))

    # Analytic radial integrals; all constants and dilation checks here are rational.
    # We choose parameters making the free-branch dilation gamma exactly rational.
    def mass(A, R, m):
        return A * R ** m / factorial(m + 1)

    def moment(A, R, m):
        return A * m * R ** (m + 1) / factorial(m + 2)

    def square(A, R, m):
        return 2 * A * A * R ** m / factorial(m + 2)

    continuum = []
    for m, theta, k in product(range(1, 6), (F(1, 3), F(2, 3)), (F(1), F(3, 2))):
        c = F(1, 6)
        for gamma in (F(1, 2), F(1), F(3, 2), F(2)):
            alpha = k * theta * gamma ** (m + 1)
            A = alpha / c
            M = mass(A, theta, m)
            T = k * theta * gamma
            if gamma <= 1:
                value = c * square(A, theta, m) / M
                expected = 2 * alpha / (m + 2)
                free_lower = (2 * (m + 1) * T - 2 * k * m * theta) / (m + 2)
                assert free_lower <= value
                if gamma < 1:
                    assert free_lower < value
            else:
                R = theta * gamma
                target_A = T / c
                assert target_A == A / gamma ** m
                assert mass(target_A, R, m) == M
                value = (c * square(target_A, R, m)
                         + 2 * k * (moment(target_A, R, m) - moment(A, theta, m))) / M
                expected = (2 * (m + 1) * T - 2 * k * m * theta) / (m + 2)
            assert value == expected
            counts['rational_continuum_formula_checks'] += 1
            if m == 2 and theta == F(1, 3) and k == 1:
                continuum.append(dict(m=m, theta=str(theta), k=str(k), gamma=str(gamma),
                                      alpha=str(alpha), exact_value=str(value)))

    # Numerical illustrations only: the parameter-uniform obstruction is proved
    # in README, using inequalities rather than a parameter grid.
    c = log(2) / 4
    largest = -float('inf')
    for m, theta, rho, alpha, K in product(
            (1, 2, 4, 8), (0.05, 0.3, 0.7, 0.95), (0.1, 0.5, 0.9),
            (0.001, 0.1, 1.0, 10.0), (1, 2, 8)):
        if alpha <= theta:
            J = 2 * alpha / (m + 2)
        else:
            T = (alpha * theta ** m) ** (1 / (m + 1))
            J = (2 * (m + 1) * T - 2 * m * theta) / (m + 2)
        b = m / (m + 2)
        E = alpha * (100 * K * theta ** m / (c * rho) + 4 / K)
        gap = J + b * theta - b - E  # replacing U by b increases the gap
        assert gap < 0
        largest = max(largest, gap)
        counts['illustrative_parameter_checks_not_a_proof'] += 1

    print(json.dumps(dict(status='passed', counts=counts, finite_examples=examples,
                          rational_continuum_examples=continuum,
                          largest_illustrative_gap=largest,
                          scope='Finite convex transport and standalone continuum theorem; '
                                'not a pi approximation proof, not a Lean check.'), indent=2))


if __name__ == '__main__':
    main()
