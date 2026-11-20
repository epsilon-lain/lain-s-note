#!/usr/bin/env python3
"""Verify rational stops and the shortest x-monotone route in the half-plane."""
from fractions import Fraction
from itertools import product
from math import acosh, asinh, gcd, hypot, isclose, log, pi


def farey_order(cap):
    a, b, c, d = 0, 1, 1, cap
    yield a, b
    while c <= cap:
        yield c, d
        jump = (cap + b) // d
        a, b, c, d = c, d, jump*c - a, jump*d - b


def stack_endpoints(a, b, cap):
    high = Fraction(a, b), Fraction(1, b*b)
    bottom_denominator = (cap // b) * b
    low = high[0], Fraction(1, bottom_denominator**2)
    return high, low


def distance(left, right):
    x, y = left
    u, v = right
    cosh_minus_one = ((x-u)**2 + (y-v)**2) / (2*y*v)
    return acosh(1 + float(cosh_minus_one))


def optimized_route(stacks):
    costs = [distance(*stacks[0])] * 2
    choices = []
    for previous, current in zip(stacks, stacks[1:]):
        trial = [
            [costs[old] + distance(previous[1-old], current[new])
             for old in (0, 1)]
            for new in (0, 1)
        ]
        parent = [min((0, 1), key=row.__getitem__) for row in trial]
        choices.append(parent)
        costs = [distance(*current) + row[p] for row, p in zip(trial, parent)]
    state = min((0, 1), key=costs.__getitem__)
    best = costs[state]
    order = [state]
    for parent in reversed(choices):
        state = parent[state]
        order.append(state)
    order.reverse()
    path = [(stack[s], stack[1-s]) for stack, s in zip(stacks, order)]
    return best, path


def check_coverage(fractions, path, cap):
    seen = set()
    for (a, b), (entry, exit) in zip(fractions, path):
        assert entry[0] == exit[0] == Fraction(a, b)
        lowest, highest = sorted((entry[1], exit[1]))
        for multiple in range(1, cap//b + 1):
            height = Fraction(1, (multiple*b)**2)
            assert lowest <= height <= highest
            seen.add((multiple*a, multiple*b))
    assert seen == {(p, q) for q in range(1, cap+1) for p in range(q+1)}


def main():
    for cap in (4, 6, 16, 64, 128):
        fractions = list(farey_order(cap))
        assert fractions == sorted(
            ((p, q) for q in range(1, cap+1) for p in range(q+1)
             if gcd(p, q) == 1), key=lambda pq: Fraction(*pq))
        assert all(b*c-a*d == 1 and b+d > cap
                   for (a, b), (c, d) in zip(fractions, fractions[1:]))
        stacks = [stack_endpoints(a, b, cap) for a, b in fractions]
        best, path = optimized_route(stacks)
        check_coverage(fractions, path, cap)
        measured = sum(distance(a, b) for a, b in path)
        measured += sum(distance(before[1], after[0])
                        for before, after in zip(path, path[1:]))
        assert isclose(measured, best, rel_tol=1e-13)
        return_route = sum(2*distance(a, b) for a, b in stacks)
        return_route += sum(distance(a[0], b[0])
                            for a, b in zip(stacks, stacks[1:]))
        assert best <= return_route and best >= (len(stacks)-1)*acosh(1.5)
        if cap <= 6:
            vertical = sum(distance(*stack) for stack in stacks)
            brute = min(
                vertical + sum(
                    distance(stacks[i][1-angles[i]], stacks[i+1][angles[i+1]])
                    for i in range(len(stacks)-1))
                for angles in product((0, 1), repeat=len(stacks))
            )
            assert isclose(best, brute, rel_tol=1e-13)
        print(f"Q={cap:3d}: {len(fractions):4d} stacks, "
              f"best monotone length {best:.6f}, return route {return_route:.6f}")
    for p, q in [(3, 1), (22, 7), (355, 113), (710, 226)]:
        x, y = Fraction(p, q), Fraction(1, q*q)
        h = hypot(float(x)-pi, float(y))
        direct = acosh(((float(x)-pi)**2 + float(y)**2 + h*h)/(2*float(y)*h))
        assert isclose(direct, asinh(q*q*abs(pi-float(x))), rel_tol=1e-8)
    print("Mark coverage, Bézout adjacency, tiny route optima, and the ray distance checked.")


if __name__ == "__main__":
    main()
