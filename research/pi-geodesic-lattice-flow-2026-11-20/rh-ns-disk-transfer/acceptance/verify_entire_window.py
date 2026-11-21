"""Exact fixed test, X=4096 and signed bandwidth L=40."""
from collections import Counter
from fractions import Fraction
from math import gcd, isqrt

X = 4096
Q = 256
L = 40


def twice_grid_area(n):
    # First-quadrant upper hull from north to east, with collinear midpoints dropped.
    chain = []
    for x in range(n + 1):
        y = isqrt(n*n - x*x)
        while len(chain) >= 2:
            ax, ay = chain[-2]
            bx, by = chain[-1]
            if (bx-ax)*(y-ay) < (by-ay)*(x-ax):
                break
            chain.pop()
        chain.append((x, y))
    return 4 * sum(y*u - x*v for (x, y), (u, v) in zip(chain, chain[1:]))


def main():
    first = {q: next(r for r in range(1, q + 1) if 2*r*r % q == 0)
             for q in range(1, Q + 1)}
    hits = Counter()
    positive = Counter()
    for n in range(X, 2*X):
        area = twice_grid_area(n)
        denominator = 2*n*n
        for q in range(1, Q + 1):
            if denominator % q:
                continue
            lower = q * area // denominator
            # Q*(L/2)*X^(-4/3)=5/64<1, so these are the only possible p.
            for p in (lower, lower + 1):
                if gcd(p, q) != 1:
                    continue
                offset = (denominator // q)*p - area
                if abs(offset)**3 > L**3 * n*n:
                    continue
                hits[p, q] += 1
                positive[p, q] += offset > 0
    total = Fraction(0)
    for (p,q), count in sorted(hits.items(),key=lambda e: Fraction(*e[0])):
        weight = Fraction(count * first[q] * Q*Q, X*q*q)
        print(f"{p}/{q}: {count} hits, {positive[p,q]} positive offsets, weight {weight}")
        total += weight
    print("distinct",len(hits),"hits",sum(hits.values()),"positive",sum(positive.values()),"F",total)


if __name__ == "__main__":
    main()
