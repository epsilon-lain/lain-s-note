"""Exact finite checks for PRIMITIVE_INTERCEPT_ATLAS.zh.md. No external inputs."""
from collections import deque
from fractions import Fraction
from math import gcd


def prime_factors(n):
    prime = 2
    while prime * prime <= n:
        if n % prime == 0:
            yield prime
            while n % prime == 0:
                n //= prime
        prime += 1
    if n > 1:
        yield n


def projective_classes(r):
    units = [t for t in range(1, r) if gcd(t, r) == 1]
    index = {}
    representatives = []
    for c in range(r):
        for d in range(r):
            pair = (c, d)
            if gcd(c, d, r) != 1 or pair in index:
                continue
            orbit = {(t*c % r, t*d % r) for t in units}
            assert len(orbit) == len(units)
            for label in orbit:
                assert label not in index
                index[label] = pair
            representatives.append(pair)
    count = Fraction(r)
    for prime in prime_factors(r):
        count *= Fraction(prime + 1, prime)
    assert len(representatives) == count
    return representatives, index


def centered(value, r):
    return (value + r // 2) % r - r // 2


def det(u, v):
    return u[0]*v[1] - u[1]*v[0]


def lin(b, u, f, v):
    return (b*u[0] + f*v[0], b*u[1] + f*v[1])


def aligned(label, u, u_A):
    r = u[1]
    w_plus = lin(1, u, 1, u_A)
    orientation = det(w_plus, u_A)
    assert abs(orientation) == 1
    b = centered(orientation * det(label, u_A), r)
    f = centered(orientation * det(w_plus, label), r)
    c, d = lin(b, w_plus, f, u_A)
    assert (c % r, d % r) == label
    return c, d


def verify_layers(u, u_prev, A, endpoints=()):
    r = u[1]
    u_A = lin(A, u, 1, u_prev)
    labels, index = projective_classes(r)
    for g in (g for g in range(1, r + 1) if r % g == 0):
        expected = {label for label in labels if gcd(label[1], r) == g}
        seen = []
        for k in range(r):
            if gcd(k, g) == 1:
                p, q = lin(g, u_A, k, u)
                seen.append(index[p % r, q % r])
        assert set(seen) == expected
        if g == 1:
            assert len(seen) == len(expected) == r

    for label in labels:
        c, d = aligned(label, u, u_A)
        assert 2 * abs(d) <= (2*A + 3)*r*r
        # Independent reduction of an integer input with positive denominator.
        P = 2*r + 1
        Q = max(1, (abs(d) // r) + 7)
        numerator, denominator = r*P + c, r*Q + d
        T = gcd(numerator, denominator)
        v, m = numerator // T, denominator // T
        assert m > 0 and gcd(v, m) == 1 and gcd(T, r) == 1
        assert index[v % r, m % r] == label
        assert (T*v-c)//r == P and (T*m-d)//r == Q
        for alpha in endpoints:
            e = r * abs(u[0] - r*alpha)
            assert 2*abs(c - d*alpha) <= e
            assert m * abs(v - m*alpha) == abs(r*(P-Q*alpha) + c-d*alpha) / Fraction(T, m)
    return len(labels)


def atan_inverse_interval(d, terms):
    partial = sum((Fraction((-1)**k, (2*k+1)*d**(2*k+1)) for k in range(terms)), Fraction())
    other = partial + Fraction((-1)**terms, (2*terms+1)*d**(2*terms+1))
    return min(partial, other), max(partial, other)


def reciprocal_pi_interval():
    a, b = atan_inverse_interval(5, 20)
    c, d = atan_inverse_interval(239, 8)
    return 1 / (16*b - 4*c), 1 / (16*a - 4*d)


def step(matrix, digit, modulus):
    a, b, c, d = matrix
    return ((digit*a+b) % modulus, a, (digit*c+d) % modulus, c)


def check_omitted_classes(r):
    labels, index = projective_classes(r)
    base = (1, 0, 0, 1)
    for digit in (0, 3, 7, 15, 1, 292):
        base = step(base, digit, r)
    queue = deque([(base, ())])
    visited = {base}
    connectors = {}
    while queue and len(connectors) < len(labels):
        matrix, word = queue.popleft()
        connectors.setdefault(index[matrix[0], matrix[2]], word)
        for digit in range(1, r+1):
            nxt = step(matrix, digit, r)
            if nxt not in visited:
                visited.add(nxt)
                queue.append((nxt, word+(digit,)))
    assert len(connectors) == len(labels)
    for target, word in connectors.items():
        matrix = base
        for digit in word:
            matrix = step(matrix, digit, r)
        checkpoint = matrix
        for j in range(1, 21):
            assert index[matrix[0], matrix[2]] == target
            matrix = step(matrix, r*(j+1), r)
            assert index[matrix[0], matrix[2]] != target
            matrix = step(matrix, r, r)
            assert matrix == checkpoint
    return len(labels)


if __name__ == "__main__":
    moduli = (2, 3, 4, 6, 8, 9, 10, 12, 15, 30, 35)
    for modulus in moduli:
        # 1/modulus and 0/1 are adjacent convergents.
        verify_layers((1, modulus), (0, 1), 17)
    lower, upper = reciprocal_pi_interval()
    u, u_prev, A = (113, 355), (106, 333), 292
    edge0 = lin(A, u, 1, u_prev)
    edge1 = lin(A+1, u, 1, u_prev)
    endpoints = sorted((Fraction(*edge0), Fraction(*edge1)))
    assert endpoints[0] < lower < upper < endpoints[1]
    assert (u[0]-lower*u[1]) * (u[0]-upper*u[1]) > 0
    assert Fraction(1, A+2) < min(u[1]*abs(u[0]-x*u[1]) for x in (lower, upper))
    assert max(u[1]*abs(u[0]-x*u[1]) for x in (lower, upper)) < Fraction(1, A)
    assert verify_layers(u, u_prev, A, (lower, upper)) == 432
    assert check_omitted_classes(6) == 12
    print("Passed 11 small moduli; 432 historical charts, 355 on first line; 12 omitted-class controls x 20 pairs.")
