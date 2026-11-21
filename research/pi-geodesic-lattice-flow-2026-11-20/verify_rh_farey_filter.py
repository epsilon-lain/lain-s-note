"""Checks for the two fixed-phase Farey-window research tasks. Python 3."""
from cmath import exp as complex_exp
from fractions import Fraction as F
from math import cos, exp, gcd, log, log1p, pi, sin, sinh


def mobius(n):
    sign = 1
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            sign = -sign
            if n % p == 0:
                return 0
        p += 1
    return -sign if n > 1 else sign


def totient(q):
    return sum(gcd(a, q) == 1 for a in range(q))


def arithmetic_coeffs(Q, H):
    """Integer Fourier numerators. Divide by L^2 for their values."""
    L = H * Q * Q
    out = [0] * L
    for d in range(1, 2 * Q):
        weight = sum(mobius(r) for r in range((Q + d - 1) // d, (2 * Q - 1) // d + 1))
        for k in range(0, L, d):
            out[k] += (L - k) * d * weight
    return out


def cycle_traces(q, count):
    """Recover primitive q-th root power sums by removing proper suborders."""
    rows = {}
    for d in range(1, q + 1):
        if q % d:
            continue
        rows[d] = [d if k % d == 0 else 0 for k in range(count)]
        for b, row in rows.items():
            if b != d and d % b == 0:
                rows[d] = [a - c for a, c in zip(rows[d], row)]
    return rows[q]


def direct_shell(Q, H):
    L = H * Q * Q
    total = 0.0
    for q in range(Q, 2 * Q):
        for a in range(q):
            if gcd(q, a) == 1:
                angle = 1 - pi * a / q
                total += (sin(L * angle) / (L * sin(angle))) ** 2
    return total


def folded_shell(Q, H):
    L = H * Q * Q
    total = 0.0
    for d in range(1, 2 * Q):
        m, rem = divmod(L, d)
        t = rem / d
        weight = sum(mobius(r) for r in range((Q+d-1)//d, (2*Q-1)//d+1))
        numerator = (1-t)*sin(m*d)**2 + t*sin((m+1)*d)**2
        total += weight * (d / L)**2 * numerator / sin(d)**2
    return total


def atan_inverse_bounds(b, terms=16):
    partial = sum((F((-1)**k, (2*k+1) * b**(2*k+1)) for k in range(terms)), F(0))
    nxt = F((-1)**terms, (2*terms+1) * b**(2*terms+1))
    return min(partial, partial + nxt), max(partial, partial + nxt)


def certify_known_peak():
    a0, a1 = atan_inverse_bounds(5)
    b0, b1 = atan_inverse_bounds(239)
    pi0, pi1 = 16*a0 - 4*b1, 16*a1 - 4*b0
    L = 32 * 256**2
    u = L * max(abs(1-F(113,355)*v) for v in (pi0, pi1))
    assert u < F(179,1000)
    assert (1-u*u/6)**2 > F(989,1000)
    print("exact Machin interval and sine inequality: F(256,32) > 989/1000")


def verify():
    for Q, H in ((1, 4), (2, 3), (4, 5), (8, 2)):
        L = H * Q * Q
        coeffs = arithmetic_coeffs(Q, H)
        traces = [cycle_traces(q,L) for q in range(Q, 2*Q)]
        assert coeffs == [(L-k)*sum(row[k] for row in traces) for k in range(L)]
        weights = {d: sum(mobius(r) for r in range((Q+d-1)//d, (2*Q-1)//d+1))
                   for d in range(1,2*Q)}
        for e in range(1,2*Q):
            assert sum(weights[d] for d in range(e,2*Q,e)) == int(Q <= e < 2*Q)
        pole_direct = sum(1/sin(1-pi*a/q)**2 for q in range(Q,2*Q)
                          for a in range(q) if gcd(a,q) == 1)
        pole_folded = sum(
            d*d/sin(d)**2 * sum(mobius(r) for r in range((Q+d-1)//d,(2*Q-1)//d+1))
            for d in range(1,2*Q))
        assert abs(pole_direct-pole_folded) < 1e-10 * pole_direct
        nodes = sum(totient(q) for q in range(Q, 2*Q))
        t = 1 / Q**2
        z = complex_exp(2j)
        angles = [pi*a/q for q in range(Q,2*Q) for a in range(q) if gcd(a,q)==1]
        ratio_log = 2*nodes*t + 2*sum(
            log(abs((exp(-2*t)*z - complex_exp(2j*a)) / (z - complex_exp(2j*a))))
            for a in angles)
        pole_log = sum(log1p((sinh(t)/sin(1-a))**2) for a in angles)
        assert abs(ratio_log - pole_log) < 1e-10
        mean = F(coeffs[0], L**2)
        second = F(coeffs[0]**2 + 2*sum(c*c for c in coeffs[1:]), L**4)
        assert mean == F(nodes, L)
        assert second >= mean * F(2*L*L+1, 3*L*L)
        integral_square = F(L*L + 2*sum(t*t for t in range(1,L)), L**4)
        assert integral_square == F(2*L*L+1,3*L**3)
        fourier = (coeffs[0]+2*sum(c*cos(2*k) for k,c in enumerate(coeffs[1:],1))) / L**2
        assert abs(direct_shell(Q,H) - fourier) < 1e-10
        assert abs(folded_shell(Q,H) - fourier) < 1e-10
        print("exact Fourier coefficients and moment tests passed:", Q, H)
    certify_known_peak()


if __name__ == "__main__":
    verify()
