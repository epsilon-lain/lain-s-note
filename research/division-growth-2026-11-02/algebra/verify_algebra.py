#!/usr/bin/env python3
"""Exact implementation checks for NOTE.md; no floating-point π input required."""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json
import random

COUNTS = {}

def check(group, condition):
    assert condition, group
    COUNTS[group] = COUNTS.get(group, 0) + 1

def trim(p):
    p = list(map(F, p))
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p

def add(p, q):
    return trim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0)
                 for i in range(max(len(p), len(q)))])

def scale(p, c):
    return trim([F(c) * x for x in p])

def mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return trim(out)

def deriv(p):
    return trim([i * p[i] for i in range(1, len(p))] or [0])

def evaluate(p, t):
    out = F(0)
    for a in reversed(p):
        out = out * t + a
    return out

def wronskian(p, q):
    return add(mul(deriv(p), q), scale(mul(p, deriv(q)), -1))

def rank(matrix):
    m = [list(map(F, row)) for row in matrix]
    r = 0
    for c in range(len(m[0])):
        pivot = next((i for i in range(r, len(m)) if m[i][c]), None)
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        pivot_value = m[r][c]
        m[r] = [x / pivot_value for x in m[r]]
        for i in range(len(m)):
            if i != r:
                factor = m[i][c]
                m[i] = [x - factor * y for x, y in zip(m[i], m[r])]
        r += 1
    return r

def beta(t):
    r = t - (t.numerator // t.denominator)
    return r * (1 - r)

def gains(t):
    bt = beta(t)
    u = 1 + bt * t**2
    return 3 * u + bt * t**4, 2 * u

rng = random.Random(20261102)
for _ in range(600):
    a = F(rng.choice([i for i in range(-7, 8) if i]))
    b = F(rng.choice([i for i in range(-7, 8) if i]))
    s, u, v, h = [F(rng.randrange(-9, 10), rng.randrange(1, 11)) for _ in range(4)]
    aa, bb = a + s*u, b + s*v
    if bb:
        check('transverse_quotient', aa / bb - a / b == s*(b*u-a*v)/(b*bb))
        check('normalized_gauge', b*aa/bb == a+s*(b*u-a*v)/bb)
        w = h
        target = a/b+s*w
        compatible_u = target*v+b*w
        check('variable_target_compatibility', (a+s*compatible_u)/bb == target)
    shared_a, shared_b = a*(1+s*h), b*(1+s*h)
    if shared_b:
        check('common_mode', shared_a/shared_b == a/b)
    eta = F(3, 4)
    if abs(s*v) <= eta*abs(b):
        check('denominator_bound', abs(bb) >= (1-eta)*abs(b))
        check('error_bound', abs(aa/bb-a/b) <= abs(s)*abs(b*u-a*v)/((1-eta)*b*b))

for d in range(1, 13):
    for a, b in [(3, 2), (-2, 7), (1, -1), (5, 3)]:
        matrix = []
        for i in range(d):
            row = [0] * (2*d)
            row[i], row[d+i] = b, -a
            matrix.append(row)
        check('finite_rank', rank(matrix) == d)
    # For V=polynomials of degree <=d-1 and R=t, the image is V+tV.
    variable_matrix = [[0]*(2*d) for _ in range(d+1)]
    for i in range(d):
        variable_matrix[i][i] = 1
        variable_matrix[i+1][d+i] = -1
    check('variable_polynomial_rank', rank(variable_matrix) == d+1)

for _ in range(200):
    polynomials = [trim([rng.randrange(-4, 5) for _ in range(rng.randrange(1, 7))]) for _ in range(3)]
    a, b, h = polynomials
    check('wronskian_scaling', wronskian(mul(h,a), mul(h,b)) == mul(mul(h,h), wronskian(a,b)))
    t = F(rng.randrange(-12, 13), rng.randrange(1, 12))
    av, bv, hv = [evaluate(p,t) for p in polynomials]
    if bv and hv:
        check('polynomial_factor_value', evaluate(mul(h,a),t)/evaluate(mul(h,b),t) == av/bv)

for n in range(-40, 41):
    t = F(n)
    aa, bb = gains(t)
    check('integer_clock', t*aa == 3*n and t*bb == 2*n)

for denominator in range(1, 25):
    for numerator in range(-3*denominator, 3*denominator+1):
        t = F(numerator, denominator)
        aa, bb = gains(t)
        check('global_positive_denominator', bb >= 2)
        check('rational_cancellation', aa/bb-F(3,2) == beta(t)*t**4/(2*(1+beta(t)*t**2)))

records = []
for n in range(2, 201):
    t = F(1,n)
    aa, bb = gains(t)
    m = n**4+n-1
    q = 2*n*n*m
    p = 3*n*n*m+n-1
    rr = aa/bb
    error = rr-F(3,2)
    check('fraction_formula', rr == F(p,q))
    check('coprime_exact', gcd(p,q) == 1 and rr.denominator == q)
    check('error_exact', error == F(n-1,q))
    check('quality_exact', q*q*error == q*(n-1))
    # For 0<t<=1/2, beta=t(1-t); the leading gains are order 3,
    # while quotient error is order 5. These inequalities also certify the orders.
    check('gain_order_bounds', F(1,2)*t**3 <= bb-2 <= 2*t**3)
    check('quotient_order_bounds', t**5/F(8) <= error <= t**5/F(2))
    if n in (2, 10, 100, 200):
        records.append({'N':n, 'p':p, 'q':q, 'error':str(error), 'q_squared_error':str(q*q*error)})

report = {
    'status':'PASS',
    'method':'Python standard library Fraction arithmetic and exact polynomial coefficients',
    'counts':COUNTS,
    'total_assertions':sum(COUNTS.values()),
    'representative_rational_readouts':records,
    'scope':'Finite implementation checks support the self-contained general proofs in NOTE.md; no theorem about π approximation follows from these examples.'
}
Path(__file__).with_name('checks.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
print(json.dumps({'status':report['status'], 'total_assertions':report['total_assertions'], 'counts':COUNTS}, ensure_ascii=False, indent=2))
