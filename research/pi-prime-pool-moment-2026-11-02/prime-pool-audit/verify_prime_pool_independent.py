"""Independent exact audit of coefficient, pairing and midpoint reductions.

Standard library only; no inference about the unproved cross-prime density.
"""
from math import comb
from pathlib import Path
import json


def primes_to(n):
    return [p for p in range(3, n+1, 2)
            if all(p % d for d in range(2, int(p**0.5)+1))]


def chi(a, p):
    r = pow(a % p, (p-1)//2, p)
    return -1 if r == p-1 else r


def recurrence_row(p):
    u = [1, 2]
    for n in range(1, p-1):
        u.append((2*(2*n+1)*u[-1]+4*n*u[-2]) * pow(n+1, -1, p) % p)
    return [a % p for a in u]


def power_row(p):
    a = [1]
    for _ in range((p-1)//2):
        b = [0] * (len(a)+2)
        for j, v in enumerate(a):
            b[j] += v
            b[j+1] -= 4*v
            b[j+2] -= 4*v
        a = [v % p for v in b]
    return a


def trace(p):
    points = 1 + sum(1 + chi(x*(x*x+2*x+2), p) for x in range(p))
    return p+1-points, points


def squarefree_part(n):
    sign = -1 if n < 0 else 1
    n = abs(n)
    d = 2
    while d*d <= n:
        while n % (d*d) == 0:
            n //= d*d
        d += 1
    return sign*n


def check():
    coefficient_checks = reflection_checks = moment_checks = 0
    root_pair_checks = midpoint_checks = 0
    midpoint_primes = []
    for p in primes_to(503):
        u = recurrence_row(p)
        a = power_row(p)
        assert len(u) == len(a) == p
        assert u == a
        coefficient_checks += p
        h = (p-1)//2
        for r in range(h+1):
            assert u[p-1-r] == pow(-4, h-r, p)*u[r] % p
            reflection_checks += 1
        roots = [r for r in range(p) if u[r] == 0]
        k = sum(u[r] == 0 for r in range(h))
        eps = int(u[h] == 0)
        assert len(roots) == 2*k+eps
        assert 2*int(bool(roots)) <= len(roots)+eps
        root_pair_checks += 1
        ap, count = trace(p)
        assert ap % p == u[h]
        if p >= 5:
            assert (ap == 0) == bool(eps)
        if eps:
            midpoint_primes.append(p)
        midpoint_checks += 1
        if p <= 127:
            for r in range(1, p):
                moment = -sum(pow(x, p-1-r, p)*chi(1-4*x-4*x*x, p)
                              for x in range(p))
                assert moment % p == u[r]
                moment_checks += 1
    reductions = []
    for p, target_ap, target_points, field_d in [(3, 2, 2, -2), (5, -2, 8, -1)]:
        ap, points = trace(p)
        assert (ap, points) == (target_ap, target_points)
        assert (-256) % p != 0 and ap % p != 0
        discriminant = ap*ap-4*p
        assert squarefree_part(discriminant) == field_d
        reductions.append({'p':p, 'points':points, 'trace':ap,
                           'frobenius_discriminant':discriminant,
                           'squarefree_field_parameter':field_d,
                           'ordinary':True})
    assert reductions[0]['squarefree_field_parameter'] != reductions[1]['squarefree_field_parameter']
    return {'status':'passed', 'maximum_prime':503,
            'coefficient_identity_checks':coefficient_checks,
            'reflection_identity_checks':reflection_checks,
            'root_pair_capacity_checks':root_pair_checks,
            'midpoint_point_count_checks':midpoint_checks,
            'character_moment_checks':moment_checks,
            'character_moment_maximum_prime':127,
            'midpoint_primes':midpoint_primes,
            'non_CM_ordinary_reduction_certificates':reductions,
            'scope':'Exact finite algebra audit. Non-CM and midpoint zero density use the written proof and cited theorems; no proof of the global mean root-capacity bound.'}


if __name__ == '__main__':
    result = check()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
