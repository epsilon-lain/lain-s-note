"""Bounded exact audit of zero pairing and midpoint Frobenius identification.

The infinite midpoint density conclusion uses the cited Serre theorem; this
program does not prove a prime-pool or mean-zero asymptotic upper bound.
"""
import json
from pathlib import Path


def primes_to(n):
    flags = bytearray(b'\1') * (n+1)
    flags[:2] = b'\0\0'
    for p in range(2, n+1):
        if flags[p]:
            flags[p*p:n+1:p] = b'\0' * len(range(p*p, n+1, p))
    return [p for p in range(3, n+1) if flags[p]]


def roots_mod(p):
    values = [1, 2]
    for n in range(1, p-1):
        values.append((2*(2*n+1)*values[-1] + 4*n*values[-2])
                      * pow(n+1, -1, p) % p)
    return values, {r for r, v in enumerate(values) if v % p == 0}


def trace(p):
    def chi(t):
        t %= p
        if t == 0:
            return 0
        return 1 if pow(t, (p-1)//2, p) == 1 else -1
    return -sum(chi(x*(x*x+2*x+2)) for x in range(p))


def main():
    limit = 2000
    pairing = point_counts = 0
    sole_midpoints = []
    total_active = total_zeros = total_midpoints = 0
    for p in primes_to(limit):
        values, roots = roots_mod(p)
        h = (p-1)//2
        midpoint = int(h in roots)
        assert {p-1-r for r in roots} == roots
        assert (len(roots)-midpoint) % 2 == 0
        assert int(bool(roots))*2 <= len(roots)+midpoint
        for r in range(h+1):
            # Avoid negative powers by using the first half only.
            assert values[p-1-r] % p == pow(-4, h-r, p)*values[r] % p
            pairing += 1
        ap = trace(p)
        assert ap % p == values[h] % p
        if p >= 5:
            assert midpoint == int(ap == 0)
            point_counts += 1
        total_active += bool(roots)
        total_zeros += len(roots)
        total_midpoints += midpoint
        if roots == {h}:
            sole_midpoints.append(p)
    # The non-CM test is based on actual point counts, not j-table lookup.
    assert trace(3) == 2 and trace(5) == -2
    assert 2*2-4*3 == -8 and (-2)**2-4*5 == -16
    result = {
        'status': 'passed', 'limit': limit,
        'reflection_identity_checks': pairing,
        'midpoint_independent_point_counts': point_counts,
        'total_active_primes': total_active,
        'total_digit_zeros': total_zeros,
        'total_midpoint_zeros': total_midpoints,
        'sole_midpoint_examples': sole_midpoints,
        'non_cm_certificate': [
            {'p': 3, 'points': 2, 'trace': 2, 'frobenius_discriminant': -8},
            {'p': 5, 'points': 8, 'trace': -2, 'frobenius_discriminant': -16}],
        'scope': 'Exact finite pairing/trace checks; no asymptotic first-moment '
                 'or active-prime capacity bound inferred from samples.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
