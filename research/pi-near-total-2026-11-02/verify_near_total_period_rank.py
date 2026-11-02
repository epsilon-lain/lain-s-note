"""Exact finite certificates for a 1D full-packet near-total escape family.

Standard library only.  The general asymptotic rank theorem is proved in the
companion note, not inferred from these finite checks.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial
from pathlib import Path
import json

PRIME = 65537


def ceilq(x):
    return -(-x.numerator // x.denominator)


def packet(hmax, k, v0, v1):
    return [(j, s, b) for j in range(k)
            for b in range(ceilq(Fraction(hmax, v1)))
            for s in range(ceilq(Fraction(hmax-v1*b, v0)))]


def selection(hmax, k, v0, v1):
    rows = packet(hmax, k, v0, v1)
    counts = Counter(s+b for j, s, b in rows if j == 0)
    packed = {n for n, count in counts.items() if 2*count-1 <= n}
    targets = []
    for j, s, b in rows:
        n = s+b
        a = b+counts[n] if n in packed else b
        targets.append((j, n-a, a))
    return rows, targets, counts, packed


@lru_cache(None)
def stirling(n, h):
    if n == 0:
        return int(h == 0)
    if h <= 0 or h > n:
        return 0
    return stirling(n-1, h-1)+h*stirling(n-1, h)


@lru_cache(None)
def entry(row, column, center, shifted):
    """Factorial-scaled Taylor row, exact integral formula reduced mod PRIME."""
    j, s, b = row
    h, a = column
    if b > a:
        return 0
    out = 0
    for ell in range(min(s, a-b)+1):
        n = s-ell
        time = factorial(h)*stirling(n, h) if shifted else pow(h, n, PRIME)
        out += (comb(a, b)*comb(a-b, ell)
                * (factorial(s)//factorial(n))*time
                * pow(j*center % PRIME, a-b-ell, PRIME))
    return out % PRIME


def pivots(matrix):
    a = [r[:] for r in matrix]
    out = []
    row = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(row, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        inverse = pow(a[row][c], -1, PRIME)
        a[row] = [(x*inverse) % PRIME for x in a[row]]
        for i in range(row+1, len(a)):
            factor = a[i][c]
            if factor:
                a[i] = [(x-factor*y) % PRIME for x, y in zip(a[i], a[row])]
        out.append(c)
        row += 1
        if row == len(a):
            break
    return out


def det_mod(matrix):
    a = [r[:] for r in matrix]
    out = 1
    for c in range(len(a)):
        pivot = next((i for i in range(c, len(a)) if a[i][c]), None)
        if pivot is None:
            return 0
        if pivot != c:
            a[c], a[pivot] = a[pivot], a[c]
            out = -out
        diagonal = a[c][c]
        out = out*diagonal % PRIME
        inverse = pow(diagonal, -1, PRIME)
        for i in range(c+1, len(a)):
            factor = a[i][c]*inverse % PRIME
            if factor:
                for k in range(c+1, len(a)):
                    a[i][k] = (a[i][k]-factor*a[c][k]) % PRIME
    return out % PRIME


def det_exact(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    out = Fraction(1)
    for c in range(len(a)):
        pivot = next((i for i in range(c, len(a)) if a[i][c]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != c:
            a[c], a[pivot] = a[pivot], a[c]
            out = -out
        value = a[c][c]
        out *= value
        for i in range(c+1, len(a)):
            factor = a[i][c]/value
            for k in range(c+1, len(a)):
                a[i][k] -= factor*a[c][k]
    return out


def finite_case(hmax, w1, p, q):
    k, w0, v0, v1 = 2, 1, 2, 8*w1
    u0, u1 = Fraction(3, 2), Fraction(3*w1)
    delta = u1-u0
    kappa = max(u0/v0, 2*delta/(v1-v0))
    assert kappa < 1 and hmax > delta/(1-kappa)
    assert hmax > w1*(k-1)/(1-max(Fraction(w0)/u0, k*w1/u1))
    rows, targets, counts, packed = selection(hmax, k, v0, v1)
    auxiliary = set(packet(hmax, k, u0, u1))
    assert set(rows) <= auxiliary and set(targets) <= auxiliary
    assert len(set(targets)) == len(rows)
    assert all(t >= 0 and w1*a <= hmax for j, t, a in targets)
    assert all(j == jj and s+b == t+a and a >= b
               for (j, s, b), (jj, t, a) in zip(rows, targets))
    outside = len(set(targets)-set(rows))
    assert outside == k*sum(counts[n] for n in packed)
    gamma = sum(a-b for (_, s, b), (_, t, a) in zip(rows, targets))
    assert gamma == k*sum(counts[n]**2 for n in packed)
    for count in set(counts[n] for n in packed):
        assert det_exact([[comb(a, b) for a in range(count, 2*count)]
                          for b in range(count)]) == 1
    columns = [(h, a) for a in range(hmax//w1+1)
               for h in range(hmax-w1*a+1)]
    center = 2*p*256*pow(q, -1, PRIME) % PRIME
    bases = []
    for shifted in (False, True):
        period = [[entry(row, col, 1, shifted) for col in columns]
                  for row in targets]
        chosen = pivots(period)
        assert len(chosen) == len(rows)
        target_det = det_mod([[line[i] for i in chosen] for line in period])
        source_det = det_mod([[entry(row, columns[i], center, shifted)
                              for i in chosen] for row in rows])
        assert target_det and source_det
        deficit = sum(s for j, s, b in rows)-sum(columns[i][0] for i in chosen)
        degree_surplus = (sum(columns[i][1] for i in chosen)
                          - sum(b for j, s, b in rows))
        assert degree_surplus >= gamma
        if shifted:
            assert deficit >= gamma
        bases.append({
            'basis': 'shifted (Y-1)^h X^a' if shifted else 'raw Y^h X^a',
            'selected_column_indices': chosen,
            'selected_columns': [columns[i] for i in chosen],
            'target_period_at_Z1_determinant_mod_prime': target_det,
            'rational_source_determinant_mod_prime': source_det,
            'source_time_sum_minus_column_h_sum': deficit,
            'transverse_degree_surplus_D': degree_surplus,
            'residual_degree_E': degree_surplus-gamma,
            'time_grade_filter_applies': shifted,
        })
    return {'H': hmax, 'K': k, 'w0': w0, 'w1': w1, 'v0': v0, 'v1': v1,
            'rational_approximation': [p, q], 'source_rows': len(rows),
            'ambient_columns': len(columns), 'outside_rows': outside,
            'outside_fraction': str(Fraction(outside, len(rows))),
            'gamma': gamma, 'normalized_transport_integer_coefficient': 1,
            'packed_buckets': sorted(packed), 'common_bases': bases}


def asymptotic_checks():
    out = []
    for w1 in [2, 5, 20, 100]:
        v0, v1, k, hmax = 2, 8*w1, 2, 10000
        counts = [min(n+1, ceilq(Fraction(hmax-v0*n, v1-v0)))
                  for n in range(ceilq(Fraction(hmax, v0)))]
        mass = k*sum(counts)
        escaped = k*sum(b for n, b in enumerate(counts) if 2*b-1 <= n)
        gamma = k*sum(b*b for n, b in enumerate(counts) if 2*b-1 <= n)
        rho = Fraction(v1*(v1-v0), (v1+v0)**2)
        cost = Fraction(2*w1*v1*(v1-v0), 3*(v1+v0)**3)
        out.append({'w1': w1, 'H': hmax, 'rows': mass,
                    'outside_fraction': float(Fraction(escaped, mass)),
                    'limiting_outside_fraction': str(rho),
                    'w1_gamma_over_MH': float(Fraction(w1*gamma, mass*hmax)),
                    'limiting_w1_gamma_over_MH': str(cost)})
    return out


def main():
    assert all(PRIME % d for d in range(2, 257))
    assert 256*256 % PRIME == PRIME-1
    cases = [finite_case(24, 2, 22, 7), finite_case(40, 2, 22, 7),
             finite_case(60, 5, 355, 113)]
    result = {'status': 'passed', 'prime': PRIME, 'image_of_i': 256,
              'finite_cases': cases, 'asymptotic_samples': asymptotic_checks(),
              'scope': 'Exact finite common raw and shifted columns. General '
                       'theorem proves full target row rank for sufficiently '
                       'large H in a one-dimensional subcritical regime; '
                       'common raw-column existence for every H is not proved.'}
    Path(__file__).with_name('near_total_verification.json').write_text(
        json.dumps(result, indent=2)+'\n')
    for case in cases:
        print(json.dumps({key: value for key, value in case.items()
                          if key not in {'common_bases', 'packed_buckets'}}))
        for base in case['common_bases']:
            print(json.dumps({key: value for key, value in base.items()
                              if key not in {'selected_column_indices', 'selected_columns'}}))
    print(json.dumps(result['asymptotic_samples'], indent=2))


if __name__ == '__main__':
    main()
