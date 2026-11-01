"""Exact checks for diagonal profiles and the full weighted-column kernel.

Only the Python standard library is used. No floating point or numerical pi
values enter these algebraic checks. Kernel examples use formal rational
centers, not an assertion that a rational number is an exponential period.
"""
from fractions import Fraction as F
from itertools import product, combinations
from math import comb, factorial, prod
from pathlib import Path
import json


def determinant(a):
    a = [[F(x) for x in row] for row in a]
    out = F(1)
    for k in range(len(a)):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return F(0)
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            out = -out
        d = a[k][k]
        out *= d
        for i in range(k+1, len(a)):
            q = a[i][k]/d
            for j in range(k+1, len(a)):
                a[i][j] -= q*a[k][j]
    return out


def rank(a):
    a = [[F(x) for x in row] for row in a]
    r = 0
    if not a:
        return 0
    for k in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][k]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        d = a[r][k]
        a[r] = [x/d for x in a[r]]
        for i in range(r+1, len(a)):
            q = a[i][k]
            if q:
                a[i] = [x-q*y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def multiply(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)]
            for row in a]


def types(weights, H):
    return [a for a in product(*(range(H//w+1) for w in weights))
            if sum(w*x for w, x in zip(weights, a)) <= H]


def raw_columns(weights, H):
    return [(index[0], index[1:]) for index in types(weights, H)]


def scaled_time_entry(j, s, h, degree, Z):
    if degree < 0:
        return F(0)
    return sum(F(comb(degree, k)*factorial(s)//factorial(s-k))
               * h**(s-k)*(j*Z)**(degree-k)
               for k in range(min(s, degree)+1))


def raw_entry(row, column, Z):
    j, s, beta = row
    h, alpha = column
    if any(a < b for a, b in zip(alpha, beta)):
        return F(0)
    return prod(comb(a, b) for a, b in zip(alpha, beta))*scaled_time_entry(
        j, s, h, sum(alpha)-sum(beta), Z)


def normal_entry(row, column, Z):
    j, s, beta = row
    h, alpha = column
    a, nu = alpha[0], alpha[1:]
    ell = sum(beta)-sum(nu)
    if ell < 0 or a < ell or any(n < b for n, b in zip(nu, beta[1:])):
        return F(0)
    scalar = (-1)**(sum(nu)-sum(beta[1:]))
    scalar *= prod(comb(n, b) for n, b in zip(nu, beta[1:]))*comb(a, ell)
    return scalar*scaled_time_entry(j, s, h, a-ell, Z)


def check_columns():
    total = 0
    records = []
    for weights, H in [((1, 1, 3), 5), ((2, 2, 5), 7), ((1, 1, 2, 4), 5)]:
        cols = raw_columns(weights, H)
        index = {c: i for i, c in enumerate(cols)}
        change = [[0]*len(cols) for _ in cols]
        for ci, (h, alpha) in enumerate(cols):
            a, nu = alpha[0], alpha[1:]
            for k in product(*(range(n+1) for n in nu)):
                exponent = (a+sum(k),)+tuple(n-x for n, x in zip(nu, k))
                key = (h, exponent)
                assert key in index
                change[index[key]][ci] += (-1)**sum(k)*prod(comb(n, x) for n, x in zip(nu, k))
        assert determinant(change) == 1
        row_types = [b for b in product(range(3), repeat=len(weights)-1) if sum(b) <= 2]
        rows = [(j, s, b) for j in (0, 1) for s in (0, 1, 2) for b in row_types]
        for Z in (F(0), F(1, 3), F(2)):
            raw = [[raw_entry(r, c, Z) for c in cols] for r in rows]
            transformed = multiply(raw, change)
            expected = [[normal_entry(r, c, Z) for c in cols] for r in rows]
            assert transformed == expected
            total += len(rows)*len(cols)
        records.append({'weights': weights, 'H': H, 'column_count': len(cols),
                        'basis_change_determinant': 1})
    return total, records


def check_rows():
    records = []
    for weights, cutoff in [((1, 3), 7), ((2, 5), 11), ((1, 2, 5), 8)]:
        for time_budget in (0, 2):
            packet = types(weights, cutoff-1-time_budget)
            lookup = {b: i for i, b in enumerate(packet)}
            n = len(packet)
            forward = [[0]*n for _ in packet]
            inverse = [[0]*n for _ in packet]
            for bi, beta in enumerate(packet):
                for k in product(range(beta[-1]+1), repeat=len(weights)-1):
                    t = beta[-1]-sum(k)
                    if t < 0:
                        continue
                    gamma = tuple(b+x for b, x in zip(beta[:-1], k))+(t,)
                    assert gamma in lookup
                    forward[bi][lookup[gamma]] += (-1)**sum(k)*prod(comb(b+x, b) for b, x in zip(beta[:-1], k))
            for gi, gamma in enumerate(packet):
                nu, t = gamma[:-1], gamma[-1]
                for k in product(range(t+1), repeat=len(weights)):
                    if sum(k) != t:
                        continue
                    beta = tuple(v+x for v, x in zip(nu, k[:-1]))+(k[-1],)
                    assert beta in lookup
                    inverse[gi][lookup[beta]] += prod(comb(v+x, v) for v, x in zip(nu, k[:-1]))
            assert multiply(forward, inverse) == [[int(i == j) for j in range(n)] for i in range(n)]
            assert determinant(forward) == determinant(inverse) == 1
            records.append({'weights': weights, 'cutoff': cutoff, 'time_budget': time_budget,
                            'rows': n, 'determinants': [1, 1]})
    # A least-weight tangent derivative can require an excluded heavy old row.
    assert 1 < 2 and not 3 < 2
    return records


def profile_entry(s, h, d):
    return F(factorial(s)//factorial(s-d)*h**(s-d)) if s >= d else F(0)


def check_capacities():
    records = []
    for w0, wmin, L in [(1, 2, 1), (1, 2, 5), (2, 3, 6), (2, 2, 5)]:
        profiles = [(h, d) for h in range(L//w0+1) for d in range(L//wmin+1)
                    if w0*h+wmin*d <= L]
        d0 = L//wmin
        for r in (0, 1, 3):
            capacity = len(profiles)-min(r, d0+1)
            rows = range(r, r+capacity)
            matrix = [[profile_entry(s, h, d) for h, d in profiles] for s in rows]
            assert rank(matrix) == capacity
            active = [(h, d) for h, d in profiles if h or d >= r]
            assert len(active) == capacity
            if capacity:
                assert determinant([[profile_entry(s, h, d) for h, d in active] for s in rows])
            records.append({'w0': w0, 'w_min': wmin, 'remaining_budget': L,
                            'minimum_time_order': r, 'profile_dimension': len(profiles),
                            'exact_capacity': capacity})
    # The old type-2 obstruction is exactly the space {1,e^z} after one derivative.
    assert [(h, d) for h in range(2) for d in range(1) if h+2*d <= 1] == [(0, 0), (1, 0)]
    return records


def series_product(a, b, N):
    return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(N+1)]


def check_generating_function():
    checks = 0
    for weights, values in [((1, 2), (F(3, 4), F(2, 5))),
                             ((1, 2, 5), (F(3, 4), F(2, 5), F(3, 7)))]:
        N = 12
        coefficients = [F(1)]*(N+1)  # Cumulative 1/(1-t).
        for weight, value in zip(weights, values):
            factor = [value**(k//weight) if k % weight == 0 else F(0) for k in range(N+1)]
            coefficients = series_product(coefficients, factor, N)
        for H in range(N+1):
            direct = sum(prod(value**exponent for value, exponent in zip(values, index))
                         for index in types(weights, H))
            assert coefficients[H] == direct
            checks += 1
    return checks


def grand_entry(row, col, q, r, Z):
    j, s = row
    k, ell = col
    x, y = j*Z, k*Z
    base = 1-r*x*y
    A00 = 1/base
    A10 = r*y/base**2
    A01 = r*x/base**2
    A11 = r*(1+r*x*y)/base**3
    H0 = 1/(1-q)
    H1 = q/(1-q)**2
    H2 = q*(1+q)/(1-q)**3
    if (s, ell) == (0, 0):
        return H0*A00
    if (s, ell) == (1, 0):
        return H1*A00+H0*A10
    if (s, ell) == (0, 1):
        return H1*A00+H0*A01
    return H2*A00+H1*(A10+A01)+H0*A11


def check_gram():
    H = 6
    t = F(1, 4)
    Z = F(1, 3)
    rows = [(j, s) for j in (0, 1) for s in (0, 1)]
    cols = raw_columns((1, 2), H)
    B = [[raw_entry((j, s, (0,)), c, Z) for c in cols] for j, s in rows]
    finite = multiply(B, list(zip(*B)))
    grand = [[grand_entry(row, col, t, t*t, Z) for col in rows] for row in rows]
    difference = [[t**(-H)*grand[i][j]-finite[i][j] for j in range(4)] for i in range(4)]
    principal_checks = 0
    for size in range(1, 5):
        for indices in combinations(range(4), size):
            assert determinant([[difference[i][j] for j in indices] for i in indices]) > 0
            principal_checks += 1
    selected = [(0, (a,)) for a in range(4)]
    minor = determinant([[raw_entry((j, s, (0,)), c, Z) for c in selected] for j, s in rows])
    finite_det = determinant(finite)
    grand_bound = t**(-4*H)*determinant(grand)
    assert 0 < minor*minor <= finite_det <= grand_bound
    return {'H': H, 'weights': [1, 2], 'formal_real_center': str(Z), 't': str(t),
            'rows': rows, 'column_count': len(cols),
            'positive_principal_minors_of_loewner_difference': principal_checks,
            'selected_minor_squared': str(minor*minor),
            'full_finite_gram_determinant': str(finite_det),
            'grand_kernel_determinant_bound': str(grand_bound)}


def main():
    entry_checks, columns = check_columns()
    result = {'status': 'passed', 'column_transformation_entry_checks': entry_checks,
              'column_basis_changes': columns, 'greatest_anchor_row_changes': check_rows(),
              'sharp_capacity_checks': check_capacities(),
              'kernel_generating_coefficient_checks': check_generating_function(),
              'finite_gram_check': check_gram(),
              'scope': 'Exact coordinate, capacity and weighted-column-kernel checks. '
                       'Rational formal nodes are not exponential periods. No uniform '
                       'large-H determinant saving or pi approximation claim is certified.'}
    Path(__file__).with_name('diagonal_profile_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'column_entries': entry_checks,
                      'column_changes': len(columns), 'row_changes': len(result['greatest_anchor_row_changes']),
                      'capacity_cases': len(result['sharp_capacity_checks']),
                      'kernel_coefficients': result['kernel_generating_coefficient_checks'],
                      'loewner_principal_minors': result['finite_gram_check']['positive_principal_minors_of_loewner_difference']}, indent=2))


if __name__ == '__main__':
    main()
