"""Standalone exact certificates for genuine period-polynomial root geometry.

Standard library only. All roots claimed in the note follow from the exact
linear/quadratic factors checked here, not numerical root approximations.
"""
from fractions import Fraction
from itertools import combinations
from math import comb, factorial, prod
from pathlib import Path
import json


def trim(poly):
    while len(poly)>1 and not poly[-1]:
        poly.pop()
    return poly


def add(a, b):
    out = [Fraction(0)]*max(len(a), len(b))
    for i, value in enumerate(a):
        out[i] += value
    for i, value in enumerate(b):
        out[i] += value
    return trim(out)


def scale(poly, value):
    return trim([x*value for x in poly])


def evaluate(poly, value):
    out = 0
    for coefficient in reversed(poly):
        out = out*value+coefficient
    return out


def interpolate(values):
    basis, out = [Fraction(1)], [Fraction(0)]
    differences = list(map(Fraction, values))
    for k in range(len(values)):
        out = add(out, scale(basis, differences[0]))
        differences = [b-a for a, b in zip(differences, differences[1:])]
        next_basis = [Fraction(0)]*(len(basis)+1)
        for i, value in enumerate(basis):
            next_basis[i] -= value*k/(k+1)
            next_basis[i+1] += value/(k+1)
        basis = next_basis
    return trim(out)


def entry(row, column, z):
    j, s, beta = row
    h, alpha = column
    if beta>alpha:
        return 0
    return sum(comb(alpha, beta)*comb(alpha-beta, k)
               *(factorial(s)//factorial(s-k))*h**(s-k)
               *(j*z)**(alpha-beta-k)
               for k in range(min(s, alpha-beta)+1))


def determinant(matrix):
    a = [row[:] for row in matrix]
    sign, previous = 1, 1
    for k in range(len(a)-1):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        diagonal = a[k][k]
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                numerator = a[i][j]*diagonal-a[i][k]*a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator//previous
            a[i][k] = 0
        previous = diagonal
    return sign*a[-1][-1]


def polynomial(rows, columns):
    bound = sum(a for h, a in columns)-sum(b for j, s, b in rows)
    denominator = prod(factorial(s) for j, s, b in rows)
    def det_at(z):
        return Fraction(determinant([[entry(r, c, z) for c in columns] for r in rows]),
                        denominator)
    poly = interpolate([det_at(z) for z in range(bound+1)])
    assert evaluate(poly, -1) == det_at(-1)
    assert evaluate(poly, bound+2) == det_at(bound+2)
    return poly


def sparse(poly):
    return {str(i): str(value) for i, value in enumerate(poly) if value}


def from_sparse(values):
    out = [Fraction(0)]*(max(values)+1)
    for i, value in values.items():
        out[i] = Fraction(value)
    return out


def at_rational_imaginary(poly, numerator=44, denominator=7):
    x = Fraction(numerator, denominator)
    real = sum(value*x**i*(-1)**(i//2)
               for i, value in enumerate(poly) if i % 2 == 0)
    imaginary = sum(value*x**i*(-1)**((i-1)//2)
                    for i, value in enumerate(poly) if i % 2)
    return [str(real), str(imaginary)], bool(real or imaginary)


def certificate(name, hmax, v0, v1, frequencies, expected, root_certificate):
    rows = [(j, s, b) for j in range(2) for b in range(hmax//v1+1)
            for s in range(hmax//v0+1) if v0*s+v1*b<hmax]
    columns = [(h, alpha) for alpha, heights in enumerate(frequencies) for h in heights]
    assert len(rows) == len(columns)
    assert len(set(columns)) == len(columns)
    assert all(h+2*alpha<=hmax for h, alpha in columns)
    computed = polynomial(rows, columns)
    assert computed == from_sparse(expected)
    rational_value, nonzero = at_rational_imaginary(computed)
    assert nonzero
    degree_surplus = sum(a for h, a in columns)
    reciprocal = list(reversed(computed+[Fraction(0)]*(degree_surplus+1-len(computed))))
    return {'name': name, 'H': hmax, 'K': 2, 'w0': 1, 'w1': 2,
            'v0': v0, 'v1': v1, 'theta': str(Fraction(2, v1)),
            'source_rows': rows, 'raw_columns': columns,
            'D': degree_surplus, 'period_polynomial': sparse(computed),
            'reciprocal_polynomial': sparse(reciprocal),
            'period_polynomial_at_44i_over_7': rational_value,
            'root_certificate': root_certificate}


def exhaustive_six_row_family():
    rows = [(j, s, 0) for j in range(2) for s in range(3)]
    cases = []
    for k, l in combinations(range(5), 2):
        for a in range(5):
            columns = [(0, 0), (1, 0), (2, 0), (k, 1), (l, 1), (a, 2)]
            assert all(h+2*alpha<=8 for h, alpha in columns)
            computed = polynomial(rows, columns)
            expected = from_sparse({3: Fraction((l-k)*(4*a-2*k-2*l), 2),
                                    4: Fraction((l-k)*(a-k)*(a-l), 2)})
            assert computed == trim(expected)
            maximum = max(h for h, alpha in columns)
            reflected = [(maximum-h, alpha) for h, alpha in columns]
            assert all(h+2*alpha<=8 for h, alpha in reflected)
            reflected_poly = polynomial(rows, reflected)
            assert reflected_poly == [value*(-1)**i for i, value in enumerate(computed)]
            cases.append({'k': k, 'l': l, 'a': a, 'polynomial': sparse(computed),
                          'reflection_checked': True})
    assert len(cases) == 50
    return {'complete_packet': {'K': 2, 'H': 8, 'w0': 1, 'w1': 2, 'v0': 3, 'v1': 10},
            'case_count': len(cases), 'cases': cases}


def main():
    # The sign certificates below use only integer comparisons.
    assert 57**2 < 7665 < 88**2
    assert 148**2-4*9*1604 == -35840 < 0
    assert 335**2-4*36*380 == 57505 > 0
    assert 57505 < 335**2
    cases = [
        certificate('negative_real_root', 6, 2, 8, [[0, 1, 2], [0, 1], [2]],
                    {3: 3, 4: 1}, {'nonzero_roots': ['-3']}),
        certificate('positive_real_root', 6, 2, 8, [[0, 1, 2], [1, 2], [0]],
                    {3: -3, 4: 1}, {'nonzero_roots': ['3']}),
        certificate('six_row_imaginary_axis', 8, 3, 10, [[1, 3], [0, 4], [1, 3]],
                    {3: -32, 5: -24}, {'nonzero_roots': ['2i/sqrt(3)', '-2i/sqrt(3)']}),
        certificate('eight_row_imaginary_axis', 8, 2, 10, [[0, 1, 2, 3], [0, 3], [1, 2]],
                    {4: 3, 6: 1}, {'nonzero_roots': ['i*sqrt(3)', '-i*sqrt(3)']}),
        certificate('eight_row_opposite_real', 8, 2, 10, [[0, 1, 2, 3], [1, 2], [0, 3]],
                    {4: -13, 6: 1}, {'nonzero_roots': ['sqrt(13)', '-sqrt(13)']}),
        certificate('twelve_row_four_axis_rays', 12, 2, 14,
                    [list(range(6)), [1, 4], [0, 2, 3, 5]],
                    {6: -69, 8: Fraction(-57, 4), 10: 1},
                    {'squared_roots': ['(57+sqrt(7665))/8', '(57-sqrt(7665))/8'],
                     'sign_certificate': '57^2<7665, so one squared root is positive and one negative',
                     'halfplane_consequence': 'Nonzero roots meet all four coordinate half-rays'}),
        certificate('twelve_row_four_open_quadrants', 12, 2, 14,
                    [list(range(6)), [2, 3], [0, 1, 4, 5]],
                    {6: Fraction(1604, 9), 8: Fraction(-148, 9), 10: 1},
                    {'quadratic_in_Z_squared_discriminant': -35840,
                     'consequence': 'Two nonreal conjugate squared roots; their square roots lie in all four open quadrants'}),
        certificate('twelve_row_two_imaginary_pairs', 12, 2, 14,
                    [list(range(6)), [0, 5], [1, 2, 3, 4]],
                    {6: Fraction(95, 9), 8: Fraction(335, 36), 10: 1},
                    {'squared_roots': ['(-335+sqrt(57505))/72', '(-335-sqrt(57505))/72'],
                     'sign_certificate': '0<57505<335^2, so both squared roots are negative'}),
    ]
    family = exhaustive_six_row_family()
    result = {'status': 'passed', 'explicit_cases': cases, 'exhaustive_six_row_family': family,
              'scope': 'Exact genuine complete packets and raw weighted columns. '
                       'Obstructs universal half-plane, narrow-sector, and imaginary-axis '
                       'zero-free assertions; does not show roots approach 2*pi*i or '
                       'that actual-period determinants become arbitrarily small.'}
    Path(__file__).with_name('period_root_geometry_verification.json').write_text(
        json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'explicit_certificates': len(cases),
                      'exhaustive_family_size': family['case_count'],
                      'polynomials': {c['name']: c['period_polynomial'] for c in cases}}, indent=2))


if __name__ == '__main__':
    main()
