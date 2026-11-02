"""Exact checks for the primitive rational-center comparison (standard library).

Input: a small coefficient fixture copied byte-for-value from the preceding
published raw-polynomial certificates. The general theorem is in the note.
"""
from fractions import Fraction as F
from math import factorial, gcd, lcm, prod
from pathlib import Path
import hashlib
import json


def A(n):
    return prod(factorial(j) for j in range(n))


def at_imaginary(poly, x):
    assert all(k % 2 == 0 for k in poly)
    return sum(c * (-1)**(k//2) * x**k for k, c in poly.items())


def derivative_majorant(poly, radius):
    return sum(k * abs(c) * radius**(k-1) for k, c in poly.items() if k)


def exact_sufficient_height(b):
    # Choosing 24b does not prove a uniform threshold; each recorded case
    # is certified by its exact rational comparison below.
    n = max(4, 24*b)
    n += n % 2
    while True:
        d = A(n)**2 // A(n//2)**4
        N = n*n//2
        ratio = 2 * b**N * F(44, 7)**(n*n) * F(d*A(n)**2, A(2*n))
        if ratio <= 1:
            break
        n += 2
    return {
        'b': b, 'certified_even_n': n, 'N_n': N,
        'twice_baseline_to_rational_floor_at_most_one': ratio <= 1,
        'ratio_numerator_bits': ratio.numerator.bit_length(),
        'ratio_denominator_bits': ratio.denominator.bit_length(),
        'primitive_denominator_bits': d.bit_length(),
        'rational_value_denominator_bits': (b**N).bit_length(),
    }


def main():
    fixture_path = Path(__file__).with_name('even_period_coefficient_fixture.json')
    fixture = json.loads(fixture_path.read_text())
    records = []
    checks = 0
    for case in fixture['cases']:
        n = case['n']
        poly = {int(k): F(v) for k, v in case['polynomial'].items()}
        d = A(n)**2 // A(n//2)**4
        N = n*n//2
        D = n*(n-1)
        assert d == lcm(*(c.denominator for c in poly.values()))
        assert max(poly) == N and poly[N]*d == (-1)**(n//2)
        assert all(k % 2 == 0 for k in poly)
        integer_poly = {k: int(c*d) for k, c in poly.items()}
        assert all(c*d == integer_poly[k] for k, c in poly.items())
        assert gcd(*integer_poly.values()) == 1
        for b in range(2, 16):
            for a in range(6*b, 7*b+1):
                if gcd(a, b) != 1:
                    continue
                value = at_imaginary(integer_poly, F(a, b))
                assert value and value.denominator == b**N
                assert abs(value) >= F(1, b**N)
                checks += 1
        named = []
        for p, q in ((22,7), (355,113), (19,6), (25,8)):
            x = F(2*p, q)
            value = at_imaginary(integer_poly, x)
            assert value.denominator == x.denominator**N
            named.append({'p':p, 'q':q, 'b':x.denominator,
                          'value_denominator':str(value.denominator)})
            checks += 1
        # Coefficient majorant proves the finite derivative bound throughout
        # the disk, rather than sampling a handful of derivative values.
        majorant = derivative_majorant(poly, F(8))
        determinant_bound = factorial(2*n)*D*9**(D-1)
        assert majorant <= determinant_bound <= factorial(2*n)*10**D
        records.append({'n':n, 'exact_least_coefficient_denominator':str(d),
                        'degree':N, 'primitive_unit_leading_coefficient':True,
                        'derivative_coefficient_majorant_check':True,
                        'named_rational_centers':named})
    heights = [exact_sufficient_height(b) for b in (2,3,7,10)]
    result = {
        'status':'passed',
        'fixture_sha256':hashlib.sha256(fixture_path.read_bytes()).hexdigest(),
        'coefficient_families':len(records),
        'exact_rational_value_denominator_checks':checks,
        'polynomial_checks':records,
        'exact_sufficient_height_certificates':heights,
        'scope':'Finite exact identities and sufficient comparisons. General denominator, '
                'derivative, asymptotic threshold and error bound are proved in the note. '
                'No bad-approximability or counterexample claim.'
    }
    Path(__file__).with_name('rational_specialization_verification.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('polynomial_checks',)},indent=2))


if __name__ == '__main__':
    main()
