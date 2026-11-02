"""Exact finite irreducibility certificates for the existing n<=12 factors.

At a good prime a monic Q-factor must reduce to a product of modular
irreducible factors of the same degree. Intersect subset-sum degree sets.
An empty proper-degree set certifies Q-irreducibility. This makes no
assertion for dimensions beyond the stored coefficient fixtures.
"""
import json
import warnings
from pathlib import Path
import sympy as s
from sympy.utilities.exceptions import SymPyDeprecationWarning

HERE = Path(__file__).resolve().parent


def main():
    warnings.filterwarnings('ignore', category=SymPyDeprecationWarning)
    z = s.symbols('z')
    data = json.loads((HERE / 'factors.json').read_text())
    records = []
    for row in data['records']:
        for factor in row['factors']:
            polynomial = s.Poly(sum(int(v) * z ** int(k)
                                   for k, v in factor['coefficients'].items()), z)
            possible = set(range(1, polynomial.degree()))
            certificates = []
            for prime in s.primerange(2, 300):
                if s.gcd(polynomial, polynomial.diff(), modulus=prime).degree():
                    continue
                unit, factors = s.factor_list(polynomial, modulus=prime)
                assert unit % prime == 1
                product = s.Poly(1, z, modulus=prime)
                degrees = []
                modular_factors = []
                for g, multiplicity in factors:
                    assert multiplicity == 1 and g.is_irreducible
                    product *= g
                    degrees.append(g.degree())
                    modular_factors.append([int(c) % prime for c in g.all_coeffs()])
                assert product == s.Poly(polynomial, modulus=prime)
                sums = {0}
                for degree in degrees:
                    sums |= {value + degree for value in list(sums)}
                previous = set(possible)
                possible &= sums
                if previous != possible:
                    certificates.append({'prime': prime, 'degrees': degrees,
                                         'monic_factors_coefficients_descending': modular_factors,
                                         'remaining_proper_degrees': sorted(possible)})
                if not possible:
                    break
            assert not possible
            records.append({'n': row['n'], 'sign': factor['sign'],
                            'degree': polynomial.degree(), 'Q_irreducible': True,
                            'modular_certificates': certificates})
    result = {'status': 'passed', 'records': records,
              'scope': 'Finite exact irreducibility only; no general irreducibility theorem.'}
    (HERE / 'modular_irreducibility.json').write_text(json.dumps(result, indent=2) + '\n')
    print('Passed: all ten reflection factors n=4,6,8,10,12 are Q-irreducible.')


if __name__ == '__main__':
    main()
