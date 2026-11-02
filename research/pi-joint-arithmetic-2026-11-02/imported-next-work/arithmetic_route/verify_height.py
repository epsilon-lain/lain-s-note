"""Finite exact checks supplementing the general REFLECTION_HEIGHT.md proof."""
import json
import math
from pathlib import Path
import sympy as s

HERE = Path(__file__).resolve().parent


def main():
    dimensions = []
    for n in range(4, 42, 2):
        r = n // 2
        pascal = s.Matrix([[s.binomial(i, j) for j in range(n)] for i in range(n)])
        shift = s.zeros(n)
        for i in range(n - 1):
            shift[i, i + 1] = 1
        conjugate = pascal * shift * pascal.inv()
        assert all(conjugate[i, j] == int(j == i + 1) - int(j == i)
                   for i in range(n - 1) for j in range(n))
        assert conjugate[n - 1, n - 1] == n - 1
        constants = []
        for sign in [1, -1]:
            matrix = s.Matrix([
                [s.binomial(i, j) if (i % 2 == 0) == (sign == 1)
                 else s.binomial(i, j - 1) for j in range(n)]
                for i in range(n)
            ])
            value = int(matrix.det())
            expected = (-1) ** (r - 1) * (n - 1) if sign == 1 else (-1) ** r
            assert value == expected
            constants.append(value)
        dimensions.append({'n': n, 'c_plus': constants[0], 'c_minus': constants[1]})

    fixtures = json.loads((HERE / 'factors.json').read_text())
    for row in fixtures['records']:
        n = row['n']
        B = int(row['normalization_B'])
        S = n * (n - 1) // 2
        for factor in row['factors']:
            coefficients = {int(k): int(v) for k, v in factor['coefficients'].items()}
            assert abs(coefficients[0]) == B * ((n - 1) if factor['sign'] == '+' else 1)
            # Clear the exact rational factor (3/2)^S.
            assert sum(abs(v) for v in coefficients.values()) * 2 ** S <= (
                B * math.factorial(n) * 2 ** n * 3 ** S
            )
    result = {'status': 'passed', 'pascal_dimensions': dimensions,
              'coefficient_fixtures_checked': [row['n'] for row in fixtures['records']],
              'scope': 'Finite identities and coefficient checks; general height theorem is analytic.'}
    (HERE / 'height_verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print('Passed: Pascal constants n=4..40; coefficient fixtures and height bounds n=4..12.')


if __name__ == '__main__':
    main()
