"""Finite exact checks; the accompanying note supplies the general proofs.

Run: python3 residual_counterexample_checks.py > residual_counterexample_results.json
Requires SymPy. No network or numerical floating-point inference is used.
"""

import json
from fractions import Fraction
from math import gcd

import sympy as sp


def atan_bounds(denominator, terms):
    total = sum(
        (Fraction((-1) ** k, (2 * k + 1) * denominator ** (2 * k + 1))
         for k in range(terms)), Fraction(0)
    )
    adjacent = total + Fraction(
        (-1) ** terms, (2 * terms + 1) * denominator ** (2 * terms + 1)
    )
    return min(total, adjacent), max(total, adjacent)


def run():
    x, y, z, t = sp.symbols("x y z t")
    matrices = []
    for degree in range(1, 8):
        highest_jet = degree + 1
        logarithm = sum(
            (-1) ** (k + 1) * t**k / sp.Integer(k)
            for k in range(1, highest_jet + 1)
        )
        columns = [("power", k) for k in range(highest_jet + 1)]
        columns += [("Y", 0), ("X1", 0)]
        rows = [(0, k) for k in range(highest_jet + 1)] + [(1, 1), (1, 0)]

        def entry(row, column):
            node, order = row
            kind, power = column
            if kind == "power":
                expression = (node * y + logarithm) ** power
            elif kind == "Y":
                expression = 1 + t
            else:
                expression = node * x + logarithm
            return sp.expand(expression).coeff(t, order)

        matrix = sp.Matrix([[entry(row, column) for column in columns] for row in rows])
        quotient = sum(y**k / sp.factorial(k) for k in range(1, degree + 1))
        determinant = sp.expand(matrix.det())
        assert sp.expand(determinant - (y - x) * quotient) == 0
        assert sp.rem(determinant, (y - x) ** 2, y) != 0
        assert sp.diff(quotient, y) == sum(
            y**k / sp.factorial(k) for k in range(degree)
        )
        matrices.append({"N": degree, "rows": len(rows), "identity": True,
                         "difference_multiplicity": 1})

    rational_centers = []
    for numerator, denominator in [(22, 7), (355, 113)]:
        center = Fraction(2 * numerator, denominator)
        a, b = center.numerator, center.denominator
        raw_numerator = 39 * b * b - a * a
        assert gcd(abs(raw_numerator), b) == 1
        base = Fraction(raw_numerator, b * b)
        assert 0 < abs(base) < 1
        for power in range(1, 9):
            assert (base**power).denominator == b ** (2 * power)
        rational_centers.append({"p": numerator, "q": denominator,
                                 "base": str(base), "checked_powers": list(range(1, 9)),
                                 "exact_reduced_denominator": True})

    lower5, upper5 = atan_bounds(5, 40)
    lower239, upper239 = atan_bounds(239, 40)
    lower_pi = 16 * lower5 - 4 * upper239
    upper_pi = 16 * upper5 - 4 * lower239
    assert 39 < 4 * lower_pi**2 < 4 * upper_pi**2 < Fraction(79, 2)

    matrix3 = sp.Matrix([[1, 0, 0], [0, 1, 0], [1, x, x + y - 2 * z]])
    polynomial3 = matrix3.det()
    assert polynomial3 == x + y - 2 * z
    assert polynomial3.subs({y: x, z: x}) == 0
    pair_checks = []
    for first, second in [(x, y), (x, z), (y, z)]:
        assert sp.expand(polynomial3.subs(first, second)) != 0
        pair_checks.append(str(first - second))

    return {
        "finite_interpolation_matrices": matrices,
        "primitive_polynomial_rational_centers": rational_centers,
        "machin_certificate": {"terms_per_arctangent": 40,
                               "certified": "0 < 4*pi^2 - 39 < 1/2",
                               "arithmetic": "Fraction, exact"},
        "three_variable_counterexample": {
            "polynomial": str(polynomial3), "diagonal_zero": True,
            "nondividing_pair_differences": pair_checks},
        "scope": "Finite checks only; not a full-row asymptotic construction or a P2848 proof."
    }


if __name__ == "__main__":
    print(json.dumps(run(), ensure_ascii=False, indent=2))
