#!/usr/bin/env python3
"""Check exact pole-degree predictions by independent rational-function arithmetic."""
from fractions import Fraction
from math import gcd

from sympy import QQ
from sympy.polys.fields import field

from eds_primitive_denominator_experiment import (
    ALT,
    binomial_weights,
    point_mul,
    xi_point,
)


def pole_orders(weights):
    last = max(weights)
    return [
        r for r in range(4, 2 * last + 1, 2)
        if sum(Fraction(weight, j) for j, weight in weights.items()
               if j % r == r // 2)
    ]


def jordan_two(n):
    return sum(gcd(gcd(a, b), n) == 1 for a in range(n) for b in range(n))


def symbolic_ratios(largest):
    # For the generic point (x,y), store [n](x,y) as (f_n(x), y*g_n(x)).
    rational_functions, x = field("x", QQ)
    cubic = x**3 - 2
    coordinates = {1: (x, rational_functions.one)}
    for j in range(2, largest + 1):
        previous_x, previous_y_factor = coordinates[j - 1]
        if j == 2:
            slope_factor = 3 * x**2 / (2 * cubic)
        else:
            slope_factor = (previous_y_factor - 1) / (previous_x - x)
        next_x = cubic * slope_factor**2 - previous_x - x
        next_y_factor = slope_factor * (x - next_x) - 1
        coordinates[j] = next_x, next_y_factor
    return {
        j: f_x / (x * f_y)
        for j, (f_x, f_y) in coordinates.items()
    }


def main():
    special = {2: 1, 6: -3}
    sparse = {2: -364, 3: 291, 4: -36}
    kernel = {2: 962, 3: -768, 4: 95}
    for moment in (1, 7):
        assert sum(weight * j**moment for j, weight in sparse.items()) == 1
        assert sum(weight * j**moment for j, weight in kernel.items()) == 0
    first_coefficient = Fraction(48, 91) * (
        sum(weight * j**13 for j, weight in sparse.items()) - 1
    )
    step = Fraction(48, 91) * sum(weight * j**13 for j, weight in kernel.items())
    assert (first_coefficient, step) == (-1031184000, 2721116160)
    assert abs(first_coefficient) < step / 2

    examples = [
        (m, label, weights)
        for m in (4, 5, 6)
        for label, weights in (
            ("binomial", binomial_weights(m)),
            ("alternative", ALT[m]),
        )
    ]
    symbolic = symbolic_ratios(max(max(weights) for _, _, weights in examples))
    examples += [(7, "binomial", binomial_weights(7)),
                 (12, "CM_sparse", sparse),
                 (0, "cancel_order_4", special)]
    for m, label, weights in examples:
        poles = pole_orders(weights)
        predicted = sum(jordan_two(r) for r in poles)
        value = sum(weight * symbolic[j] for j, weight in weights.items())
        map_degree = max(value.numer.degree(), value.denom.degree())
        assert 2 * map_degree == predicted
        if label == "binomial":
            assert poles == list(range(4, 2 * (m + 2) + 1, 2))
        if label == "CM_sparse":
            assert predicted == 84
        if label == "cancel_order_4":
            assert poles == [12]
        print(f"m={m} {label}: {poles}, height coefficient {predicted}")

    # At p=5, P=(3,5) reduces to a nonzero 2-torsion point.
    # The ratio for multiplier 5 nevertheless has a 5 in its denominator.
    ratio = xi_point(point_mul(5)) / xi_point(point_mul(1))
    assert ratio.numerator % 5 and ratio.denominator % 5 == 0
    assert ratio.denominator % 25
    print(f"p=5, s=1, j=5: reduced denominator has 5-adic valuation 1 "
          f"(denominator={ratio.denominator})")


if __name__ == "__main__":
    main()
