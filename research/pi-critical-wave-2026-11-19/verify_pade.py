"""Exact arithmetic checks for the diagonal exponential Padé detector."""

from fractions import Fraction
from math import factorial, gcd


def values(degree, numerator, denominator):
    real, imaginary = 0, 0
    for order in range(degree + 1):
        coefficient = factorial(2 * degree - order) // (
            factorial(degree - order) * factorial(order)
        )
        term = coefficient * numerator**order * denominator**(degree - order)
        if order % 4 == 0:
            real += term
        elif order % 4 == 1:
            imaginary += term
        elif order % 4 == 2:
            real -= term
        else:
            imaginary -= term
    return real, imaginary


def check_primitive_area():
    tested = 0
    for denominator in range(2, 16):
        for numerator in range(3 * denominator, 4 * denominator + 1):
            if gcd(numerator, denominator) != 1:
                continue
            left = values(0, numerator, denominator)
            for degree in range(1, 19):
                right = values(degree, numerator, denominator)
                a, b = left
                c, d = right
                assert a * d - b * c == numerator ** (2 * degree - 1)
                if degree % 2 == 0:
                    assert Fraction(c, denominator**degree).denominator == denominator**degree
                left = right
                tested += 1
    return tested


def check_positive_period_bounds():
    # A decimal enclosure from the alternating arctangent series and Machin's identity.
    def atan_bounds(reciprocal, terms):
        partial = sum(
            Fraction((-1)**k, (2*k + 1) * reciprocal**(2*k + 1))
            for k in range(terms)
        )
        next_partial = partial + Fraction(
            (-1)**terms, (2*terms + 1) * reciprocal**(2*terms + 1)
        )
        return sorted((partial, next_partial))

    five_lo, five_hi = atan_bounds(5, 180)
    other_lo, other_hi = atan_bounds(239, 80)
    pi_lo, pi_hi = 16*five_lo - 4*other_hi, 16*five_hi - 4*other_lo
    scale = 10**230
    lo = Fraction(int(pi_lo*scale), scale)
    hi = Fraction(int(pi_hi*scale) + 1, scale)
    assert 3 < lo < pi_lo < pi_hi < hi < 4
    for degree in range(2, 42, 2):
        lower, upper = Fraction(), Fraction()
        constant = factorial(degree)
        for order in range(0, degree + 1, 2):
            coefficient = factorial(2 * degree - order) // (
                factorial(degree - order) * factorial(order)
            )
            term_a = (-1)**(order // 2) * coefficient * lo**order
            term_b = (-1)**(order // 2) * coefficient * hi**order
            lower += min(term_a, term_b)
            upper += max(term_a, term_b)
        assert lower > lo**(2*degree + 1) * constant / (4 * factorial(2*degree + 1))
        assert upper < hi**(2*degree + 1) * constant / (2 * factorial(2*degree + 1))
    return 20


if __name__ == "__main__":
    print(f"Checked {check_primitive_area()} adjacent integer Padé determinants and denominators.")
    print(f"Verified {check_positive_period_bounds()} positive period residuals with a rational π enclosure.")
