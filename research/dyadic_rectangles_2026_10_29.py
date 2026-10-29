"""Exact checks for the dyadic factorial rectangle construction.

The pi interval is certified using alternating rational arctangent sums in
Machin's identity.  The finite-sum identity and all checks use integers and
fractions.  Run with: python dyadic_rectangles_2026_10_29.py
"""

from decimal import Decimal, localcontext
from fractions import Fraction
from math import comb, factorial, gcd


def atan_bounds(inverse, terms):
    assert inverse > 1 and terms > 0
    partial = sum(
        (Fraction((-1) ** j, (2 * j + 1) * inverse ** (2 * j + 1))
         for j in range(terms)), Fraction()
    )
    next_term = Fraction(1, (2 * terms + 1) * inverse ** (2 * terms + 1))
    return (partial, partial + next_term) if terms % 2 == 0 else (partial - next_term, partial)


def pi_bounds():
    # arctan(1) = 4 arctan(1/5) - arctan(1/239).
    a, b = atan_bounds(5, 200)
    c, d = atan_bounds(239, 80)
    return 16 * a - 4 * d, 16 * b - 4 * c


def v2(n):
    assert n > 0
    return (n & -n).bit_length() - 1


def bala_sum(n):
    """Return the integer sum from the finite identity used in the proof."""
    terms = []
    for k in range(1, n + 1):
        numerator = 2 ** (k - 1) * factorial(k - 1) * factorial(n + k - 1)
        value, rem = divmod(numerator, factorial(2 * k - 1))
        assert rem == 0
        assert v2(value) >= v2(factorial(n))
        terms.append(value)
    return sum(terms)


def rectangle_series_checks(limit=100):
    lower_pi, upper_pi = pi_bounds()
    assert lower_pi < upper_pi and upper_pi - lower_pi < Fraction(1, 10**200)
    raw_num, raw_den = 0, 1
    partial, term = Fraction(), Fraction(2)
    table = []
    for n in range(1, limit + 1):
        raw_num = (2 * n - 1) * raw_num + 2 * factorial(n - 1)
        raw_den *= 2 * n - 1
        partial += term
        term *= Fraction(n, 2 * n + 1)
        assert partial == Fraction(raw_num, raw_den)
        assert raw_den % 2 == 1 and raw_num == 2 * bala_sum(n)
        p, q = partial.numerator, partial.denominator
        assert p % (2 ** (1 + v2(factorial(n)))) == 0
        assert q > Fraction(2 ** (n - 1), n + 1)

        # The positive tail is between the next term and twice that term.
        error_lo, error_hi = lower_pi - partial, upper_pi - partial
        assert term < error_lo < error_hi < 2 * term
        lower = Fraction(2 ** (n - 1), (n + 1) ** 2 * (2 * n + 1))
        assert q * q * error_lo > lower

        common_factor = gcd(raw_num, raw_den)
        analytic = Fraction(2 * factorial(2 * n), 2**n * (2 * n + 1) * common_factor**2)
        assert analytic == q * q * term
        if n in (1, 2, 3, 4, 8, 16, 32):
            with localcontext() as ctx:
                ctx.prec = 25
                err = Decimal(error_lo.numerator) / Decimal(error_lo.denominator)
                table.append((n, str(partial), f'{err:.8E}', f'{q*q*err:.8E}'))
    return table


def integral_rational(n):
    """Polynomial division for x^(4n)(1-x)^(4n)/(1+x^2)."""
    polynomial = [0] * (8 * n + 1)
    for j in range(4 * n + 1):
        polynomial[4 * n + j] = (-1) ** j * comb(4 * n, j)
    quotient = [0] * (8 * n - 1)
    for degree in range(8 * n, 1, -1):
        coefficient = polynomial[degree]
        quotient[degree - 2] = coefficient
        polynomial[degree] = 0
        polynomial[degree - 2] -= coefficient
    assert polynomial[1] == 0 and polynomial[0] == (-4) ** n
    rational_integral = sum(
        (Fraction(coefficient, degree + 1) for degree, coefficient in enumerate(quotient)),
        Fraction(),
    )
    return -4 * rational_integral / (-4) ** n


def integral_checks():
    lower_pi, upper_pi = pi_bounds()
    for n in range(1, 9):
        approximation = integral_rational(n)
        if n % 2:
            error_lo, error_hi = approximation - upper_pi, approximation - lower_pi
        else:
            error_lo, error_hi = lower_pi - approximation, upper_pi - approximation
        assert 0 < error_lo < error_hi < Fraction(4, 1024**n)
    assert integral_rational(1) == Fraction(22, 7)
    assert integral_rational(2) == Fraction(47171, 15015)
    assert integral_rational(3) == Fraction(431302721, 137287920)


if __name__ == '__main__':
    rows = rectangle_series_checks()
    integral_checks()
    print('All exact checks passed; n=1..100 for the binary series, n=1..8 for the 1024 family.')
    print('N | P_N | certified-error approximation | q^2 times error approximation')
    for row in rows:
        print(' | '.join(map(str, row)))
