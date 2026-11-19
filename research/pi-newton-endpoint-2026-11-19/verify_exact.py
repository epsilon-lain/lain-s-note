"""Independent exact full-determinant audit. Requires SymPy.

Usage: python verify_exact.py [path/to/verify_newton.py]
Only the ordered row/column labels and selection of the common square minor
are shared with the modular test. The entries here use the original t variable.
"""
from functools import lru_cache
from math import comb, factorial, gcd, prod
from pathlib import Path
from runpy import run_path
from sys import argv

from sympy import GF, ZZ_I
from sympy.polys.matrices import DomainMatrix


@lru_cache(None)
def transverse_coefficients(center, degree, derivative_order, coordinate, shifted, first_approximant):
    """Coefficient of v^b in the q-scaled Gaussian-integer product in T+v."""
    p, q = first_approximant if coordinate == 1 else (355, 113)
    poly = [ZZ_I.one]
    for index in range(degree):
        origin = index % 2 if coordinate == 2 and shifted else 0
        constant = ZZ_I(0, 2 * p * (center - origin))
        poly = [x + y for x, y in zip([a * constant for a in poly] + [ZZ_I.zero],
                                      [ZZ_I.zero] + poly)]
    return tuple(a * (comb(d, derivative_order) * q**(d - derivative_order))
                 for d, a in enumerate(poly) if d >= derivative_order)


@lru_cache(None)
def log_moment(s, h, order):
    """s! [t^s] (1+t)^h log(1+t)^order equals the order-th derivative of (h)_s."""
    if order < 0 or order > s:
        return 0
    if s == 0:
        return 1
    return ((h - s + 1) * log_moment(s - 1, h, order)
            + order * log_moment(s - 1, h, order - 1))


def gaussian_entry(row, col, shifted, first_approximant=(22, 7)):
    j, s, b1, b2 = row
    h, a1, a2 = col
    first = transverse_coefficients(j, a1, b1, 1, shifted, first_approximant)
    second = transverse_coefficients(j, a2, b2, 2, shifted, first_approximant)
    return sum((x * y * log_moment(s, h, i + k) for i, x in enumerate(first[:s+1])
                for k, y in enumerate(second[:s+1-i])), ZZ_I.zero)


def exponent_in(n, prime):
    assert n > 0
    value = 0
    while n % prime == 0:
        n //= prime
        value += 1
    return value


def reduced_data(rows, cols, shifted):
    data = [[gaussian_entry(row, col, shifted) for col in cols] for row in rows]
    determinant = DomainMatrix(data, (len(rows), len(rows)), ZZ_I).det()
    scale = (prod(factorial(row[1]) for row in rows)
             * 7**(sum(col[1] for col in cols) - sum(row[2] for row in rows))
             * 113**(sum(col[2] for col in cols) - sum(row[3] for row in rows)))
    x, y = int(determinant.x), int(determinant.y)
    assert x or y
    common = gcd(scale, gcd(x, y))
    numerator_norm_squared = (x // common)**2 + (y // common)**2
    return numerator_norm_squared, scale // common


def diagonal_residue(rows, cols, shifted, difference_degree=None):
    # Move only the first coordinate to 2i*355/113, so all Newton column
    # coefficients keep exactly the same values used in the selected minor.
    field = GF(65537)
    max_transverse = max(b1 + b2 for _, _, b1, b2 in rows)

    def project(row, col):
        if difference_degree is not None and col == (0, 0, 1):
            # Divide the new column by (r1-r2)^(L-B) before specializing.
            j, s, b1, b2 = row
            if s or b1 + b2 != max_transverse:
                return field(0)
            degree = difference_degree
            coefficient = comb(degree, b1) * comb(degree - b1, b2) * (-1)**b2
            scaled = 113**(2 * degree - b1 - b2) * j**(degree - max_transverse)
            return field(coefficient * scaled)
        value = gaussian_entry(row, col, shifted, first_approximant=(355, 113))
        return field(int(value.x) + 256 * int(value.y))
    data = [[project(row, col) for col in cols] for row in rows]
    return int(DomainMatrix(data, (len(rows), len(rows)), field).det())


def floor_log10_rational(p, q):
    if p < q:
        d = floor_log10_rational(q, p)
        return -d if q == p * 10**d else -d-1
    quotient = p // q
    d = 0
    while quotient >= 10:
        quotient //= 10
        d += 1
    return d


def experiment(height, verifier):
    rows = verifier['row_packet'](height)
    selected = verifier['greedy_basis'](rows, verifier['columns'](height))
    raw_sq, raw_d = reduced_data(rows, selected, False)
    new_sq, new_d = reduced_data(rows, selected, True)
    # The floor of ten times log10 of |raw determinant| / |Newton determinant|
    # is the floor of log10((raw_sq * new_d^2 / (new_sq * raw_d^2))^5).
    scaled_log = floor_log10_rational((raw_sq * new_d**2)**5, (new_sq * raw_d**2)**5)
    diagonal = [diagonal_residue(rows, selected, shifted) for shifted in (False, True)]
    expected_diagonal = {14: [65398, 38025], 18: [40417, 53672]}
    assert diagonal == expected_diagonal[height]
    quotient_residue = diagonal_residue(rows, selected, True, difference_degree=height // 3)
    assert quotient_residue == {14: 50958, 18: 57195}[height]
    print(f'H={height}; exchanged determinant after removing its difference factor: '
          f'diagonal residue={quotient_residue} mod 65537')
    print(f'H={height}; scaled determinants at equal centers modulo 65537: raw/Newton={diagonal}')
    result = (len(str(raw_d)), len(str(new_d)), len(str(raw_sq)), len(str(new_sq)),
              exponent_in(raw_d, 113), exponent_in(new_d, 113), scaled_log)
    print(f'H={height}; raw/Newton denominator digits={result[:2]}; '
          f'reduced numerator squared-norm digits={result[2:4]}; '
          f'113-exponents={result[4:6]}; '
          f'floor(10*log10(abs(raw)/abs(Newton)))={scaled_log}')
    return result


if __name__ == '__main__':
    verifier_file = Path(argv[1]) if len(argv) > 1 else Path(__file__).with_name('verify_newton.py')
    verifier = run_path(str(verifier_file))
    assert experiment(14, verifier) == (114, 100, 316, 282, 53, 47, 34)
    assert experiment(18, verifier) == (264, 216, 731, 618, 120, 98, 89)
