"""Exact finite checks for the normalized endpoint norm theorem.

Python standard library only. No numerical pi, matrix square roots or singular
values are used; the all-size norm theorem is proved in the companion note.
"""
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import json


ZERO = (F(0), F(0))
ONE = (F(1), F(0))


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def scale(a, c):
    return a[0]*c, a[1]*c


def multiply(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def inverse(a):
    norm = a[0]*a[0]+a[1]*a[1]
    assert norm
    return a[0]/norm, -a[1]/norm


def power(a, k):
    if k < 0:
        return power(inverse(a), -k)
    out = ONE
    for _ in range(k):
        out = multiply(out, a)
    return out


def nu_direct(k, z):
    out = ZERO
    for ell in range(k+1):
        out = add(out, scale(power(z, -ell-1),
                             F((-1)**ell*factorial(k), factorial(k-ell))))
    return out


def exp_polynomial(k, z):
    out = ZERO
    for j in range(k+1):
        out = add(out, scale(power(z, j), F(1, factorial(j))))
    return out


def multiply_matrix(a, b):
    return [[sum(x*y for x, y in zip(row, col))
             for col in zip(*b)] for row in a]


def hilbert_inverse_formula(n):
    return [[F((-1)**(i+j)*(i+j+1)*comb(n+i, n-j-1)*
               comb(n+j, n-i-1)*comb(i+j, i)**2)
             for j in range(n)] for i in range(n)]


def gram_checks():
    inverse_checks = legendre_checks = diagonal_checks = 0
    for n in range(2, 15):
        gram = [[F(1, i+j+1) for j in range(n)] for i in range(n)]
        inv = hilbert_inverse_formula(n)
        assert multiply_matrix(gram, inv) == [
            [F(i == j) for j in range(n)] for i in range(n)]
        inverse_checks += 1
        for i in range(n):
            expected = (2*i+1)*comb(n+i, n-i-1)**2*comb(2*i, i)**2
            assert inv[i][i] == expected <= 2*n*2**(8*n)
            diagonal_checks += 1
        # P_j(2u-1) coefficients in the monomial basis.
        coeff = [[F((-1)**(j-k)*comb(j, k)*comb(j+k, k))
                  if k <= j else F(0) for k in range(n)]
                 for j in range(n)]
        transformed = multiply_matrix(multiply_matrix(coeff, gram),
                                      list(map(list, zip(*coeff))))
        expected = [[F(1, 2*i+1) if i == j else F(0)
                     for j in range(n)] for i in range(n)]
        assert transformed == expected
        legendre_checks += 1
    return {'hilbert_inverse_matrix_checks': inverse_checks,
            'inverse_diagonal_checks': diagonal_checks,
            'shifted_legendre_gram_checks': legendre_checks}


def moment_checks():
    count = 0
    centers = [(F(1), F(2)), (F(-3, 2), F(1, 3)),
               (F(0), F(6)), (F(0), F(44, 7))]
    for z in centers:
        previous = None
        for k in range(33):
            direct = nu_direct(k, z)
            closed = scale(multiply(power(z, -k-1),
                                    exp_polynomial(k, scale(z, -1))),
                           (-1)**k*factorial(k))
            assert direct == closed
            if k:
                assert direct == add(inverse(z), scale(multiply(inverse(z), previous), -k))
            else:
                assert direct == inverse(z)
            previous = direct
            count += 1
    return count


def tail_bound(k):
    T = F(44, 7)
    assert k+2 > T
    return T**(k+1)/factorial(k+1)/(1-T/F(k+2))


def tail_checks():
    T = F(44, 7)
    base = tail_bound(16)
    assert base < F(1, 6) < F(1, 2)
    count = 0
    for k in range(16, 101):
        ratio = T/F(k+2)*(1-T/F(k+2))/(1-T/F(k+3))
        assert tail_bound(k+1)/tail_bound(k) == ratio < 1
        assert tail_bound(k) < F(1, 6)
        # The partial positive series tail must remain below its majorant.
        partial = sum(T**j/F(factorial(j)) for j in range(k+1, k+20))
        assert partial < tail_bound(k)
        count += 1
    records = []
    for n in (9, 10, 12, 16, 24, 32, 48, 64):
        rational_lower = F(factorial(2*n-1), 2)/T**(2*n-1)
        old_upper = 2*3**8*n**3*2**(8*n)*factorial(2*n-2)
        assert 0 < rational_lower < old_upper
        assert (2*n-1)*F(factorial(2*n-2), 2)/T**(2*n-1) == rational_lower
        records.append({'n':n,
                        'rational_lower_numerator_bits':rational_lower.numerator.bit_length(),
                        'rational_lower_denominator_bits':rational_lower.denominator.bit_length(),
                        'upper_bound_bits':old_upper.bit_length()})
    return {'exact_tail_base':str(base), 'tail_base_strict_upper':'1/6',
            'tail_monotonicity_and_majorant_checks':count,
            'rational_bound_records':records}


def main():
    result = {'status':'passed', **gram_checks(),
              'gaussian_rational_moment_identity_checks':moment_checks(),
              **tail_checks(),
              'scope':'Exact finite premises and identities supplement the proof of '
                      'log norm=2n log n+O(n). No numerical pi or matrix square root. '
                      'No determinant, optimal disk-radius or derivative sharpness claim.'}
    Path(__file__).with_name('gram_norm_sharpness_verification.json').write_text(
        json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
