"""Exact finite checks for the Hilbert-normalized local-period disk argument.

Only Python's standard library is used. The infinite norm/disk statements
are proved in LOCAL_PERIOD_DISC.md, not extrapolated from these checks.
"""
from fractions import Fraction as F
from functools import cache
from math import comb, factorial, prod
from pathlib import Path
import hashlib
import json


def determinant(matrix):
    a = [list(map(F, row)) for row in matrix]
    if not a:
        return F(1)
    result = F(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            result = -result
        diagonal = a[i][i]
        result *= diagonal
        for j in range(i+1, len(a)):
            ratio = a[j][i]/diagonal
            for k in range(i+1, len(a)):
                a[j][k] -= ratio*a[i][k]
    return result


@cache
def A(n):
    return prod(factorial(j) for j in range(n))


def moments(n, z, y):
    values = [(y-1)/z]
    for k in range(1, 2*n-1):
        values.append(y/z-k*values[-1]/z)
    return values


def nu(k, z):
    return sum(F((-1)**ell*factorial(k), factorial(k-ell))/z**(ell+1)
               for ell in range(k+1))


def raw_general(n, z, y):
    rows = [(j,s) for j in (0,1) for s in range(n)]
    cols = [(h,a) for h in (0,1) for a in range(n)]
    return determinant([
        [y**(h*j)*sum(F(comb(a,k))*F(j*z)**(a-k)*
                      F(h**(s-k),factorial(s-k))
                      for k in range(min(s,a)+1))
         for h,a in cols] for j,s in rows])


def endpoint_identity_checks():
    count = derivative_count = 0
    for n in range(2,9):
        for z in (F(1), F(-2), F(1,2)):
            for y in (F(1), F(2), F(3,2), F(-1)):
                values = moments(n,z,y)
                hankel = determinant([[values[i+j] for j in range(n)] for i in range(n)])
                assert raw_general(n,z,y) == z**(n*n)*hankel/A(n)**2
                count += 1
            at_one = moments(n,z,F(1))
            at_two = moments(n,z,F(2))
            at_half = moments(n,z,F(1,2))
            for k in range(2*n-1):
                derivative = nu(k,z)
                assert at_two[k]-at_one[k] == derivative
                assert at_half[k] == at_one[k]-derivative/2
                if k:
                    assert derivative == 1/z-k*nu(k-1,z)/z
                else:
                    assert derivative == 1/z
                derivative_count += 1
    return {'raw_hankel_identity_checks':count,
            'endpoint_derivative_moment_checks':derivative_count}


def hilbert_checks():
    checks = 0
    for n in range(2,15):
        matrix = [[F(1,i+j+1) for j in range(n)] for i in range(n)]
        delta = determinant(matrix)
        assert delta == F(A(n)**4,A(2*n)) > 0
        B = 2*n*2**(8*n)
        for i in range(n):
            cofactor = determinant([[matrix[j][k] for k in range(n) if k != i]
                                   for j in range(n) if j != i])
            diagonal = (2*i+1)*comb(n+i,n-i-1)**2*comb(2*i,i)**2
            assert cofactor/delta == diagonal <= B
            checks += 1
    return checks


def L(n):
    return 2*3**8*n**3*2**(8*n)*factorial(2*n-2)


def local_constant_checks():
    moment_checks = 0
    for n in range(2,41):
        bound = 3**8*factorial(2*n-2)
        for k in range(2*n-1):
            direct_majorant = sum(F(factorial(k),factorial(k-ell)*5**(ell+1))
                                   for ell in range(k+1))
            assert direct_majorant <= bound
            moment_checks += 1
        ell = L(n)
        assert ell >= n
        r = F(1,8*n*ell)
        assert r <= F(1,8*n*n) <= F(1,32)
        assert 5 < 6-r and F(44,7)+r < 8
        exponent = F(n*n,6)*r+n*r+2*n*r*ell
        assert exponent <= F(19,48) < F(1,2)
        error_denominator = 32*3**8*n**4*2**(8*n)*factorial(2*n-2)
        assert r/2 == F(1,error_denominator)
    return {'rational_moment_majorants':moment_checks,
            'local_radius_constant_checks':39}


def rational_center_certificates():
    cases = []
    for b in (2,3,7,10):
        n = 24*b
        N = n*n//2
        d = A(n)**2//A(n//2)**4
        assert d*A(n//2)**4 == A(n)**2
        ratio = 2*b**N*F(44,7)**(n*n)*F(d*A(n)**2,A(2*n))
        assert ratio < 1
        error_denominator = 32*3**8*n**4*2**(8*n)*factorial(2*n-2)
        text = str(error_denominator)
        digits = len(text)
        assert error_denominator < 10**digits
        cases.append({'b':b, 'n':n, 'polynomial_degree':N,
                      'twice_rational_period_upper_below_value_floor':True,
                      'comparison_numerator_bits':ratio.numerator.bit_length(),
                      'comparison_denominator_bits':ratio.denominator.bit_length(),
                      'certified_error_lower_bound_formula':
                          '1/[32*3^8*n^4*2^(8n)*(2n-2)!]',
                      'strict_error_lower_bound_power_of_ten':-digits,
                      'error_denominator_bits':error_denominator.bit_length(),
                      'error_denominator_decimal_sha256':hashlib.sha256(text.encode()).hexdigest()})
    return cases


def main():
    result = {'status':'passed',
              **endpoint_identity_checks(),
              'hilbert_inverse_diagonal_checks':hilbert_checks(),
              **local_constant_checks(),
              'rational_center_certificates':rational_center_certificates(),
              'scope':'Finite exact identities, constants, and sufficient rational-center '
                      'exclusions. General Gram norm/disk theorem is proved in the note; '
                      'no claim about pi bad approximability or a best irrationality bound.'}
    Path(__file__).with_name('local_period_disc_verification.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
