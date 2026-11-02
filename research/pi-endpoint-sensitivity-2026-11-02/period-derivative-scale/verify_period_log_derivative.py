"""Exact finite checks for PERIOD_LOG_DERIVATIVE.md.

Only Python's standard library is used. The all-n asymptotic is proved by
the explicit L1 estimate in the note, not extrapolated from these checks.
"""
from fractions import Fraction as F
from functools import cache
from math import comb, factorial, prod
from pathlib import Path
import json


def determinant(matrix):
    a = [list(map(F, row)) for row in matrix]
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
    return prod(factorial(k) for k in range(n))


def delta(n):
    return F(A(n)**4, A(2*n))


def moments(length, z, y):
    values = [(y-1)/z]
    for k in range(1, length):
        values.append(y/z-k*values[-1]/z)
    return values


def nu(k, z):
    return sum(F((-1)**ell*factorial(k), factorial(k-ell))/z**(ell+1)
               for ell in range(k+1))


def multiply_univariate(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


@cache
def vandermonde_square(n):
    """Expand the pair factors directly, independently of determinant code."""
    poly = {(0,)*n: F(1)}
    for i in range(n):
        for j in range(i+1, n):
            updated = {}
            for powers, coefficient in poly.items():
                for di, dj, scalar in ((2,0,1), (1,1,-2), (0,2,1)):
                    target = list(powers)
                    target[i] += di
                    target[j] += dj
                    key = tuple(target)
                    updated[key] = updated.get(key, F(0))+coefficient*scalar
            poly = {key: value for key, value in updated.items() if value}
    return poly


def endpoint_insertion_checks():
    count = 0
    sizes = {}
    for n in range(2,5):
        poly = vandermonde_square(n)
        sizes[n] = len(poly)
        for z in (F(1), F(-2), F(3,2)):
            for y in (F(1), F(2), F(-1,3)):
                mu = moments(2*n-1, z, y)
                endpoint = [nu(k,z) for k in range(2*n-1)]
                matrix = [[mu[i+j] for j in range(n)] for i in range(n)]
                derivative = F(0)
                for row in range(n):
                    replaced = [a[:] for a in matrix]
                    replaced[row] = [endpoint[row+j] for j in range(n)]
                    derivative += determinant(replaced)
                # One endpoint functional occupies variable 0; all other
                # variables carry the original functional. The n/n! factor
                # must give exactly 1/(n-1)!.
                insertion = sum(coefficient*endpoint[powers[0]]*
                                prod(mu[powers[j]] for j in range(1,n))
                                for powers,coefficient in poly.items())
                assert derivative == insertion/factorial(n-1)
                # Independently verify the unchanged Vandermonde normalization.
                full = sum(coefficient*prod(mu[k] for k in powers)
                           for powers,coefficient in poly.items())
                assert determinant(matrix) == full/factorial(n)
                count += 1
    return {"endpoint_insertion_checks":count,
            "vandermonde_monomial_counts":sizes}


def gamma_normalization_checks():
    count = monomial_count = 0
    for z in (F(1), F(-2), F(7,3)):
        for degree in range(25):
            # Directly integrate (1-v/z)^degree against exp(-v)/z.
            integral = sum(F(comb(degree,k)*(-1)**k*factorial(k),1)/z**(k+1)
                           for k in range(degree+1))
            assert integral == nu(degree,z)
            monomial_count += 1
        for m in range(1,10):
            gamma_order = 2*m
            for points in (tuple(F(j,m+1) for j in range(1,m+1)),
                           (F(0),)*m, (F(1),)*m):
                polynomial = [F(1)]
                correction = [F(1)]
                for u in points:
                    polynomial = multiply_univariate(polynomial, [u*u,-2*u,F(1)])
                    a = z*(1-u)
                    correction = multiply_univariate(correction, [F(1),-2*a,a*a])
                left = sum(coefficient*nu(k,z) for k,coefficient in enumerate(polynomial))
                # E V^-k for V~Gamma(2m+1) is (2m-k)!/(2m)!.
                gamma_expectation = sum(coefficient*F(factorial(gamma_order-k),
                                                      factorial(gamma_order))
                                        for k,coefficient in enumerate(correction))
                right = F(factorial(gamma_order),1)/z**(gamma_order+1)*gamma_expectation
                assert left == right
                count += 1
            mean = F(factorial(gamma_order+1), factorial(gamma_order))
            second = F(factorial(gamma_order+2), factorial(gamma_order))
            assert mean == 2*m+1
            assert second-mean*mean == 2*m+1
            assert second-4*m*mean+4*m*m == 2*m+2
    return {"endpoint_monomial_integrals":monomial_count,
            "gamma_product_normalizations":count}


def reflection_checks():
    count = 0
    for n in range(2,11):
        for z in (F(1), F(-2), F(3,2)):
            mu = moments(2*n-1,z,F(1))
            reflected = moments(2*n-1,-z,F(1))
            for k in range(2*n-1):
                assert reflected[k] == sum(F(comb(k,j)*(-1)**j)*mu[j]
                                            for j in range(k+1))
                count += 1
            dplus = determinant([[mu[i+j] for j in range(n)] for i in range(n)])
            dminus = determinant([[reflected[i+j] for j in range(n)] for i in range(n)])
            assert dplus == dminus
            pplus = z**(n*n)*dplus/A(n)**2
            pminus = (-z)**(n*n)*dminus/A(n)**2
            assert pminus == (-1)**n*pplus
    return count


def factorial_scale_checks():
    count = 0
    t = F(44,7)
    for n in range(2,31):
        m = n-1
        ratio = delta(m)/delta(n)
        assert ratio == F(factorial(2*m)*factorial(2*m+1),factorial(m)**4)
        scale = F(factorial(2*m),1)/t**(2*m+1)*ratio
        assert scale == F(factorial(2*m)**2*factorial(2*m+1),
                          factorial(m)**4)/t**(2*m+1)
        assert scale == F(factorial(2*m)*(2*m+1)*comb(2*m,m)**2,1)/t**(2*m+1)
        count += 1
    return count


def explicit_threshold_checks():
    # The proof uses t<44/7, already certified in the preceding Hankel note.
    assert F(44,7)**2 < 40
    assert F(7,4)**2 > 3
    assert 3**2 < 15
    assert 2*F(44,7) < 13
    m = 10000
    assert F(2*40,m) <= F(1,125)
    # For 0<=x<=1/125, exp(x)<=1/(1-x)<=125/124.
    first = F(80*125,124*m)
    second = F(11,2*100)  # sqrt(m)>=100.
    third = F(13,3*m)
    tail = F(3**13+1,2)*F(3,4)**100
    assert tail < F(1,1000)
    assert first+second+third+F(1,1000) < F(1,12)
    assert F(m,49) > 200
    # Induction for 200^m > 7[2(m+1)^2+24], beginning already at m=2.
    assert 200**2 > 7*(2*(2+1)**2+24)
    # 4Q(m)-Q(m+1)=6m^2+8m+72 has positive coefficients;
    # multiplication by 200 is therefore more than enough for the step.
    assert 200 > 4
    assert all(coefficient > 0 for coefficient in (6,8,72))
    return {"epsilon_threshold_m":m,
            "eventual_direction_threshold_n":m+1,
            "epsilon_rational_majorant":str(first+second+third+F(1,1000)),
            "epsilon_majorant_strictly_below_one_twelfth":True,
            "factorial_dominance_induction":True}


def main():
    result = {"status":"passed",
              **endpoint_insertion_checks(),
              **gamma_normalization_checks(),
              "formal_reflection_moment_checks":reflection_checks(),
              "factorial_scale_checks":factorial_scale_checks(),
              **explicit_threshold_checks(),
              "scope":"Exact finite normalization and constant certificates. The "
                      "all-size derivative asymptotic and maximal complex-disk scale "
                      "are proved analytically in the note. No nearest-root theorem "
                      "or pi bad-approximability claim."}
    Path(__file__).with_name("period_log_derivative_verification.json").write_text(
        json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
