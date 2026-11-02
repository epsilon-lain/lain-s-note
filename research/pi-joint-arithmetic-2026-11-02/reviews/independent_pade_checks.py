"""Exact independent checks of Padé coefficients and arithmetic bounds.

Does not alter the supplied archive. General asymptotics are audited in
INDEPENDENT_HEIGHT_PADE_REVIEW.md rather than inferred from these tests.
"""
import json
import math
from fractions import Fraction as Q
from pathlib import Path


def main():
    series_checks = 0
    denominator_checks = 0
    coefficient_checks = 0
    for m in range(1, 32):
        a = [Q(math.factorial(2*m-k)*math.factorial(m),
               math.factorial(2*m)*math.factorial(k)*math.factorial(m-k))
             for k in range(m+1)]
        for degree in range(2*m+11):
            lhs = sum((-1)**k*a[k]/math.factorial(degree-k)
                      for k in range(min(m, degree)+1))
            if degree <= m:
                lhs -= a[degree]
            ell = degree-2*m-1
            rhs = Q(0) if ell < 0 else Q(
                (-1)**m*math.factorial(m+ell)*math.factorial(m),
                math.factorial(2*m)*math.factorial(2*m+ell+1)
                *math.factorial(ell))
            assert lhs == rhs
            series_checks += 1
        for k in range(2*m+3):
            actual = a[k] if k <= m else Q(0)
            limit = Q(1, 2**k*math.factorial(k))
            assert 0 <= limit-actual <= limit*Q(k*(k-1), 2*m)
            coefficient_checks += 1
        if m < 3 or m % 2 == 0:
            continue
        coefficients = []
        for k in range(1,m+1,2):
            numerator = math.factorial(2*m-k)
            denominator = math.factorial(k)*math.factorial(m-k)
            assert numerator % denominator == 0
            coefficients.append(numerator//denominator)
        assert coefficients[-1] == 1
        assert coefficients[0] == math.factorial(2*m-1)//math.factorial(m-1)
        assert all(x > y for x,y in zip(coefficients, coefficients[1:]))
        for b in range(2,14):
            for a_num in range(1,2*b+2):
                if math.gcd(a_num,b) != 1:
                    continue
                square = Q(-a_num*a_num,b*b)
                value = sum(c*square**j for j,c in enumerate(coefficients))
                assert value and value.denominator == b**(m-1)
                denominator_checks += 1
    # Rational bounds establish the stated conservative m=10001 threshold:
    # cosh(1/2)-1 < (1/8)/(1-1/48)=6/47, e^4<3^4=81.
    threshold = 10001
    assert Q(6,47)+Q(12*81,threshold) < Q(1,4)
    assert Q(16*81,threshold) < Q(41,47)
    assert threshold+1 > 16
    result = {
        'status':'passed',
        'exact_remainder_series_coefficients':series_checks,
        'coefficient_error_bounds':coefficient_checks,
        'exact_full_denominator_checks':denominator_checks,
        'm_range':[1,31],
        'explicit_local_threshold_verified':threshold,
        'scope':'Exact finite algebra plus sufficient rational local estimates; no claim that finite tests prove the all-m remainder identity or asymptotic formulas.'
    }
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
