"""Exact coefficient certificate for S(Psi)(t) < 0 for every t >= 3.

The infinite result uses the monotone positive tail bound in the companion
proof. This program constructs its finite rational certificate, not a grid test.
"""
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('transport', HERE/'verify_transport.py')
transport = importlib.util.module_from_spec(spec)
spec.loader.exec_module(transport)


def record(x):
    return {'numerator': x.numerator, 'denominator': x.denominator}


def main():
    _, _, digits = transport.certify_digits(11)
    assert digits == [3, 0, 1, 1, 0, 2, 1, 1, 1, 0, 0]
    coefficients = {k: digits[k] for k in range(2, 11) if digits[k]}
    polynomial = {}
    for i, a in coefficients.items():
        for j, b in coefficients.items():
            power = i+j
            polynomial[power] = polynomial.get(power, 0) + (
                i*j*(i*i+j*j-3*i*j+1)*a*b)
    polynomial = {n:c for n,c in sorted(polynomial.items()) if c}
    assert polynomial[4] == -12
    assert {n:c for n,c in polynomial.items() if c > 0} == {9:48}

    r = F(1, 3)
    S1 = sum((k*a*r**k for k,a in coefficients.items()), F(0))
    S3 = sum((k**3*a*r**k for k,a in coefficients.items()), F(0))
    T1 = 3*(r/(1-r)**2 - sum((k*r**k for k in range(1,11)), F(0)))
    T3 = 3*(r*(1+4*r+r*r)/(1-r)**4
             - sum((k**3*r**k for k in range(1,11)), F(0)))
    positive_prefix = sum((max(c,0)*r**(n-4)
                           for n,c in polynomial.items()), F(0))
    tail_bound = (2*((S1+T1)*(S3+T3)-S1*S3)
                  + (S1+T1)**2-S1*S1)/r**4
    combined = positive_prefix+tail_bound
    assert positive_prefix == F(16,81)
    assert combined == F(9079787,3188646)
    assert combined < 3
    schwarzian_constant = (1-r)**4/(2*(2-r)**2)
    assert schwarzian_constant == F(8,225)
    result = {
        'native_prefix_certified': ''.join(map(str,digits)),
        'active_coefficients_through_10': coefficients,
        'N0_polynomial_coefficients': polynomial,
        'evaluation_point_r': record(r),
        'prefix_moment_1': record(S1),
        'prefix_moment_3': record(S3),
        'tail_bound_moment_1': record(T1),
        'tail_bound_moment_3': record(T3),
        'positive_prefix_bound': record(positive_prefix),
        'positive_tail_bound': record(tail_bound),
        'combined_positive_bound': record(combined),
        'combined_strictly_less_than_3': True,
        'N_over_r4_strictly_less_than': -9,
        'uniform_Schwarzian_upper_bound': '-8/(225*t^2), t>=3',
        'scope': ('Exact finite coefficient and rational tail-bound certificate. '
                  'The companion proof explains uniformity on t>=3. '
                  'No P2848 conclusion and no Lean verification.')
    }
    (HERE/'schwarzian_ode_results.json').write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'prefix':result['native_prefix_certified'],
                      'bound':str(combined), 'verified':True}))


if __name__ == '__main__':
    main()
