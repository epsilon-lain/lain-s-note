"""Exact examples for local normalizers and height-aware residual certificates.

The polynomial maps are evaluated directly over Fraction; the separate checker
uses homogeneous integer forms. No claimed integer-definition theorem is used.
"""
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent


def valuation(value, prime):
    if not isinstance(prime,int) or prime < 2:
        raise ValueError('prime must be an integer at least two')
    value = F(value)
    if not value:
        raise ValueError('valuation at zero is not a finite integer')
    def count(n):
        n, exponent = abs(n), 0
        while n % prime == 0:
            n //= prime
            exponent += 1
        return exponent
    return count(value.numerator) - count(value.denominator)


def normalizer(x, order=1):
    x = F(x)
    if not isinstance(order,int) or order < 1:
        raise ValueError('order must be a positive integer')
    beta = 2504 + 51*x*x
    r4 = beta * (beta-1) * (beta-2401) * (32*beta+343)
    result = 149**(2*order+1) * r4 / (1+149**2*r4)
    p, q = x.numerator, x.denominator
    contact_numerators = [q*q*y for y in [beta,beta-1,beta-2401,32*beta+343]]
    assert all(y.denominator == 1 for y in contact_numerators)
    product = 1
    for y in contact_numerators:
        product *= int(y)
    raw_numerator = 149**(2*order+1)*product
    raw_denominator = q**8+149**2*product
    actual_gcd = gcd(raw_numerator, raw_denominator)
    expected_gcd = (2**min(5, 2*valuation(q,2)) *
                    3**(4*(q % 3 == 0)) * 17**(4*(q % 17 == 0)) *
                    149**(2*(q % 149 == 0)))
    assert actual_gcd == expected_gcd
    assert result == F(raw_numerator, raw_denominator)
    expected_v149 = 2*order-1 if q % 149 == 0 else 2*order+1
    assert valuation(result,149) == expected_v149
    return {'x': str(x), 'order': order, 'beta': str(beta), 'R': str(result),
            'input_denominator': str(q), 'output_denominator': str(result.denominator),
            'raw_gcd': str(actual_gcd), 'expected_gcd': str(expected_gcd),
            'v149_R': expected_v149}


def quartic_example(x, b, s, t):
    coords = tuple(map(F,(x,b,s,t)))
    x,b,s,t = coords
    residuals = [x-b,s*b,s*t-1]
    L = lcm(*(z.denominator for z in coords))
    cleared = [L*L*r for r in residuals]
    assert all(r.denominator == 1 for r in cleared)
    loss = sum(r*r for r in residuals)
    normalized = L**4*loss
    assert normalized == sum(r*r for r in cleared)
    assert normalized.denominator == 1 and (not loss or normalized >= 1)
    return {'coordinates': list(map(str,coords)), 'L': str(L),
            'loss': str(loss), 'normalized_loss': str(normalized)}


def atan_interval(t, bits=192):
    # Exact alternating partial sums; no floating-point pi or sqrt input.
    partial, j = F(0), 0
    threshold = F(1,2**bits)
    while True:
        term = F((-1)**j,(2*j+1)*t**(2*j+1))
        if abs(term) < threshold:
            return tuple(sorted((partial,partial+term)))
        partial += term
        j += 1


def pi_interval():
    z = (1,0)
    for _ in range(4):
        z = (5*z[0]-z[1],z[0]+5*z[1])
    assert z == (476,480) == (2*239-2,2*239+2)
    l5,u5 = atan_interval(5)
    l239,u239 = atan_interval(239)
    assert 0 < 4*l5-u239 < 4*u5-l239 < 1
    lower,upper = 16*l5-4*u239,16*u5-4*l239
    assert 3 < lower < upper < 4
    return lower,upper


def linear_range_certificate(lower, upper, forms):
    lower,upper = F(lower),F(upper)
    assert lower < upper
    if len(forms) != 2:
        raise ValueError('two independent forms are required')
    (A0,B0),(A1,B1) = forms
    assert all(isinstance(v,int) for v in [A0,B0,A1,B1])
    if A0 == 0 or A1 == 0:
        raise ValueError('both useful linear forms require nonzero leading coefficients')
    assert A0*B1-A1*B0 != 0
    errors = [max(abs(A*lower-B),abs(A*upper-B)) for A,B in forms]
    largest = max(errors)
    cap = (1/(2*largest)).numerator // (1/(2*largest)).denominator
    if cap < 1:
        raise ValueError('the supplied error bounds certify no positive denominator')
    assert cap*largest <= F(1,2) < (cap+1)*largest
    coefficient = 2*max(abs(A0),abs(A1))
    return {'forms': forms, 'eta_upper': list(map(str,errors)),
            'all_integer_numerators': True, 'denominator_min': 1,
            'denominator_max': cap, 'error_lower_bound': f'1/({coefficient}*q)',
            'scope': 'Finite denominator range; known 22/7 and 355/113 are inputs, not discoveries.'}


def run():
    normalizers = [normalizer(x,N) for N in [1,2]
                   for x in [F(0),F(3),F(22,7),F(355,113),F(1,2),
                             F(1,149),F(3,149**2),F(149)]]
    shrinking = [quartic_example(F(1,2),F(1,2),F(1,M),M)
                 for M in [1,2,4,8,16,64]]
    sharp = [quartic_example(F(1,L),F(1,L),F(1,L),L)
             for L in [2,3,5,7,16,149]]
    assert all(F(row['loss']) == F(1,4*int(row['coordinates'][-1])**2)
               for row in shrinking)
    assert all(row['normalized_loss'] == '1' for row in sharp)
    lower,upper = pi_interval()
    finite_pi = linear_range_certificate(lower,upper,[(7,22),(113,355)])
    assert finite_pi['denominator_max'] == 56
    out = {'checks': 'PASS', 'normalizer_examples': normalizers,
           'fixed_false_input_residuals': shrinking,
           'sharp_quartic_gap_examples': sharp,
           'pi_interval': list(map(str,(lower,upper))),
           'finite_pi_linear_certificate': finite_pi,
           'scope': 'Elementary constructive certificates and examples, not a new pi approximation or proof of the source paper.'}
    (ROOT/'examples.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'checks':'PASS','normalizer_examples':len(normalizers),
                      'residual_examples':len(shrinking)+len(sharp),
                      'pi_denominator_range':[1,finite_pi['denominator_max']],
                      'pi_error_lower_bound':finite_pi['error_lower_bound']}))


if __name__ == '__main__':
    run()
