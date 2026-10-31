"""Exact finite checks for the fixed pi-word in complex bases.

Python standard library only. No stored digits of pi or continued fractions.
The natural-boundary and non-D-finite statements require the mathematical
proof and the Polya-Carlson theorem; they are NOT finite-test conclusions.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import hashlib
import json


def atan_bounds(inv: int, terms: int):
    s = sum((Q((-1)**j, (2*j+1)*inv**(2*j+1)) for j in range(terms)), Q(0))
    err = Q((-1)**terms, (2*terms+1)*inv**(2*terms+1))
    return min(s, s+err), max(s, s+err)


def pi_digits(count: int = 128, bits: int = 512):
    a, b = atan_bounds(5, 160)
    c, d = atan_bounds(239, 50)
    low, high = 16*a-4*d, 16*b-4*c
    scale = 1 << bits
    p = (low*scale).__floor__()
    P = -((-high*scale).__floor__())
    r, R = scale, scale
    digits = []
    for j in range(1, count+1):
        s = p*r//scale
        S = (P*R+scale-1)//scale
        digit = s//scale
        if digit != S//scale:
            raise ArithmeticError(f"Digit {j} ambiguous: increase precision")
        r, R = s-digit*scale, S-digit*scale
        assert 0 < r <= R < scale
        assert 0 <= digit <= 3
        digits.append(digit)
    return Q(p, scale), Q(P, scale), digits


@dataclass(frozen=True)
class C:
    re: Q = Q(0)
    im: Q = Q(0)

    def __add__(self, other):
        if not isinstance(other, C):
            other = C(Q(other))
        return C(self.re+other.re, self.im+other.im)
    __radd__ = __add__

    def __neg__(self):
        return C(-self.re, -self.im)

    def __sub__(self, other):
        return self + (-other if isinstance(other, C) else -Q(other))

    def __mul__(self, other):
        if not isinstance(other, C):
            other = C(Q(other))
        return C(self.re*other.re-self.im*other.im,
                 self.re*other.im+self.im*other.re)
    __rmul__ = __mul__

    def norm2(self):
        return self.re**2+self.im**2

    def inv(self):
        v = self.norm2()
        if not v:
            raise ZeroDivisionError
        return C(self.re/v, -self.im/v)

    def __pow__(self, k):
        if k < 0:
            return self.inv()**(-k)
        out, x = C(Q(1)), self
        while k:
            if k & 1:
                out = out*x
            x, k = x*x, k//2
        return out


def value(digits, z):
    u, acc = z.inv(), C()
    for digit in reversed(digits):
        acc = (acc+digit)*u
    return acc


def derivative(digits, z):
    u, power, out = z.inv(), z.inv(), C()
    for n, digit in enumerate(digits, 1):
        power = power*u
        out = out + (-n*digit)*power
    return out


def horner_derivative(digits, z):
    # Simultaneous Horner differentiation of g(u), then F'(z)=-u^2 g'(u).
    u, acc, grad = z.inv(), C(), C()
    for coefficient in list(reversed(digits))+[0]:
        grad, acc = grad*u+acc, acc*u+coefficient
    return -(u*u*grad)


def bound_strings(x, error, places=25):
    scale = 10**places
    a, b = ((x-error)*scale).__floor__(), -((-(x+error)*scale).__floor__())
    def dec(v):
        return ('-' if v < 0 else '')+str(abs(v)//scale)+'.'+str(abs(v)%scale).zfill(places)
    return [dec(a), dec(b)]


def main():
    lo, hi, digits = pi_digits()
    assert digits[:6] == [3,0,1,1,0,2]
    assert hi*hi < 10  # Hence ln(10/pi)>ln(pi): outside the log-time Taylor disk.
    i = C(Q(0), Q(1))
    bases = [C(Q(3)), C(Q(4)), C(Q(10)), C(Q(0),Q(3)),
             C(Q(3),Q(4)), C(Q(-3),Q(4)), C(Q(22,7),Q(1,100)),
             C(Q(0),Q(2)), C(Q(-1),Q(1))]
    for z in bases:
        assert derivative(digits,z) == horner_derivative(digits,z)
        # Angular flow: partial_theta F = i partial_ln_R F = i z F'.
        weighted = sum(((-n*digit)*z**(-n) for n,digit in enumerate(digits,1)), C())
        assert weighted*i == i*z*derivative(digits,z)
    # Exact four-direction Fourier filter selects positions modulo four.
    filters = 0
    for radius in (2,3,4):
        vals = [value(digits, C(Q(radius))*i**j) for j in range(4)]
        for residue in range(4):
            filtered = sum((vals[j]*i**(j*residue) for j in range(4)), C())*Q(1,4)
            direct = sum((Q(digit, radius**n) for n,digit in enumerate(digits,1)
                          if n % 4 == residue), Q())
            assert filtered == C(direct)
            filters += 1
    # Four-digit packing into ordinary integer base b^4.
    # It proves in the note: both rational coordinates would force eventual
    # periodicity of the entire original word. Check finite packing exactly.
    packed_checks = 0
    for b in (2,3,4,10):
        B = b**4
        real_packed = imag_packed = Q()
        for k in range(len(digits)//4):
            x1,x2,x3,x4 = digits[4*k:4*k+4]
            cr = x4-b*b*x2+3*b*b
            ci = x3-b*b*x1+3*b*b
            assert 0 <= cr < B and 0 <= ci < B
            real_packed += Q(cr,B**(k+1))
            imag_packed += Q(ci,B**(k+1))
            # Each packed coefficient uniquely recovers its two input digits.
            assert cr % (b*b) == x4 and 3-cr//(b*b) == x2
            assert ci % (b*b) == x3 and 3-ci//(b*b) == x1
            packed_checks += 2
        offset = sum((Q(3*b*b,B**(k+1)) for k in range(len(digits)//4)),Q())
        v = value(digits,C(Q(0),Q(b)))
        assert real_packed == v.re+offset
        assert imag_packed == v.im/b+offset
    # Universal analytic bound: tail derivative on |u|<=1/3 is at most 7/4.
    R = Q(1,3)
    tail_derivative = 3*(1/(1-R)**2-1-2*R)
    assert tail_derivative == Q(7,4)
    assert 3-tail_derivative == Q(5,4)
    assert 3+tail_derivative == Q(19,4)
    # Finite-prefix sanity checks; the all-series theorem is proved separately.
    outside = [z for z in bases if z.norm2() >= 9]
    pairs = 0
    for z,w in combinations(outside,2):
        delta = value(digits,z)-value(digits,w)
        relative = (z.inv()-w.inv()).norm2()
        assert Q(25,16)*relative <= delta.norm2() <= Q(361,16)*relative
        pairs += 1
    assert i**4 == C(Q(1))
    assert i**2+C(Q(1)) == C()
    beta = 2*i
    assert beta**2+3 == C(Q(-1))        # (103)_(2i) = -1.
    assert beta+2*beta.inv() == i       # (10.2)_(2i) = i.
    penney = C(Q(-1), Q(1))
    assert penney**2+2*penney+2 == C()
    result = {
        'status':'all_exact_finite_checks_passed',
        'scope':'Finite arithmetic checks, not proofs of natural boundary, non-D-finiteness, or P2848.',
        'native_digits':len(digits), 'precision_bits':512,
        'pi_method':'Machin identity; alternating rational remainder; outward dyadic propagation',
        'pi_dyadic_endpoints':[str(lo),str(hi)],
        'digit_prefix':''.join(map(str,digits[:64])),
        'digit_sha256':hashlib.sha256(bytes(digits)).hexdigest(),
        'derivative_checks':len(bases), 'angular_flow_checks':len(bases),
        'four_direction_filter_checks':filters,
        'packed_digit_reconstruction_checks':packed_checks,
        'finite_prefix_univalence_pair_checks':pairs,
        'univalence_constants':{'lower':'5/4','upper':'19/4','domain':'|z|,|w| >= 3; reciprocal-base distance'},
        'pi_squared_less_than_10_certified':True,
        'i_and_2i_finite_identities_checked':True,
        'certified_values':[]
    }
    for z, radius in [(C(Q(0),Q(2)),Q(2)),(C(Q(0),Q(3)),Q(3)),
                      (C(Q(3),Q(4)),Q(5)),(penney,Q(7,5))]:
        assert z.norm2() >= radius**2
        v = value(digits,z)
        error = Q(3)/(radius**len(digits)*(radius-1))
        result['certified_values'].append({'base':[str(z.re),str(z.im)],
          'real_interval':bound_strings(v.re,error), 'imag_interval':bound_strings(v.im,error),
          'complex_absolute_tail_bound':str(error)})
    path = Path(__file__).with_name('complex_analyticity_results.json')
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
