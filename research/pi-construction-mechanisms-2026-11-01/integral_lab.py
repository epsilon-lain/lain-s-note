"""Exact small Zeilberger–Zudilin forms and numerical rate ablations.

Run with Python containing sympy and mpmath. No network and no source files
are required. Algebraic checks use exact rational arithmetic; decimal output
is evidence for inspection, not a replacement for the paper's asymptotic proof.
"""
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as sp

t, z = sp.symbols("t z")


def zz_form(n):
    order = 3*n+1
    numerator = 5*t**(2*n)*(t**4+6*t**2+25)**(2*n)
    denominator = (25-t**2)**order
    local = sp.Poly(sp.expand(numerator.subs(t, z-5)), z)
    coeffs = []
    for j in range(3*n+1):
        degree = 3*n-j
        coeffs.append(sp.Add(*[
            local.nth(k)*sp.binomial(order+degree-k-1, degree-k)
            / sp.Integer(10)**(order+degree-k)
            for k in range(degree+1)
        ]))
    polynomial, _ = sp.div(numerator, denominator, t)
    paired = sp.Add(*[
        coeffs[j]*((5+t)**(-j-1)+(5-t)**(-j-1))
        for j in range(3*n+1)
    ])
    assert sp.cancel(numerator/denominator-polynomial-paired) == 0
    primitive = sp.integrate(polynomial, t)
    for j in range(1, 3*n+1):
        primitive += coeffs[j]/j*((5-t)**(-j)-(5+t)**(-j))
    a, b = -1-2*sp.I, -1+2*sp.I
    u = sp.simplify(sp.I*(primitive.subs(t,b)-primitive.subs(t,a)))
    v = -coeffs[0]/2
    assert u.is_Rational and v.is_Rational and v != 0
    Pz = 2+6*z+9*z**2+6*z**3+2*z**4
    # A second, independent residue extraction after the Möbius substitution.
    # Coefficients of the remaining inverse power are positive binomials.
    top = sp.Poly(sp.expand((1+z)**(2*n)*Pz**(2*n)),z)
    if n:
        positive_coefficient = sp.Add(*[
            top.nth(k)*sp.binomial(4*n+3*n-k-1,3*n-k)
            for k in range(3*n+1)
        ])
    else:
        positive_coefficient = sp.Integer(1)
    assert positive_coefficient > 0
    assert v == (-1)**(n+1)*sp.Integer(4)**(n-1)*positive_coefficient
    phi = math.prod(p for p in sp.primerange(7,2*n+1)
                    if p*p > 3*n and sp.Rational(1,2) <= sp.Rational(n,p)%1 < sp.Rational(2,3))
    lcm = sp.ilcm(*range(1,4*n+1)) if n else sp.Integer(1)
    multiplier = sp.Integer(2)**(2-(5*n//2))*lcm/phi
    U, V = sp.cancel(multiplier*u), sp.cancel(multiplier*v)
    assert U.is_Integer and V.is_Integer
    return dict(n=n, pole_coefficients=list(map(str,coeffs)), polynomial=str(polynomial),
                u=str(u),v=str(v),positive_residue_coefficient=str(positive_coefficient),
                multiplier=str(multiplier),U=str(U),V=str(V),
                J_decimal=str(sp.N(u+v*sp.pi,35)))


def rate_ablations():
    mp.mp.dps = 60
    poly = lambda x: 108*x**3-2359989*x**2+138304*x-2048
    small = mp.findroot(poly, (mp.mpc('.029','.003'),mp.mpc('.029','.004')))
    large = mp.findroot(poly,(mp.mpf(21000),mp.mpf(22000)))
    r, s = mp.log(large), mp.log(abs(small))
    prime_saving = mp.pi/(2*mp.sqrt(3))-mp.log(3*mp.sqrt(3)/4)
    costs = {"naive_lcm":mp.mpf(4),
             "2_adic_saving":4-mp.mpf(5)/2*mp.log(2),
             "2_adic_and_removable_primes":4-mp.mpf(5)/2*mp.log(2)-prime_saving}
    out = {}
    for label,c in costs.items():
        sigma,tau=r+c,-s-c
        out[label] = dict(cost=str(c),sigma=str(sigma),tau=str(tau),
                          measure_bound=str(1+sigma/tau) if tau>0 else None)
    return dict(raw_coefficient_rate=str(r),raw_integral_rate=str(s),
                removable_prime_saving=str(prime_saving),ablations=out)


def main():
    # Exact identities establish the displayed Möbius substitution for every n.
    tz = -5*(1+z)/(1-z)
    Pz = 2+6*z+9*z**2+6*z**3+2*z**4
    Sz = (1+z)**2*Pz**2/(1-z)**4
    kernel = t**2*(t**4+6*t**2+25)**2/(25-t**2)**3
    assert sp.cancel(5/(25-tz**2)*sp.diff(tz,z)-1/(2*z)) == 0
    assert sp.cancel(kernel.subs(t,tz)+4*Sz/z**3) == 0
    forms = [zz_form(n) for n in range(5)]
    for first,second in zip(forms,forms[1:]):
        det = sp.Integer(first["U"])*sp.Integer(second["V"])-sp.Integer(second["U"])*sp.Integer(first["V"])
        assert det != 0
        first["adjacent_determinant"] = str(det)
    # Symmetry ablation is exact from elementary logarithmic primitives.
    # i∫5/(25−t²)dt = −π/4; i∫5t/(25−t²)dt = 5 arctan(1/7).
    result = dict(forms=forms, rates=rate_ablations(),
                  symmetry_ablation={"even":"-pi/4", "multiply_by_t":"5*atan(1/7)"})
    destination=Path(__file__).with_name("integral_lab_results.json")
    destination.write_text(json.dumps(result,indent=2)+"\n")
    print(destination)
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
