"""Exact maximal-flatness coefficients, derivative kernel and shared anchors.
SymPy is used for finite polynomial identities, not floating-point tests.
"""
from fractions import Fraction as F
from math import comb,factorial,lcm
from pathlib import Path
import json
import sympy as s
w=s.symbols('w')

def coeffs(r):
    return [F(2*(-1)**(j+1)*factorial(r)**2,factorial(r-j)*factorial(r+j)) for j in range(1,r+1)]

def run():
    cases=[];moments=0;ranks=0
    for r in range(1,17):
        c=coeffs(r)
        assert sum(c)==1
        for k in range(1,r):
            assert sum(v*j**(2*k) for j,v in enumerate(c,1))==0;moments+=1
        C=F(4**r,comb(2*r,r))
        # w=e^(2*pi*i*y). sin(pi*y)^2=(2-w-w**(-1))/4.
        left=1-sum(s.Rational(v.numerator,v.denominator)*(w**j+w**(-j))/2 for j,v in enumerate(c,1))
        right=s.Rational(C.numerator,C.denominator)*((2-w-1/w)/4)**r
        assert s.expand((left-right)*w**r)==0
        # Exact leading term of H is C*pi^(2r)/(2r+1).
        signed_moment=sum(v*j**(2*r) for j,v in enumerate(c,1))
        assert signed_moment==(-1)**(r+1)*factorial(r)**2
        assert C<=2*r+1
        D=lcm(*(v.denominator*j for j,v in enumerate(c,1)))
        bound=comb(2*r,r)*lcm(*range(1,r+1))
        assert bound%D==0
        for R in (r,r+1,r+3):
            V=s.Matrix([[j**(2*k) for j in range(1,R+1)] for k in range(r)])
            assert V.rank()==r;ranks+=1
        cases.append({'r':r,'order':2*r+1,'frequency_max':r,'coefficients':[str(x) for x in c],
          'kernel_constant':str(C),'rational_common_denominator':str(D),'common_denominator_bits':D.bit_length(),
          'coefficient_integer_max_bits':max(abs(int(v*D/j)).bit_length() for j,v in enumerate(c,1))})
    # Symbolic mixed-scale anchor and jet tests (replace pi by a formal symbol).
    x,p=s.symbols('x p',nonzero=True)
    shared_checks=0
    for r in range(1,5):
        c=coeffs(r)
        phi=sum(s.Rational(v.numerator,v.denominator)*s.sin(2*p*j*x)/(2*p*j) for j,v in enumerate(c,1))
        H=x-phi
        jet=s.series(H,x,0,2*r+2).removeO().expand()
        assert all(jet.coeff(x,k)==0 for k in range(0,2*r+1))
        expected=s.Rational(4**r,comb(2*r,r))*p**(2*r)/(2*r+1)
        assert s.simplify(jet.coeff(x,2*r+1)-expected)==0
        for T in (s.Rational(1),s.Rational(2),s.Rational(5),s.Rational(3,2)):
            q=2 # T*q is integral for every selected T.
            for n in range(-4,5):
                assert s.simplify((H.subs({x:T*q*n,p:s.pi})/T)-q*n)==0
                shared_checks+=1
    out={'status':'PASS','coefficient_families':len(cases),'nonconstant_moment_checks':moments,
         'exact_vandermonde_rank_checks':ranks,'shared_scale_anchor_checks':shared_checks,
         'cases':cases,'scope':'Finite exact identities for separately proved all-order shared-modulator construction; no pi approximation theorem.'}
    Path(__file__).with_name('checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='cases'},indent=2))

if __name__=='__main__':run()
