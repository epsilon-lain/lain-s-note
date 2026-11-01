"""Exact explanatory models; finite checks do not prove a property of pi.

Run: python3 core_models.py > results_core_models.json
Requires SymPy and mpmath. Continued-fraction intervals are rational and
rigorous. Quadratures are explicitly numerical cross-checks only.
"""
from fractions import Fraction as Q
from math import factorial
import json
import sympy as sp
import mpmath as mp


def textq(x):
    return str(x.numerator) if x.denominator == 1 else str(x)


def niven_polynomial(a, b, n):
    """f=x^n(a-bx)^n/n! at a rational endpoint r=a/b."""
    x=sp.Symbol('x')
    f=sp.Poly(x**n*(a-b*x)**n/sp.factorial(n), x).as_expr()
    F=sum((-1)**j*sp.diff(f,x,2*j) for j in range(n+1))
    assert sp.expand(sp.diff(F,x,2)+F-f)==0
    r=sp.Rational(a,b)
    assert sp.expand(f.subs(x,r-x)-f)==0
    derivatives=[]
    for k in range(2*n+1):
        u=sp.diff(f,x,k).subs(x,0)
        v=sp.diff(f,x,k).subs(x,r)
        assert u.is_Integer and v.is_Integer and v==(-1)**k*u
        derivatives.append(int(u))
    F0=int(F.subs(x,0)); Fp0=int(sp.diff(F,x).subs(x,0))
    # This boundary value is NOT generally an integer. Period closure is
    # an additional condition; dropping it must not reproduce a false proof.
    symbolic=F0*(1-sp.cos(r))-Fp0*sp.sin(r)
    mp.mp.dps=80
    rr=mp.mpf(a)/b
    direct=mp.quad(lambda z:z**n*(a-b*z)**n/mp.factorial(n)*mp.sin(z),[0,rr])
    boundary=mp.mpf(F0)*(1-mp.cos(rr))-mp.mpf(Fp0)*mp.sin(rr)
    assert abs(direct-boundary)<mp.mpf('1e-55')*max(1,abs(boundary))
    return dict(a=a,b=b,n=n,endpoint_derivatives=derivatives,
                F0=F0,Fp0=Fp0,closed_period_integer=2*F0,
                actual_boundary_expression=str(symbolic),
                actual_integral_numeric=mp.nstr(direct,30),
                period_defect_numeric=mp.nstr(2*F0-boundary,30),
                numerical_scope='80-digit quadrature check; not interval certification')


def convergents(digits):
    pp,p=0,1; qq,q=1,0
    for n,a in enumerate(digits):
        pn=a*p+pp; qn=a*q+qq
        assert pn*q-p*qn == (-1)**(n+1)
        yield n,pn,qn
        pp,p=p,pn;qq,q=q,qn


def cf_error_interval(p,q,pp,qq,next_digit):
    """For the SAME convergent and any tail y in (a,a+1)."""
    a=next_digit
    k_lo=Q(q,(a+1)*q+qq)
    k_hi=Q(q,a*q+qq)
    # Independent Mobius endpoint calculation uses the exact determinant.
    assert abs(p*qq-pp*q)==1
    vals=[Q(p*a+pp,q*a+qq),Q(p*(a+1)+pp,q*(a+1)+qq)]
    direct=sorted(q*q*abs(v-Q(p,q)) for v in vals)
    assert direct==[k_lo,k_hi]
    return k_lo,k_hi


def cf_model(kind, length=25):
    digits=[1]+[1]*(length+1) if kind=='bounded' else [0]+list(range(1,length+2))
    conv=list(convergents(digits)); out=[]
    for n,p,q in conv[:-1]:
        if n<1:continue
        _,pp,qq=conv[n-1]
        lo,hi=cf_error_interval(p,q,pp,qq,digits[n+1])
        if kind=='bounded':
            assert lo>=Q(1,3)
        else:
            assert q>=factorial(n)
            assert lo>=Q(1,n+3) and hi<=Q(1,n+1)
        out.append(dict(n=n,a_next=digits[n+1],p=str(p),q=str(q),
                        normalized_error_lower=textq(lo),normalized_error_upper=textq(hi),
                        upper_decimal=float(hi)))
    return dict(kind=kind,definition='[1;1,1,...]' if kind=='bounded' else '[0;1,2,3,...]',
                rows=out,scope='Exact intervals for the infinite continued fraction; these are not pi data.')


def main():
    examples=[niven_polynomial(a,b,n) for a,b,n in
              [(3,1,1),(3,1,4),(22,7,1),(22,7,4),(355,113,2)]]
    # Symbolically isolate why the source's selected parameter architecture
    # cannot work at the endpoint, even with all secondary errors set to 0.
    A,d=sp.symbols('A d',positive=True)
    endpoint=sp.expand(2*A-(A*A+d)-1)
    assert sp.expand(endpoint+((1-A)**2+d))==0
    result=dict(niven=examples,continued_fractions=[cf_model('bounded'),cf_model('unbounded')],
                endpoint_identity='2A-(A^2+d)-1=-(1-A)^2-d < 0 for 0<A<1,d>0',
                exact_checks='polynomial ODE, endpoint integrality/symmetry, CF determinant and Mobius interval identities',
                limitation='No independent proof of the OpenAI interpolation theorem, irrationality exponent of pi, or pi BA/nonBA.')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
