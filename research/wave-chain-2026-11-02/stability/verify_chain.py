"""Exact witnesses for an all-index invariant barrier and integer winding chain.

The all-index conclusion is proved in NOTE.md, not inferred from this finite run.
Only Python's standard library is needed. All acceptance tests are rational.
"""
from fractions import Fraction as F
from pathlib import Path
from math import gcd
import json


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(p, q):
    out = [F(0)] * max(len(p), len(q))
    for j, c in enumerate(p): out[j] += c
    for j, c in enumerate(q): out[j] += c
    return trim(out)


def scale(p, c): return trim([c*x for x in p])


def mul(p, q):
    out = [F(0)] * (len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): out[i+j] += a*b
    return trim(out)


def power(p, n):
    result = [F(1)]
    for _ in range(n): result = mul(result, p)
    return result


def l1(p): return sum(map(abs, p), F(0))


def next_state(C, R, lam):
    return 5*C, add(mul([F(1), F(0), F(1)], R), [-lam*C, F(C)])


def model(n, lam, growing=True):
    C, R = 2, [F(1)]
    for _ in range(n): C, R = next_state(C, R, lam)
    k = n+1 if growing else 2
    carrier = [F(0)]*k + [F(C)]
    return C, R, add(carrier, R), k


def complex_mul(a,b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def evaluate(p,z):
    v=(F(0),F(0))
    for a in reversed(p):
        v=complex_mul(v,z)
        v=(v[0]+a,v[1])
    return v


def circle(t):
    t=t%1
    quadrant=int(4*t)
    u=4*t-quadrant
    a,b=(1-u*u)/(1+u*u),2*u/(1+u*u)
    return [(a,b),(-b,a),(-a,-b),(b,-a)][quadrant]


def cross(a,b): return a[0]*b[1]-a[1]*b[0]
def dot(a,b): return a[0]*b[0]+a[1]*b[1]


def winding(points):
    assert points[0]==points[-1]
    total=0
    for a,b in zip(points,points[1:]):
        if cross(a,b)==0 and dot(a,b)<=0:
            raise ValueError('polygon passes through zero')
        c=cross(a,b)
        if a[1]<=0<b[1] and c>0: total+=1
        if b[1]<=0<a[1] and c<0: total-=1
    return total


def main():
    rows=[]
    affine_checks=0
    exact_bound_checks=0
    A=[F(1),F(0),F(1)]
    DminusA=[F(4),F(0),F(-1)]
    parameters=[F(-1),F(-7,9),F(-1,3),F(0),F(2,7),F(1,2),F(1)]
    for lam in parameters:
        C,R=2,[F(1)]
        for n in range(41):
            rho=F(2,5)
            bound=rho**n/F(2)+(1+abs(lam))/3*(1-rho**n)
            assert l1(R)/C<=bound<=F(2,3)
            # Rational parameter denominators do not grow with n.
            v=lam.denominator
            assert all((v*c).denominator==1 for c in R)
            exact_bound_checks+=1
            # Independent closed identity (5-A)R_n = A^n(5-A)+2B(5^n-A^n).
            An=power(A,n)
            right=add(mul(An,DminusA),scale(mul([-lam,F(1)],add([F(5**n)],scale(An,-1))),2))
            assert mul(DminusA,R)==right
            if n in (0,1,2,5,10,20,40):
                rows.append({'n':n,'lambda':str(lam),'C':str(C),
                             'residual_l1_over_C':str(l1(R)/C),
                             'proved_upper':str(bound),
                             'uniform_boundary_margin_at_least':'1/3',
                             'growing_chain_winding':n+1,
                             'constant_chain_winding':2,
                             'parameter_clear_factor':v,
                             'max_residual_integer_bits':max(abs(int(v*x)).bit_length() for x in R)})
            C,R=next_state(C,R,lam)
    # Check exact affine parameter transport, independently of the norm bound.
    for n in range(21):
        _,Rzero,_,_=model(n,F(0))
        _,Rone,_,_=model(n,F(1))
        for lam in parameters:
            _,R,_,_=model(n,lam)
            assert R==add(Rzero,scale(add(Rone,scale(Rzero,-1)),lam))
            affine_checks+=1
    # Certified sampled winding, independently of the symbolic Rouche argument.
    winding_rows=[]
    for n in (0,2,4):
        for lam in (F(-1),F(1)):
            C,R,P,k=model(n,lam)
            P=scale(P,F(1,C))
            # Four rational circle arcs have speed at most 8.
            L=8*sum(j*abs(c) for j,c in enumerate(P))
            N=16
            while L/(2*N)>=F(1,3): N*=2
            points=[evaluate(P,circle(F(j,N))) for j in range(N)]
            points.append(points[0])
            observed=winding(points)
            assert observed==k
            winding_rows.append({'n':n,'lambda':str(lam),'N':N,'L':str(L),
                                 'interpolation_error_upper':str(L/(2*N)),
                                 'certified_winding':observed})
    result={'status':'PASS','scope':'Finite exact verification of a separately proved all-n invariant barrier. No pi approximation claim.',
            'exact_bound_and_closed_identity_checks':exact_bound_checks,
            'affine_parameter_checks':affine_checks,
            'independent_certified_winding_checks':len(winding_rows),
            'examples':rows,'winding_examples':winding_rows}
    Path(__file__).with_name('checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('examples','winding_examples')},indent=2))


if __name__=='__main__': main()
