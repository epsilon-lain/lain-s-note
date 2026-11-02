"""Exact rational interval controls for the separately proved growth dynamics.
No floating-point comparison enters a PASS decision; displays alone use float.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json


def atan_reciprocal(q, terms):
    value=sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(terms)),F(0))
    rem=F(1,(2*terms+1)*q**(2*terms+1))
    return (value,value+rem) if terms%2==0 else (value-rem,value)


def add(A,B):return A[0]+B[0],A[1]+B[1]
def neg(A):return -A[1],-A[0]
def mul(A,B):
    values=[a*b for a in A for b in B]
    return min(values),max(values)
def pt(x):return F(x),F(x)
def scale(A,c):return mul(A,pt(c))
def divide(A,B):
    assert B[0]>0
    return mul(A,(1/B[1],1/B[0]))
def abs_upper(A):return max(abs(A[0]),abs(A[1]))

A=atan_reciprocal(5,36);B=atan_reciprocal(239,12)
PI=add(scale(A,16),scale(B,-4))
assert F(314,100)<PI[0]<PI[1]<F(22,7) and PI[1]**2<10
TWO_PI=scale(PI,2)


def sin_interval(X,degree=22):
    # Taylor polynomial through degree 2*degree+1; the even next term is zero.
    Y=mul(X,X)
    V=pt(F((-1)**degree,factorial(2*degree+1)))
    for j in range(degree-1,-1,-1):
        V=add(pt(F((-1)**j,factorial(2*j+1))),mul(Y,V))
    V=mul(X,V)
    R=abs_upper(X)**(2*degree+3)/factorial(2*degree+3)
    return add(V,(-R,R))


def correction(e):return divide(sin_interval(scale(TWO_PI,e)),TWO_PI)
def error_map(e):return add(pt(e),neg(correction(e)))


def rational_backend_checks():
    def centered(t):
        n=(t+F(1,2)).numerator//(t+F(1,2)).denominator
        return n,t-n
    def Q(t):
        n,e=centered(t)
        return n+4*e**3
    def psi(t):
        n,e=centered(t)
        return e**3*(1-4*e**2)
    def update(y,m):
        return (Q(y[0])+m[0]+psi(y[1])/4,
                Q(y[1])+m[1]-psi(y[0])/4)
    C=F(17,4);R=F(1,4);noise=F(47,256)
    tube=0;integer_jumps=0
    for j in range(-4,5):
        for k in range(-4,5):
            u=(F(j,16),F(k,16));H=update(u,(0,0))
            assert max(abs(t) for t in H)<=C*max(abs(t) for t in u)**3
            for z1 in (-noise,noise):
                for z2 in (-noise,noise):
                    assert max(abs(H[0]+z1),abs(H[1]+z2))<=R
                    tube+=1
            for N in [(-23,17),(0,0),(10**8,-10**9)]:
                for m in [(-2,3),(0,0),(10**7,-11)]:
                    result=update((u[0]+N[0],u[1]+N[1]),m)
                    assert result==(H[0]+N[0]+m[0],H[1]+N[1]+m[1])
                    integer_jumps+=1
    # Same polynomial value and first derivative at the wrap boundary.
    for e in (F(-1,2),F(1,2)):
        assert e**3*(1-4*e**2)==0
        assert 3*e**2-20*e**4==F(-1,2)
        assert 12*e**2==3
    assert 4*F(1,2)**3==F(1,2)
    # Fixed rational coefficient budget; common denominator fifth-power charge.
    denominator_checks=0
    for D in range(2,31):
        for b in range(1,D):
            u=(F(b,D),F(D-b,D))
            out=update(u,(1,1))
            assert all((4*D**5)%v.denominator==0 for v in out)
            denominator_checks+=1
    return {'status':'PASS','rational_coupled_tube_checks':tube,
            'variable_integer_jump_checks':integer_jumps,
            'fifth_power_denominator_checks':denominator_checks,
            'cubic_constant':str(C),'radius':str(R),'noise_budget':str(noise),
            'first_contraction_factor':str(C*R**2)}


def run():
    shrink=0
    for j in range(1,100):
        e=F(j,200)
        H=error_map(e)
        assert 0<H[0]<=H[1]<e
        assert H[1] <= F(20,3)*e**3
        shrink+=1
    noise=F(1,12); radius=F(1,6)
    tube_checks=0
    for j in range(-100,101):
        H=error_map(F(j,600))
        for eta in (-noise,F(0),noise):
            assert abs_upper(add(H,pt(eta)))<radius
            tube_checks+=1
    # Nonzero 2D coupling B = [[0,1/4],[-1/4,0]].
    # Cubic coupling uses C(e)=sin(2*pi*e)/(2*pi), not nearest-integer rounding.
    C=F(83,12); R=F(1,4)
    noise_budget=R-C*R**3
    coupled=0
    for j in range(-4,5):
        for k in range(-4,5):
            e1,e2=F(j,16),F(k,16)
            h1,h2=error_map(e1),error_map(e2)
            c1,c2=correction(e1),correction(e2)
            cube1=mul(mul(c1,c1),c1);cube2=mul(mul(c2,c2),c2)
            H1=add(h1,scale(cube2,F(1,4)))
            H2=add(h2,scale(cube1,F(-1,4)))
            r=max(abs(e1),abs(e2))
            assert max(abs_upper(H1),abs_upper(H2))<=C*r**3
            for eta1 in (-noise_budget,noise_budget):
                for eta2 in (-noise_budget,noise_budget):
                    assert max(abs_upper(add(H1,pt(eta1))),abs_upper(add(H2,pt(eta2))))<=R
                    coupled+=1
    envelope=[]
    r=R
    for k in range(7):
        expected=R*(C*R**2)**((3**k-1)//2)
        assert r==expected
        envelope.append({'step':k,'exact':str(r),'display':float(r)})
        r=C*r**3
    out={'status':'PASS','arithmetic':'Fraction intervals, Machin identity plus Taylor remainder; floats are display-only',
         'one_dimensional_shrink_and_cubic_checks':shrink,
         'bounded_noise_tube_checks':tube_checks,
         'coupled_two_dimensional_tube_checks':coupled,
         'normalized_coupling_radius':str(R),'rational_cubic_constant':str(C),
         'rational_noise_budget':str(noise_budget),'cubic_envelope':envelope,
         'scope':'Finite controls of separately proved all-time dynamics; no theorem of improved pi approximation.',
         'rational_C1_backend':rational_backend_checks()}
    Path(__file__).with_name('checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='cubic_envelope'},indent=2))

if __name__=='__main__':run()
