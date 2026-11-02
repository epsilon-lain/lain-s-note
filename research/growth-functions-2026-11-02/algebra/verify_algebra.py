"""Exact implementation controls for a separately proved growth-function algebra."""
from fractions import Fraction as F
from pathlib import Path
import json
import random
import sympy as sy

# Polynomials in formal variables t and B, with integer coefficients.
def clean(p):return {k:v for k,v in p.items() if v}
def plus(p,q):
    out=p.copy()
    for m,c in q.items():out[m]=out.get(m,0)+c
    return clean(out)
def scale(a,p):return clean({m:a*c for m,c in p.items()})
def times(p,q):
    out={}
    for (i,j),c in p.items():
        for (k,l),d in q.items():out[(i+k,j+l)]=out.get((i+k,j+l),0)+c*d
    return clean(out)
def shiftB(p):return {(i,j+1):c for (i,j),c in p.items()}
def add(A,C):
    a,u=A;b,v=C
    return a+b,plus(u,v)
def multiply(A,C):
    a,u=A;b,v=C
    return a*b,plus(plus(scale(a,v),scale(b,u)),shiftB(times(u,v)))
def gain(A):
    a,u=A
    return plus({(0,0):a},shiftB(u))
def val(p,t,B):return sum((F(c)*t**i*B**j for (i,j),c in p.items()),F(0))
def degree(p,d0):return max((i+d0*j for i,j in p),default=0)
def basis(t,order):
    fractional=t-t.numerator//t.denominator
    return fractional**(order+1)*(1-fractional)**(order+1)
def evaluate(A,t,order):return t*val(gain(A),t,basis(t,order))

def run():
    rng=random.Random(20261102)
    trials=80
    symbol_checks=0;sample_checks=0;den_checks=0;cost_checks=0
    for _ in range(trials):
        a,b,c=[rng.randrange(-9,10) for j in range(3)]
        def randpoly():
            p={}
            for j in range(5):
                m=rng.randrange(4),rng.randrange(3)
                p[m]=p.get(m,0)+rng.randrange(-7,8)
            return clean(p)
        A=(a,randpoly());C=(b,randpoly());D=(c,randpoly())
        assert gain(add(A,C))==plus(gain(A),gain(C))
        assert gain(multiply(A,C))==times(gain(A),gain(C))
        assert multiply(A,C)==multiply(C,A)
        assert multiply(multiply(A,C),D)==multiply(A,multiply(C,D))
        assert multiply(A,add(C,D))==add(multiply(A,C),multiply(A,D))
        assert multiply(A,(1,{}))==A
        assert add(A,(-a,scale(-1,A[1])))==(0,{})
        symbol_checks+=7
        for order in (0,1,2):
            d0=2*order+2
            for pair in (A,C,add(A,C),multiply(A,C)):
                label,u=pair
                for n in range(-5,6):
                    assert evaluate(pair,F(n),order)==label*n
                    sample_checks+=1
                Udegree=degree(u,d0)
                N=len(u);H=max(map(abs,u.values()),default=0)
                for p,q in ((-19,7),(-1,3),(0,1),(2,5),(17,9),(21,4)):
                    t=F(p,q);q=t.denominator;p=t.numerator
                    E=Udegree+d0+1
                    out=evaluate(pair,t,order)
                    assert q**E % out.denominator==0
                    den_checks+=1
                    # With common denominator q**E, numerator obeys the proved bound.
                    numerator=out.numerator*(q**E//out.denominator)
                    M=max(abs(p),q,1)
                    assert abs(numerator)<=(abs(label)+N*H)*M**E
                    cost_checks+=1
                if u:
                    new=multiply(A,C)[1]
                    assert degree(new,d0)<=max(degree(A[1],d0),degree(C[1],d0),d0+degree(A[1],d0)+degree(C[1],d0))
    t=sy.symbols('t',positive=True)
    divergence=sy.limit((t+sy.sin(sy.pi*t)**sy.Rational(1,4))*(2*t+sy.sin(sy.pi*t)**sy.Rational(1,4))/t,t,0,dir='+')
    assert divergence==sy.oo
    # Stable sine member belongs to the complete normal form, including derivative 0.
    normal_checks=0
    for a in range(1,8):
        u=-a*sy.cos(sy.pi*t)
        normal=a*t+sy.sin(sy.pi*t)*u/sy.pi
        actual=a*t-a*sy.sin(2*sy.pi*t)/(2*sy.pi)
        assert sy.trigsimp(normal-actual)==0
        assert sy.limit(sy.diff(normal,t),t,0)==0
        normal_checks+=2
    # The rational B_s backend has derivatives 0 through s on both cell boundaries.
    boundary_checks=0
    for order in range(7):
        B=t**(order+1)*(1-t)**(order+1)
        for d in range(order+1):
            dB=sy.diff(B,t,d)
            assert dB.subs(t,0)==0 and dB.subs(t,1)==0
            boundary_checks+=1
    output={'status':'PASS','scope':'Exact controls for the proved smooth-clock algebra and a finite rational expression grammar; not certification of all arbitrary functions or a pi approximation construction.',
      'formal_polynomial_identity_checks':symbol_checks,'integer_readout_checks':sample_checks,
      'rational_denominator_divisibility_checks':den_checks,'rational_bit_bound_controls':cost_checks,
      'sine_normal_form_checks':normal_checks,'rational_backend_boundary_checks':boundary_checks,
      'continuous_product_counterexample_limit':'infinity',
      'dependencies':{'python':'3','sympy':sy.__version__}}
    Path(__file__).with_name('checks.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(output,ensure_ascii=False,indent=2))
if __name__=='__main__':run()
