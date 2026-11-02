"""Exact checks for the wave-language reconstruction, not a check of the full source proof."""
import json
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import sympy as s
x,z=s.symbols('x z')
checks={}
# Sharp examples of the grouped Taylor occupancy lower bound.
cases=0
for sizes in [(1,),(2,),(3,),(4,),(2,2),(3,2),(1,2,3)]:
    blocks=[]
    for n in sizes:
        blocks.append(s.Matrix([[(x*j)**d for d in range(n)] for j in range(n)]))
    D=s.factor(s.diag(*blocks).det())
    expected=sum(n*(n-1)//2 for n in sizes)
    poly=s.Poly(D,x)
    order=min(m[0] for m,c in poly.terms() if c)
    assert order==expected
    cases+=1
checks['sharp_grouped_order_examples']=cases
# The actual f_{a,P}=binomial(alpha,a)*exp(h*z)*z**(alpha-a),
# evaluated at lambda*j. A truncated exponential gives exact low coefficients.
actual=[]
columns=[(0,0),(1,0),(0,1),(1,1),(0,2),(1,2)]
for labels in [[(0,0),(0,1),(0,2),(1,0)],[(0,0),(0,1),(1,0),(1,1)],[(0,0),(0,1),(0,2)]]:
    M=len(labels);degree=M*(M-1)//2+2
    matrix=[]
    for a,j in labels:
        row=[]
        for h,alpha in columns[:M]:
            if a>alpha:row.append(s.Integer(0));continue
            f=s.binomial(alpha,a)*sum(F(h**d,s.factorial(d))*z**(d+alpha-a) for d in range(degree+1))
            row.append(s.expand(f.subs(z,x*j)))
        matrix.append(row)
    D=s.Poly(s.Matrix(matrix).det(),x)
    ns={a:sum(aa==a for aa,j in labels) for a,j in labels}
    expected=sum(n*(n-1)//2 for n in ns.values())
    assert all(D.nth(d)==0 for d in range(expected))
    actual.append({'groups':list(ns.values()),'guaranteed_order':expected,'identically_zero':D.is_zero})
checks['actual_exponential_polynomial_checks']=actual
# The exact nu=2 gap is negative whenever A^2 < theta < 1.
barrier=0
for den in range(2,45):
    for num in range(1,den):
        theta=F(num,den)
        for anum in range(1,den):
            A=F(anum,den)
            if theta<A<1 and A*A<theta:
                assert 2*(A-theta)-(1-theta)<0
                barrier+=1
checks['exact_boundary_checks']=barrier
# A feasible rational example strictly above the critical exponent.
nu=F(3);theta=F(9,10);A=F(47,50);C=F(21,20);B=F(943,1000);eta=F(1,1000)
assert 0<theta<A<B<1 and C>1 and B<1/C and B<C*theta<1
G=nu*(A*(1-eta)-theta)-(1-theta)
assert G>0 and A*A<theta
checks['feasible_above_two']={'nu':str(nu),'theta':str(theta),'A':str(A),'B':str(B),'C':str(C),'eta':str(eta),'gap':str(G)}
checks['status']='PASS'
checks['scope']='Algebraic mechanism and finite exact controls only; interpolation and the full source theorem are not independently certified.'
Path(__file__).with_name('checks.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(checks,ensure_ascii=False,indent=2))
