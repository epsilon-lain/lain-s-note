from endpoint_design import contour
from fractions import Fraction as F
from pathlib import Path
import json
import mpmath as mp
import sympy as sy
from symmetric_contour import integral
mp.mp.dps=100
x=sy.symbols('x')
cases=[(1,0,0,1,1),(1,1,2,1,2),(1,2,3,1,1),(1,1,3,2,3),(2,1,3,2,1),(2,2,4,3,1),(1,2,4,1,3),(2,1,3,5,7)]
for a,b,c,s,t in cases:
 A,B,_,(L,u,v),N,D=contour(a,b,c,s,t)
 f=sy.Poly.from_list(list(reversed(N)),x).as_expr()/sy.Poly.from_list(list(reversed(D)),x).as_expr()
 # Independent symbolic antiderivative/quadrature check of the pi branch.
 def fun(y):
  z=mp.mpc(u,y)
  q=((z-u)**2+v*v)*((z-(L-u))**2+v*v)
  return q**a*(2*z-L)**(2*b)/(z**(c+1)*(L-z)**(c+1))
 numerical=mp.quad(fun,[-v,0,v])
 exact=mp.mpf(A.numerator)/A.denominator*mp.pi+mp.mpf(B.numerator)/B.denominator
 assert abs(numerical-exact)<mp.mpf('1e-80')*max(1,abs(exact))
 # Reflection parameter equivalence: s,t -> t,2s.
 A2,B2,*_=contour(a,b,c,t,2*s)
 scale=F(2)**(4*a+2*b-2*c-1)
 assert A2==scale*A and B2==scale*B
# Compare fixed original geometry, accounting for center factor 4^b.
for a,b,c in [(1,1,2),(1,2,3),(2,1,3),(3,2,4)]:
 A,B,*_=integral(a,b,c)
 A2,B2,*_=contour(a,b,c,1,2)
 assert (A2,B2)==(4**b*A,4**b*B)
result={'status':'PASS','independent_100_decimal_quadratures':len(cases),'parameter_reflection_exact_checks':len(cases),'original_geometry_exact_checks':4,'scope':'Exact arithmetic identities supplemented by independent high-precision numerical contour checks; numerical quadrature does not prove the general formulas.'}
Path(__file__).with_name('contour_identity_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
