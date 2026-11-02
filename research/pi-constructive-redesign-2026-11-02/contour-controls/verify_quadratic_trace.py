from quadratic_trace_contour import normalized_trace,K
from pathlib import Path
import sympy as sy
import mpmath as mp
import json
mp.mp.dps=100
root=mp.sqrt(2)

def evaluate(expr,r):
 return mp.mpf(str(sy.N(K.to_sympy(expr).subs(sy.sqrt(2),sy.Float(str(r),110)),105)))
cases=[(1,0,0),(1,1,2),(1,2,3),(2,1,3),(2,0,4),(3,1,4)]
for a,b,c in cases:
 value,r0,B=normalized_trace(a,b,c)
 integrals=[]
 for rho in [root,-root]:
  w=rho-1
  residue=evaluate(r0,rho)
  def integrand(t):
   z=mp.mpc(1,t)
   return ((z-1)**2+w*w)**(2*a)*(z-1)**(2*b)/(z**(c+1)*(2-z)**(c+1))
  J=mp.quad(integrand,[-w,0,w])
  expected=(mp.pi/2 if rho>0 else -3*mp.pi/2)*residue+evaluate(B,rho)
  assert abs(J-expected)<mp.mpf('1e-75')*max(1,abs(expected))
  integrals.append(J/residue)
 assert abs(sum(integrals)-(-mp.pi+mp.mpf(value.numerator)/value.denominator))<mp.mpf('1e-70')
result={'status':'PASS','parameter_cases':cases,'embedding_quadratures':2*len(cases),'scope':'Independent 100-decimal numerical checks of exact field formulas, including both signed conjugate branches.'}
Path(__file__).with_name('quadratic_trace_identity_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
