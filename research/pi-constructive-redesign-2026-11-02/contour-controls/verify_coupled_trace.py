from coupled_trace import coupled
from quadratic_trace_contour import K
from pathlib import Path
import sympy as sy
import mpmath as mp
import json
mp.mp.dps=100
root=mp.sqrt(2)
def evaluate(expr,rho):
 return mp.mpf(str(sy.N(K.to_sympy(expr).subs(sy.sqrt(2),sy.Float(str(rho),110)),105)))
cases=[(0,2,1,2),(2,0,1,2),(1,2,1,3),(2,1,1,4),(3,3,2,6),(4,1,0,5)]
for a,b,d,c in cases:
 value,r0,B=coupled(a,b,d,c)
 total=mp.mpc(0)
 for rho in [root,-root]:
  w=rho-1;other=-rho-1
  residue=evaluate(r0,rho)
  def integrand(t):
   z=mp.mpc(1,t)
   return ((z-1)**2+w*w)**a*((z-1)**2+other*other)**b*(z-1)**(2*d)/(z**(c+1)*(2-z)**(c+1))
  J=mp.quad(integrand,[-w,0,w])
  expected=(mp.pi/2 if rho>0 else -3*mp.pi/2)*residue+evaluate(B,rho)
  assert abs(J-expected)<mp.mpf('1e-75')*max(1,abs(expected))
  total+=J/residue
 assert abs(total-(-mp.pi+mp.mpf(value.numerator)/value.denominator))<mp.mpf('1e-70')
result={'status':'PASS','cases':cases,'embedding_quadratures':2*len(cases),'scope':'Independent 100-decimal conjugate-contour verification of coupled exact field formulas; not a general proof or an infinite quality result.'}
Path(__file__).with_name('coupled_trace_identity_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
