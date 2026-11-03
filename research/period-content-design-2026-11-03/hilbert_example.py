"""Exact rational Schur lifts, globally integerized without a quality penalty."""
import json
from pathlib import Path
from fractions import Fraction as F
from math import gcd,lcm,comb
import sympy as s
ROOT=Path(__file__).parent

def inner(a,b):
 return sum((F(x)*F(y)/F(i+j+1) for i,x in enumerate(a) if x for j,y in enumerate(b) if y),F(0))
def dotrows(rows,v):
 return [sum(F(a)*b[j] for a,b in zip(v,rows)) for j in range(len(rows[0]))]
def frac(x):return F(int(x.p),int(x.q))

def verify(block,row):
 m=block['m'];v=row['quotient_coordinates']; CH=block['period_lift_H']; CB=block['period_lift_coordinates']; ZH=block['zero_period_H'];ZB=block['zero_period_lift_coordinates']
 h0=dotrows(CH,v);r0=dotrows(CB,v)
 A=s.Matrix([[s.Rational(inner(a,b)) for b in ZH] for a in ZH]); b=s.Matrix([s.Rational(inner(a,h0)) for a in ZH]); u=list(map(frac,A.inv()*(-b)))
 h=[x+y for x,y in zip(h0,dotrows(ZH,u))];r=[x+y for x,y in zip(r0,dotrows(ZB,u))]
 assert all(inner(h,z)==0 for z in ZH)
 projected=inner(h,h)
 assert projected==F(row['projected_hilbert_norm_squared'])
 assert r[-1]==int(row['alpha']);assert sum(r[:-1])==sum(map(int,row['R_coefficients']))
 scale=lcm(*(x.denominator for x in h+r)); hi=[int(x*scale) for x in h];ri=[int(x*scale) for x in r]
 content=gcd(*hi,*ri);hi=[x//content for x in hi];ri=[x//content for x in ri]
 alpha=ri[-1];z=sum(ri[:-1]);pfrac=F(-4*z,2**(4*m)*alpha)
 assert pfrac==F(int(row['p']),int(row['q']))
 # Rebuild the polynomial identity from integer primitive coordinates.
 c=4*m+1; R=[0]+ri[:-1]; N=[0]*(2*c+1)
 for k,a in enumerate(R):
  if k:
   N[k-1]+=k*a;N[k+1]+=(k-2*(c-1))*a
 for j in range(c):N[2*j]+=alpha*comb(c-1,j)
 target=[0]*len(N)
 for i,a in enumerate(hi):
  for j in range(2*m+1):target[i+2*m+j]+=a*(-1)**j*comb(2*m,j)
 assert N==target
 assert gcd(*ri)==1
 norm=inner(hi,hi);q=pfrac.denominator
 certificate=F(72,7)*F(q**4,alpha**2)*F(4,121)**(2*m)*norm
 expected=F(row['joint_quality_squared_upper'])*projected/F(row['hilbert_norm_squared'])
 assert certificate==expected
 assert F(row['quality_upper'])**2<certificate
 return {'m':m,'p':str(pfrac.numerator),'q':str(q),'old_bound_squared':row['joint_quality_squared_upper'],'new_bound_squared':str(certificate),'new_bound_display':float(certificate)**.5,'old_bound_display':float(F(row['joint_quality_squared_upper']))**.5,'integer_scale':str(scale),'removed_content':str(content),'alpha':str(alpha),'z':str(z),'period_gcd':str(gcd(alpha,z)),'raw_fraction_gcd':str(gcd(4*z,2**(4*m)*alpha)),'R_coefficients':list(map(str,ri[:-1])),'H_coefficients':list(map(str,hi)),'scope':'Same fraction; exact rational optimum attained by global integerization. Actual quality unchanged.'}

def main():
 data=json.loads((ROOT.parent/'pi-period-quotient-design-2026-11-02'/'period-main-kernel'/'hilbert_forced_period_lift.json').read_text()); results=[]
 for b in data['per_m']:
  results.append(verify(b,b['best']))
 print(json.dumps([{k:v for k,v in r.items() if k in ['m','p','q','old_bound_display','new_bound_display','integer_scale','removed_content']} for r in results],indent=2))
 (ROOT/'schur_exact_results.json').write_text(json.dumps({'status':'passed','results':results},indent=2)+'\n')
if __name__=='__main__':main()
