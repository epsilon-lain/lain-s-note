from endpoint_design import contour,pibounds
from fractions import Fraction as F
from pathlib import Path
from math import log,lcm,gcd
import json
pl,ph=pibounds(500)
rows=[]
for family in [(1,1,2),(1,2,3),(2,1,3)]:
 for s,t in [(1,2),(1,1),(2,3),(5,7)]:
  for n in [1,2,4,8]:
   A,B,U,coords,*_=contour(*(n*v for v in family),s,t)
   if not A:continue
   r=-B/A;p,q=r.numerator,r.denominator
   if r<pl:lo,hi=(pl-r)*q*q,(ph-r)*q*q
   elif r>ph:lo,hi=(r-ph)*q*q,(r-pl)*q*q
   else:raise AssertionError('pi interval too wide')
   assert hi/lo<F(1000000000001,1000000000000)
   D=lcm(A.denominator,B.denominator);g=gcd(int(D*A),int(D*B))
   errorlo=lo/(q*q);errorhi=hi/(q*q)
   rows.append({'family':family,'s':s,'t':t,'n':n,'kappa':str(F(t,s)+F(2*s,t)),
                'logQ_over_n_summary':(log(lo.numerator)-log(lo.denominator))/n,
                'log_error_over_n_summary':(log(errorlo.numerator)-log(errorlo.denominator))/n,
                'log_q_over_n_summary':log(q)/n,'q_digits':len(str(q)),
                'D_bits':D.bit_length(),'g_bits':g.bit_length(),'p':str(p),'q':str(q)})
for fam in [(1,1,2),(1,2,3),(2,1,3)]:
 for n in [1,2,4,8]:
  a=next(r for r in rows if r['family']==fam and (r['s'],r['t'],r['n'])==(1,1,n))
  b=next(r for r in rows if r['family']==fam and (r['s'],r['t'],r['n'])==(1,2,n))
  assert (a['p'],a['q'])==(b['p'],b['q'])
result={'status':'PASS','scope':'Finite exact contour-coefficient arithmetic with rational Machin error bounds; all logarithms are approximate summaries, not rigorous enclosures or asymptotic claims.','records':rows}
Path(__file__).with_name('geometry_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':'PASS','n8':[{k:v for k,v in r.items() if k not in ['p','q']} for r in rows if r['n']==8]},indent=2))
