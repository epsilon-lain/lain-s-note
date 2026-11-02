from symmetric_contour import integral
from fractions import Fraction as F
from math import log,ceil,lcm,gcd
from pathlib import Path
import json

def atan(a,N):
 s=sum((F((-1)**k,(2*k+1)*a**(2*k+1)) for k in range(N)),F(0));t=F(1,(2*N+1)*a**(2*N+1))
 return (s,s+t) if N%2==0 else (s-t,s)

def pi_interval(N):
 l,h=atan(5,N);u,v=atan(239,N);return 16*l-4*v,16*h-4*u

def run():
 families=[(1,2,3),(2,1,3),(6,3,10),(1,1,2),(2,2,3),(3,1,4),(3,2,4),(2,3,4),(1,0,1),(1,0,2),(1,1,3),(1,2,4),(2,1,4)]
 rows=[]
 for family in families:
  for n in range(1,13):
   a,b,c=[n*t for t in family]
   A,B,*_,univ=integral(a,b,c)
   if not A:continue
   r=-B/A;p,q=r.numerator,r.denominator
   D=lcm(A.denominator,B.denominator);alpha=int(D*A);beta=int(-D*B);g=gcd(alpha,beta)
   N=max(120,ceil((q.bit_length()*log(2)*2+80)/(2*log(5))))
   pl,ph=pi_interval(N)
   if r<pl:lo,hi=(pl-r)*q*q,(ph-r)*q*q
   elif r>ph:lo,hi=(r-ph)*q*q,(r-pl)*q*q
   else:raise AssertionError('pi interval too wide')
   assert lo>0 and hi/lo<F(1000000000001,1000000000000)
   rows.append({'family':family,'n':n,'logQ_over_n_endpoint_summaries':[(log(lo.numerator)-log(lo.denominator))/n,(log(hi.numerator)-log(hi.denominator))/n],
                'q_digits':len(str(q)),'actual_D_bits':D.bit_length(),'universal_D_bits':univ.bit_length(),'g_bits':g.bit_length(),
                'p':str(p),'q':str(q),'quality_below_one':hi<1})
 result={'status':'PASS','scope':'Finite scaling evidence, exact rationals and Machin error intervals; displayed logs are floating summaries, no asymptotic theorem.', 'families':families,'n_range':[1,12],'records':rows}
 Path(__file__).with_name('scale_checks.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'status':'PASS','n12':[{k:v for k,v in r.items() if k not in ('p','q')} for r in rows if r['n']==12]},indent=2))
if __name__=='__main__':run()
