from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import comb,gcd
import json
@lru_cache(None)
def J(k,c):
 if not c:return F(0),F(1,2*k+1)
 if k:
  a,b=J(k-1,c-1);u,v=J(k-1,c);return a-u,b-v
 if c==1:return F(1,4),F(0)
 a,b=J(0,c-1);return a*F(2*c-3,2*c-2),b*F(2*c-3,2*c-2)+F(1,(2*c-2)*2**(c-1))
def integral(a,b,c):
 return tuple(sum(((-1)**j*comb(b,j)*J(a+j,c)[i] for j in range(b+1)),F(0)) for i in (0,1))
def atan(a):
 s=sum((F((-1)**k,(2*k+1)*a**(2*k+1)) for k in range(100)),F(0));return s,s+F(1,201*a**201)
l,h=atan(5);u,v=atan(239);lo,hi=16*l-4*v,16*h-4*u
best=None;count=0;record=[]
for budget in range(1,22):
 for a in range(8):
  for b in range(8):
   c=budget-a-b
   if not 1<=c<=8:continue
   A,B=integral(a,b,c);C,D=integral(a+1,b,c)
   for u in range(-12,13):
    for v in range(1,13):
     if gcd(u,v)!=1:continue
     alpha=u*A+v*C;beta=u*B+v*D
     if not alpha:continue
     r=-beta/alpha
     if r<=0:continue
     count+=1;p,q=r.numerator,r.denominator
     if r<lo:L,H=(lo-r)*q*q,(hi-r)*q*q
     elif r>hi:L,H=(r-hi)*q*q,(r-lo)*q*q
     else:raise AssertionError('pi interval too wide')
     if best is None or H<best:
      best=H;row={'a':a,'b':b,'c':c,'u':u,'v':v,'alpha':str(alpha),'beta':str(beta),'p':p,'q':q,'quality_interval':[str(L),str(H)],'quality_display':float((L+H)/2)};record.append(row)
result={'count':count,'coeff_bounds':{'a_b':[0,7],'c':[1,8],'u':[-12,12],'v':[1,12],'a_b_c_sum_max':21},'records':record,'best':record[-1],'scope':'Bounded exact experiment; recovered known approximants, not a new pi record or infinite construction.'}
open(Path(__file__).with_name('combination_results.json'),'w').write(json.dumps(result,indent=2)+'\n')
print(json.dumps({**result,'records':[{k:v for k,v in x.items() if k!='quality_interval'} for x in record],'best':{k:v for k,v in record[-1].items() if k!='quality_interval'}},indent=2))
