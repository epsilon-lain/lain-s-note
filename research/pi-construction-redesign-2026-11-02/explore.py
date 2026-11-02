from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import comb
import json,time
@lru_cache(None)
def moment(k,c):
 if c==0:return F(0),F(1,2*k+1)
 if k:
  a,b=moment(k-1,c-1);u,v=moment(k-1,c)
  return a-u,b-v
 if c==1:return F(1,4),F(0)
 a,b=moment(0,c-1)
 return a*F(2*c-3,2*c-2),b*F(2*c-3,2*c-2)+F(1,(2*c-2)*2**(c-1))
def atan(a,K=200):
 s=sum((F((-1)**k,(2*k+1)*a**(2*k+1)) for k in range(K)),F(0));t=F(1,(2*K+1)*a**(2*K+1));return s,s+t
l,h=atan(5);u,v=atan(239);lo,hi=16*l-4*v,16*h-4*u
rows=[];best=F(10**200);records=[];start=time.monotonic()
for budget in range(1,49):
 for a in range(0,min(16,budget)+1):
  for b in range(0,min(16,budget-a)+1):
   c=budget-a-b
   if not 1<=c<=24:continue
   A=B=F(0)
   for j in range(b+1):
    x,y=moment(a+j,c);coef=(-1)**j*comb(b,j);A+=coef*x;B+=coef*y
   if not A:continue
   r=-B/A
   if r<=0:continue
   p,q=r.numerator,r.denominator
   if r<lo:L,H=(lo-r)*q*q,(hi-r)*q*q
   elif r>hi:L,H=(r-hi)*q*q,(r-lo)*q*q
   else:raise AssertionError('pi interval too wide')
   row={'a':a,'b':b,'c':c,'A':str(A),'B':str(B),'p':p,'q':q,'quality_interval':[str(L),str(H)],'quality_display':float((L+H)/2)}
   rows.append(row)
   if H<best:
    best=H;records.append(row)
rows.sort(key=lambda r:r['quality_display'])
result={'scope':'Finite redesign experiment, not new pi record or asymptotic theorem','count':len(rows),'elapsed_seconds':time.monotonic()-start,'record_improvements':records,'best_10':rows[:10]}
open(Path(__file__).with_name('search_results.json'),'w').write(json.dumps(result,indent=2)+'\n')
print(json.dumps({**result,'best_10':[{k:v for k,v in r.items() if k!='quality_interval'} for r in rows[:10]],'record_improvements':[{k:v for k,v in r.items() if k!='quality_interval'} for r in records]},indent=2))
