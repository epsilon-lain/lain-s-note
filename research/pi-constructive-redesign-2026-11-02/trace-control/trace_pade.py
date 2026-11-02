"""Fixed Q(sqrt(2)) design; no continued fractions or target fractions as input."""
from fractions import Fraction as F
from pathlib import Path
from math import gcd
import json

def add(a,b):return a[0]+b[0],a[1]+b[1]
def mul(a,b):return a[0]*b[0]+2*a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def scale(a,k):return a[0]*k,a[1]*k

def atan(a,K=500):
 s=sum((F((-1)**k,(2*k+1)*a**(2*k+1)) for k in range(K)),F(0))
 return s,s+F(1,(2*K+1)*a**(2*K+1))
def run(N=100):
 l,h=atan(5);u,v=atan(239);lo,hi=16*l-4*v,16*h-4*u
 x=(-1,1);x2=mul(x,x);Ap,A=(1,0),(1,0);Bp,B=(0,0),x
 records=[];rows=[]
 for n in range(1,N+1):
  if n>1:
   k=n-1;Ap,A=A,add(scale(A,2*k+1),scale(mul(x2,Ap),k*k));Bp,B=B,add(scale(B,2*k+1),scale(mul(x2,Bp),k*k))
  norm=A[0]**2-2*A[1]**2
  assert norm
  r=F(-8*(B[0]*A[0]-2*B[1]*A[1]),norm);p,q=r.numerator,r.denominator
  if r<lo:L,H=(lo-r)*q*q,(hi-r)*q*q
  elif r>hi:L,H=(r-hi)*q*q,(r-lo)*q*q
  else:raise AssertionError('need more pi terms')
  row={'n':n,'p':str(p),'q':str(q),'denominator_digits':len(str(q)),'quality_display':float(H) if H<10**300 else '>1e300','Q_below_one':H<1}
  rows.append(row)
  if n in (1,2,3,4,5,10,20,50,100):records.append(row)
 out={'scope':'Constructive quadratic-trace baseline, not an infinite quality-zero theorem','indices':[1,N],'successful_indices':[r['n'] for r in rows if r['Q_below_one']],'samples':records,'minimum':min(rows,key=lambda r:float(r['quality_display']))}
 Path(__file__).with_name('trace_pade_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':run()
