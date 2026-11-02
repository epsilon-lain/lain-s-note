from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json
from antiderivative_search import data,atan
l,h=atan(5);u,v=atan(239);lo,hi=16*l-4*v,16*h-4*u
rows=[]
for c in range(4,33):
 best=None
 for a in range(1,c):
  b=c-a;M,B=data(c,a,b);assert B.rows==1
  r=list(map(int,B[0,:]));alpha=r[-1];r1=sum(r[:-1])
  if not alpha:continue
  value=F(-4*r1,alpha*2**(c-1));p,q=value.numerator,value.denominator
  if value<lo:L,H=(lo-value)*q*q,(hi-value)*q*q
  elif value>hi:L,H=(value-hi)*q*q,(value-lo)*q*q
  else:raise AssertionError('narrow enclosure')
  if best is None or H<best[0]:best=(H,{'c':c,'a':a,'b':b,'p':str(p),'q':str(q),'Q_display':float(H),'alpha':str(alpha),'R_at_one':str(r1),'primitive_numerator_scale':str(next(int(x) for x in M*__import__('sympy').Matrix(r) if x))})
 rows.append(best[1])
Path(__file__).with_name('single_kernel_control.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
