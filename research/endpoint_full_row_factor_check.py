"""Exact finite full-row test for the diagonal-factor construction.

This does not check the separated-weight hypotheses or prove the H-to-infinity
interpolation theorem. It uses rational centres and standard-library fractions.
"""

from fractions import Fraction as F
from math import comb

def mul(p,q,cut):
 o=[F(0)]*(cut+1)
 for i,a in enumerate(p):
  for j,b in enumerate(q):
   if i+j<=cut:o[i+j]+=a*b
 return o

def ent(col,row,r1,r2):
 h,a,b=col;j,s,da,db=row
 if a<da or b<db:return F(0)
 p=[F(comb(h,k)) if k<=h else F(0) for k in range(s+1)]
 for r,e in ((r1,a-da),(r2,b-db)):
  G=[F(j*r)]+[F((-1)**(k+1),k) for k in range(1,s+1)]
  for _ in range(e):p=mul(p,G,s)
 return comb(a,da)*comb(b,db)*p[s]

def add_basis(basis,column):
 x=column[:]
 for pivot,vec in basis:
  if x[pivot]:
   z=x[pivot]
   x=[a-z*b for a,b in zip(x,vec)]
 try:pivot=next(i for i,a in enumerate(x) if a)
 except StopIteration:return False
 scale=x[pivot];x=[a/scale for a in x]
 basis.append((pivot,x));basis.sort(key=lambda t:t[0])
 return True

H=12; S=(H-1)//2; h=S-1
rows=[(j,s,a,b) for j in range(2) for s in range(H) for a in range(H) for b in range(H) if 2*s+4*a+8*b<H]
cols=[(z,a,b) for z in range(H+1) for a in range(H//2+1) for b in range(H//4+1) if z+2*a+4*b<=H]

def special(i,r1,r2):
 out=[]
 for row in rows:
  total=F(0)
  for power in range(h+1):
   coefficient=F((-1)**(h-power)*comb(h,power))
   col=(power,1,0) if i==0 else (power,0,1)
   total+=coefficient*ent(col,row,r1,r2)
  out.append(total)
 return out

p1=special(0,1,2);p2=special(1,1,2)
basis=[]
assert add_basis(basis,p1) and add_basis(basis,p2)
selected=[]
for col in cols:
 vec=[ent(col,row,1,2) for row in rows]
 if add_basis(basis,vec):selected.append(col)
 if len(basis)==len(rows):break
assert len(basis)==len(rows)
assert len(selected)==len(rows)-2
assert special(0,3,3)==special(1,3,3)
for r1 in range(1,4):
 for r2 in range(1,4):
  a=special(0,r1,r2);b=special(1,r1,r2)
  assert all(y-x == F((r2-r1)*(j if s==h and da==0 and db==0 else 0))
             for (j,s,da,db),x,y in zip(rows,a,b))
assert all(b-a == ((j if s==h and da==0 and db==0 else 0)) for (j,s,da,db),a,b in zip(rows,p1,p2))
print('H=',H,'h=',h,'rows=',len(rows),'columns=',len(cols))
print('Full row rank with two engineered columns; diagonal factor verified; monomial extensions=',len(selected))
