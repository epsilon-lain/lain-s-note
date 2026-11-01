"""Finite truncated-series and modular matrix arithmetic, standard library only."""
from fractions import Fraction as F
from math import comb
from itertools import product
from collections import Counter
import random, json
from pathlib import Path

PRIME=1000003

def mul(x,y,N):
 z=[0]*(N+1)
 for i,u in enumerate(x):
  for j,v in enumerate(y[:N+1-i]):z[i+j]=(z[i+j]+u*v)%PRIME
 return z

def powers(x,k,N):
 out=[[1]+[0]*N]
 for _ in range(k):out.append(mul(out[-1],x,N))
 return out

def columns(m,H,weights,w0=F(1)):
 return [(h,a) for a in product(*(range(H//w+1) for w in weights))
         for h in range(int((H-sum(w*b for w,b in zip(weights,a)))/w0)+1)]

def rows(m,H,weights,K=2,theta=F(1,2)):
 def ceil_fraction(x):
  x=F(x)
  return -((-x.numerator)//x.denominator)
 return [(j,s,b) for j in range(K)
         for b in product(*(range(int(H*theta/w)+1) for w in weights))
         for s in range(max(0,ceil_fraction(H-sum(w*x for w,x in zip(weights,b))/theta)))]

def vectors(cols,rowlist,assignment,W=2):
 N=max(r[1] for r in rowlist)
 maxpow=max(sum(a) for h,a in cols)
 log=[0]+[((-1)**(k+1)*pow(k,-1,PRIME))%PRIME for k in range(1,N+1)]
 centers={}
 for j in set(r[0] for r in rowlist):
  base=log[:];base[0]=j*W%PRIME
  centers[j]=powers(base,maxpow,N)
 ans=[]
 for (j,s,beta),a in zip(rowlist,assignment):
  line=[]
  for h,alpha in cols:
   if any(aa>al for aa,al in zip(a,alpha)):
    line.append(0);continue
   factor=1
   for al,aa in zip(alpha,a):factor=factor*comb(al,aa)%PRIME
   poly=centers[j][sum(alpha)-sum(a)]
   line.append(factor*sum(comb(h,k)*poly[s-k] for k in range(min(h,s)+1))%PRIME)
  ans.append(line)
 return ans

def rank_pivots(matrix):
 a=[r[:] for r in matrix]; n=len(a); r=0;piv=[]
 for c in range(len(a[0])):
  k=next((i for i in range(r,n) if a[i][c]),None)
  if k is None:continue
  a[r],a[k]=a[k],a[r]
  inv=pow(a[r][c],-1,PRIME)
  a[r]=[x*inv%PRIME for x in a[r]]
  for i in range(r+1,n):
   f=a[i][c]
   if f:a[i]=[(x-f*y)%PRIME for x,y in zip(a[i],a[r])]
  piv.append(c);r+=1
  if r==n:break
 return r,piv

def cost_parts(rowlist,assignment,H,weights):
 linear=2*sum(sum(w*(x-b) for w,x,b in zip(weights,a,row[2]))
              for row,a in zip(rowlist,assignment))
 counts=Counter(assignment)
 collisions=sum(n*n for n in counts.values())
 return F(linear,len(rowlist)*H),F(collisions,len(rowlist)*H)

