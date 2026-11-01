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

def main():
 rng=random.Random(2864)
 results=[]
 for m,H,ws,K,theta,w0 in [(1,4,(1,),4,F(1,8),F(1,4)),(1,8,(1,),2,F(1,4),F(1,4)),(2,6,(1,1),4,F(1,3),F(1,4))]:
  assert K*theta**m<1 and K*w0*theta**m<1
  cols=columns(m,H,ws,w0);rr=rows(m,H,ws,K,theta)
  print('case',m,H,ws,len(rr),len(cols),flush=True)
  base=[r[2] for r in rr]
  rank,piv=rank_pivots(vectors(cols,rr,base))
  if rank<len(rr):
   print('base rank deficient',rank,flush=True);continue
  best=base; cp=cost_parts(rr,base,H,ws)
  # Seek positive transverse shifts while preserving full rank.
  changes=0
  for trial in range(300):
   candidate=best[:]
   idx=rng.randrange(len(rr));j,s,beta=rr[idx]
   if j==0:continue
   i=rng.randrange(m)
   new=list(candidate[idx]);new[i]+=1
   if sum(w*a for w,a in zip(ws,new))>H:continue
   candidate[idx]=tuple(new)
   newcp=cost_parts(rr,candidate,H,ws)
   # Penalty uses 0.17<c so accept only robust cost improvement.
   if newcp[0]+F(17,100)*newcp[1] >=cp[0]+F(17,100)*cp[1]:continue
   nr,np=rank_pivots(vectors(cols,rr,candidate))
   if nr==len(rr):best,cp,piv=candidate,newcp,np;changes+=1
  shifted=sum(a!=r[2] for r,a in zip(rr,best))
  results.append({'m':m,'H':H,'weights':ws,'K':K,'theta':str(theta),'w0':str(w0),'rows':rr,'columns':cols,'assignment':best,'rank':rank,
                  'minor_column_indices':piv,'normalized_linear':str(cp[0]),'normalized_collision':str(cp[1]),
                  'shifted_rows':shifted,'accepted_changes':changes})
  print('result',rank,shifted,cp,flush=True)
 Path(__file__).with_name('experiment.json').write_text(json.dumps(results,indent=2))

if __name__=='__main__':main()
