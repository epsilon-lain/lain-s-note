"""Independent polynomial/period/quality checks for saved constrained candidates.
Only the Python standard library; no imported generation or certificate implementation.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,gcd
import json
HERE=Path(__file__).resolve().parent

def trim(p):
 p=list(p)
 while len(p)>1 and p[-1]==0:p.pop()
 return p

def mul(p,q):
 out=[0]*(len(p)+len(q)-1)
 for i,a in enumerate(p):
  for j,b in enumerate(q):out[i+j]+=a*b
 return trim(out)

def divrem(p,q):
 work=list(map(F,p));q=list(map(F,trim(q)));out=[F(0)]*max(1,len(work)-len(q)+1)
 for j in range(len(work)-len(q),-1,-1):
  val=work[j+len(q)-1]/q[-1];out[j]=val
  for k,c in enumerate(q):work[j+k]-=val*c
 return trim(out),trim(work[:len(q)-1] or [0])

def atan(q,n):
 x=sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(n)),F(0))
 y=x+F((-1)**n,(2*n+1)*q**(2*n+1));return min(x,y),max(x,y)

def run():
 data=json.loads((HERE/'results.json').read_text());checked=0;methods={}
 a,b=atan(5,100),atan(239,36);PI=(16*a[0]-4*b[1],16*a[1]-4*b[0])
 S=[5,-10,1];Ssq=mul(S,S)
 P=[-1,-1,0,0,0,1,1]
 assert divrem(P,[-1,0,1])[1]==[0] and divrem(P,[-1,0,0,0,0,1])[1]==[0]
 # Pull back P((1+i*x)/(1-i*x))*(1-i*x)^6 without complex floats.
 cyc=[(1,0),(0,1),(-1,0),(0,-1)]
 result=[(0,0) for j in range(7)]
 for j,c in enumerate(P):
  for r in range(j+1):
   for k in range(6-j+1):
    re,im=cyc[(r+k)%4];v=c*comb(j,r)*comb(6-j,k)*(-1)**k
    old=result[r+k];result[r+k]=(old[0]+v*re,old[1]+v*im)
 assert result==[(0,0),(0,20),(0,0),(0,-40),(0,0),(0,4),(0,0)]
 for case in data['rows']:
  m=case['m'];c=case['c'];a,b=m,2*m
  assert case['basis_rank']==(m if case['method']=='core_S5' else m+2)
  pc=case['basis_period_coordinates']
  actual_rank=(0 if all(v==[0,0] for v in pc) else 1)
  if any(a*d-b*c for a,b in pc for c,d in pc):actual_rank=2
  assert actual_rank==case['period_map_rank']
  assert actual_rank==(1 if case['method']=='core_S5' and m==1 else 2)
  for row in case['best']:
   R=list(map(int,row['R']));alpha=int(row['alpha']);N=[0]*(c+1)
   for j,r in enumerate(R):
    N[j]+=(2*j+1)*r;N[j+1]+=(2*j-2*c+3)*r
   for j in range(c):N[j]+=alpha*comb(c-1,j)
   H=list(map(int,row['H']))
   endpoint=[0]*a+[(-1)**j*comb(b,j) for j in range(b+1)]
   assert trim(N)==mul(endpoint,H)
   rat=F(-sum(R),2**(c-3)*alpha)
   assert str(rat.numerator)==row['p'] and str(rat.denominator)==row['q']
   if case['method']=='core_S5':assert divrem(H,S)[1]==[0]
   if case['method']=='literal_shared_g':
    lam=F(case['lambda']);W=list(map(F,[comb(6,j) for j in range(7)]))
    for j,v in enumerate(Ssq):W[j+1]-=4*lam*v
    assert divrem(H,W)[1]==[0]
   energy=sum((F(h*k,2*i+2*j+1) for i,h in enumerate(H) for j,k in enumerate(H)),F(0))
   assert energy==F(row['energy']) and energy>0
   bound=(F(4*rat.denominator**2,abs(alpha)*16**m))**2*energy
   assert bound==F(row['gram_Q_bound_squared'])
   lo,hi=PI[0]-rat,PI[1]-rat;assert lo*hi>0
   qlo,qhi=min(abs(lo),abs(hi))*rat.denominator**2,max(abs(lo),abs(hi))*rat.denominator**2
   stored_lo,stored_hi=map(F,row['quality_interval'])
   assert stored_lo<=qlo<=qhi<=stored_hi and qhi*qhi<bound
   checked+=1;methods[case['method']]=methods.get(case['method'],0)+1
 out={'status':'PASS','independent_standard_library':True,'cayley_gaussian_integer_identity':True,'shared_unit_root_divisions':2,'cases':len(data['rows']),'selected_candidates_replayed':checked,'methods':methods,'scope':'All saved selected candidates independently checked; no claim of global search or a new pi approximation record.'}
 (HERE/'independent_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':run()
