from fractions import Fraction as F
from math import comb,factorial,isqrt
from pathlib import Path
from functools import lru_cache
import json,sys
import matrix_core as ex
import search

P=1_000_000_009
I=430_477_711

def det(matrix):
 a=[r[:] for r in matrix];n=len(a);d=1
 assert all(len(r)==n for r in a)
 for c in range(n):
  k=next((r for r in range(c,n) if a[r][c]),None)
  if k is None:return 0
  if k!=c:a[k],a[c]=a[c],a[k];d=-d
  v=a[c][c];d=d*v%P;inv=pow(v,-1,P)
  for r in range(c+1,n):
   f=a[r][c]*inv%P
   if f:
    for j in range(c+1,n):a[r][j]=(a[r][j]-f*a[c][j])%P
 return d%P

def original_matrix(columns,rows,centers):
 N=max(r[1] for r in rows)
 maxpow=max(sum(a) for h,a in columns)
 log=[0]+[((-1)**(k+1)*pow(k,-1,P))%P for k in range(1,N+1)]
 powcache={}
 for j in set(r[0] for r in rows):
  for i,R in enumerate(centers):
   f=log[:];f[0]=j*R%P
   powcache[j,i]=ex.powers(f,maxpow,N)
 @lru_cache(None)
 def poly(j,degree):
  f=[1]+[0]*N
  for i,d in enumerate(degree):f=ex.mul(f,powcache[j,i][d],N)
  return f
 matrix=[]
 for j,s,beta in rows:
  line=[]
  for h,alpha in columns:
   if any(b>a for b,a in zip(beta,alpha)):line.append(0);continue
   factor=1
   for a,b in zip(alpha,beta):factor=factor*comb(a,b)%P
   f=poly(j,tuple(a-b for a,b in zip(alpha,beta)))
   line.append(factor*sum(comb(h,k)*f[s-k] for k in range(min(s,h)+1))%P)
  matrix.append(line)
 return matrix

def exp_bounds(x,N=60):
 x=F(x)
 lo=sum((x**k/factorial(k) for k in range(N+1)),F())
 return lo,lo+x**(N+1)/factorial(N+1)/(1-x/(N+2))

def atan_bounds(d,N):
 s=sum((F((-1)**k,(2*k+1)*d**(2*k+1)) for k in range(N)),F())
 t=F((-1)**N,(2*N+1)*d**(2*N+1))
 return min(s,s+t),max(s,s+t)

def main():
 assert all(P%d for d in range(2,isqrt(P)+1))
 assert I*I%P==P-1
 ex.PRIME=P
 al,au=atan_bounds(5,30);bl,bu=atan_bounds(239,10)
 pil,piu=16*al-4*bu,16*au-4*bl
 approximants={1:(6,2),2:(22,7),5:(355,113)}
 for weight,(p,q) in approximants.items():
  assert exp_bounds(weight-1)[1]<q<exp_bounds(weight)[0]
  assert max(abs(pil-F(p,q)),abs(piu-F(p,q)))<F(1,q*q)
 fixtures=json.loads(Path(__file__).with_name('fixtures.json').read_text())
 results=[]
 for x in fixtures:
  H=x['H'];weights=x['weights'];K=x['K'];theta=F(x['theta']);v0=F(x['v0']);w0=F(x['w0'])
  m=len(weights)
  assert v0<w0 and K*theta**m<1 and K*w0/v0*theta**m<1
  rr=search.rows(H,weights,K,theta,v0)
  cc=ex.columns(m,H,weights,w0)
  assert len(rr)==x['rows'] and len(cc)==x['columns']
  assert len(set(rr))==len(rr) and len(set(cc))==len(cc)
  assert all(v0*s+sum(w*b for w,b in zip(weights,beta))/theta<H for j,s,beta in rr)
  assert all(w0*h+sum(w*a for w,a in zip(weights,alpha))<=H for h,alpha in cc)
  aa=[tuple(a) for a in x['assignment']];ix=x['minor_indices']
  assert rr==search.rows(H,weights,K,theta,v0)
  assert cc==ex.columns(m,H,weights,w0)
  assert len(ix)==len(rr)==len(aa) and len(set(ix))==len(ix)
  assert all(0<=i<len(cc) for i in ix)
  assert all(all(a>=b for a,b in zip(a,r[2])) and sum(w*v for w,v in zip(weights,a))<=H for r,a in zip(rr,aa))
  assert all(a==r[2] for r,a in zip(rr,aa) if r[0]==0)
  assert max(r[1] for r in rr)<34
  test=ex.vectors(cc,rr,aa,W=2)
  td=det([[row[i] for i in ix] for row in test]);assert td
  centers=[2*I*p*pow(q,-1,P)%P for p,q in (approximants[w] for w in weights)]
  orig=original_matrix(cc,rr,centers)
  od=det([[row[i] for i in ix] for row in orig]);assert od
  base=ex.cost_parts(rr,[r[2] for r in rr],H,weights)
  cost=ex.cost_parts(rr,aa,H,weights)
  assert cost[1]<base[1] and cost[0]-base[0]+F(1,6)*(cost[1]-base[1])<0
  bbar=F(sum(sum(w*b for w,b in zip(weights,r[2])) for r in rr),len(rr)*H)
  assert cost[0]+cost[1]/4<1-bbar
  result={'m':m,'H':H,'weights':weights,'K':K,'theta':str(theta),'v0':str(v0),'w0':str(w0),
          'rows':len(rr),'columns':len(cc),'shifted_rows':sum(a!=r[2] for r,a in zip(rr,aa)),
          'positive_type_steps':x['accepted_moves'],'base_cost':[str(v) for v in base],
          'new_cost':[str(v) for v in cost],'bbar':str(bbar),
          'translated_det_at_W2_mod_prime':td,'original_det_same_columns_mod_prime':od,
          'cost_strictly_decreases':True,'cost_below_1_minus_bbar':True}
  results.append(result)
  print(result,flush=True)
 # Capacity failure: 40 a=0 rows but only 36 different functions e^(hz)z^d.
 H=10;weights=(2,5);theta=F(1,4);v0=F(1,2);w0=F(1);K=2
 rr=search.rows(H,weights,K,theta,v0);cc=ex.columns(2,H,weights,w0)
 pure=sum(beta==(0,0) for j,s,beta in rr)
 function_dim=sum(int((H-min(weights)*d)/w0)+1 for d in range(H//min(weights)+1))
 matrix=ex.vectors(cc,rr,[r[2] for r in rr],W=2)
 rank,piv=ex.rank_pivots(matrix)
 assert pure==40 and function_dim==36 and rank==44==len(rr)-(pure-function_dim)
 report={'status':'passed','modulus':P,'sqrt_minus_one':I,'positive_certificates':results,
         'unshifted_capacity_failure':{'m':2,'H':H,'weights':weights,'theta':str(theta),'v0':str(v0),'w0':str(w0),'K':K,
             'rows':len(rr),'columns':len(cc),'pure_a0_rows':pure,'a0_function_dimension':function_dim,
             'rank_lower_bound_mod_prime':rank,'rank_upper_bound_analytic':44,'exact_rank_at_period':44},
         'scope':'Finite structural certificates. Actual distinct rational centers and nu=2 error hypotheses verified. Unknown separated-weight thresholds and all-large-H conclusions are not certified.'}
 Path(__file__).with_name('verification.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
