from fractions import Fraction as F
from pathlib import Path
from math import lcm,gcd
import json
from verify_family import J
from antiderivative_search import atan
r=json.loads(Path(__file__).with_name('new_family_results.json').read_text());l,h=atan(5,200);u,v=atan(239,200);lo,hi=16*l-4*v,16*h-4*u
checks=0
for group in r['per_m']:
 m=group['m'];a,b,c=m,2*m,4*m+1;L=lcm(*range(1,2*m+1));L//=L&-L
 for row in group['candidates']:
  alpha=int(row['alpha']);rv=list(map(int,row['R_coefficients']));N=list(map(int,row['N_coefficients']));H=list(map(int,row['H_coefficients']));q=int(row['q']);p=int(row['p'])
  assert alpha%L==0
  assert sum(rv)==int(row['R_at_one'])
  assert F(p,q)==F(-4*sum(rv),alpha*2**(c-1)) and gcd(p,q)==1
  A=sum((v*J(j,c)[0] for j,v in enumerate(N)),F(0));B=sum((v*J(j,c)[1] for j,v in enumerate(N)),F(0))
  assert A==F(alpha,4) and B==F(sum(rv),2**(c-1))
  errlo,errhi=((lo-F(p,q))*q*q,(hi-F(p,q))*q*q) if F(p,q)<lo else ((F(p,q)-hi)*q*q,(F(p,q)-lo)*q*q)
  blo,bhi=map(F,row['quality_rational_bounds']);assert blo<errlo<errhi<bhi
  shape=F(4*q*q*sum(abs(x) for x in H),abs(alpha))*F(2,27)**m
  assert shape==F(row['shape_quality_upper']) and errhi<shape
  checks+=1
result={'status':'PASS','completed_m':r['completed_m'],'candidate_identity_and_quality_checks':checks,'scope':'Independent moment identities, reduced denominators, mandatory odd factors, shape upper bounds and independently truncated Machin error intervals. No infinite-subsequence theorem.'}
Path(__file__).with_name('new_family_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
