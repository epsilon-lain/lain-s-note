from fractions import Fraction as F
from functools import lru_cache
from math import comb,factorial
from pathlib import Path
from random import Random
import json
import sympy as s
@lru_cache(None)
def J(k,c):
 if c==0:return F(0),F(1,2*k+1)
 if k:
  a,b=J(k-1,c-1);u,v=J(k-1,c);return a-u,b-v
 if c==1:return F(1,4),F(0)
 a,b=J(0,c-1);u=F(2*c-3,2*c-2);return u*a,u*b+F(1,(2*c-2)*2**(c-1))
def matrix(c):
 M=s.zeros(c+1,c+1)
 for j in range(c):M[j,j]=2*j+1;M[j+1,j]=2*j-2*c+3
 for j in range(c):M[j,c]=comb(c-1,j)
 return M
rng=Random(20261102);tests=0
for c in range(1,19):
 M=matrix(c);assert abs(M.det())==2**(c-1)*factorial(c-1)
 for rep in range(20):
  r=[rng.randrange(-1000,1001) for j in range(c+1)]
  N=list(map(int,M*s.Matrix(r)));A=sum((v*J(j,c)[0] for j,v in enumerate(N)),F(0));B=sum((v*J(j,c)[1] for j,v in enumerate(N)),F(0))
  assert A==F(r[-1],4) and B==F(sum(r[:-1]),2**(c-1));tests+=1
result={'status':'PASS','independent_moment_checks':tests,'determinant_degrees':[1,18],'seed':20261102,'scope':'Exact full-parameter algebraic identity and arithmetic-index audits; not an infinite approximation theorem.'}
Path(__file__).with_name('family_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
