"""Independent finite Smith/Hermite checks of the proved odd period image."""
from math import comb,lcm
from pathlib import Path
import json
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
from sympy.matrices.normalforms import hermite_normal_form

def odd(n):
 while n%2==0:n//=2
 return n

def image(c,a,b):
 M=s.zeros(c+1,c+1)
 for j in range(c):M[j,j]=2*j+1;M[j+1,j]=2*j-2*c+3
 for j in range(c):M[j,c]=comb(c-1,j)
 rows=[list(M[j,:]) for j in range(a)]
 rows += [[sum(comb(j,k)*M[j,v] for j in range(k,c+1)) for v in range(c+1)] for k in range(b)]
 C=s.Matrix(rows);D,U,V=smith_normal_decomp(DomainMatrix.from_Matrix(C).convert_to(s.ZZ));V=V.to_Matrix()
 assert abs(V.det())==1
 B=V[:,a+b:].T
 assert C*B.T==s.zeros(a+b,c+1-a-b)
 P=s.Matrix([[int(v[-1]),sum(map(int,v[:-1]))] for v in B.tolist()]).T
 return hermite_normal_form(P)
checks=0
for total in range(2,11):
 for a in range(1,total):
  b=total-a
  for c in [total+1,total+2]:
   H=image(c,a,b);L=lcm(*range(1,max(2*a-1,b)+1));assert odd(int(H.det()))==odd(L);checks+=1
fair=[]
for a,b in [(7,7),(5,9)]:
 H=image(16,a,b);expected=s.diag(45045,1024) if a==7 else s.diag(315,32768);assert H==expected
 fair.append({'c':16,'a':a,'b':b,'H':[[str(v) for v in row] for row in H.tolist()],'index':str(H.det())})
assert int(fair[0]['index'])*32==int(fair[1]['index'])*143
result={'status':'PASS','odd_period_lattice_checks':checks,'same_c':16,'same_a_plus_b':14,'period_lattices':fair,'odd_gain_factor':143,'dyadic_cost_factor':32,'complete_covolume_gain':'143/32','scope':'Exact arithmetic period-image improvement; no automatic Q improvement.'}
Path(__file__).with_name('fair_endpoint_comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
