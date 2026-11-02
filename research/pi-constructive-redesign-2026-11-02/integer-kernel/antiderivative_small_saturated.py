"""Integer antiderivative coordinates; no pi or target rational enters LLL."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, gcd, lcm
from pathlib import Path
import json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
import sympy.polys.matrices.lll as lll_module
# SymPy 1.14 PythonMPQ math.floor rounds through float for large values.
# Keep its LLL algorithm but supply exact nearest-integer floor.
lll_module.mfloor=lambda value: value.numerator//value.denominator


def data(c,a,b):
 # N(t) coefficient matrix, t=x^2; last coordinate is alpha.
 M=s.zeros(c+1,c+1)
 for j in range(c):
  M[j,j]=2*j+1
  M[j+1,j]=2*j-2*c+3
 for j in range(c):M[j,c]=comb(c-1,j)
 rows=[list(M[j,:]) for j in range(a)]
 rows += [[sum(comb(j,k)*M[j,v] for j in range(k,c+1)) for v in range(c+1)] for k in range(b)]
 C=s.Matrix(rows)
 normal,left,right=smith_normal_decomp(DomainMatrix.from_Matrix(C).convert_to(s.ZZ))
 T=right.to_Matrix();rank=a+b
 assert abs(T.det())==1
 B=T[:,rank:].T
 B=B.lll(delta=F(3,4))
 assert C*B.T==s.zeros(rank,c+1-rank)
 return M,B


def samples(c):
 out=[]
 # Rational samples are selection features only, not error certificates.
 for i in range(1,10):
  t=F(i*i,100);den=(1+t)**c
  out.append([((2*j+1)*t**j+(2*j-2*c+3)*t**(j+1))/den for j in range(c)] + [1/(1+t)])
 return out


def atan(a,K=120):
 v=sum((F((-1)**k,(2*k+1)*a**(2*k+1)) for k in range(K)),F(0))
 return v,v+F(1,(2*K+1)*a**(2*K+1))


def run():
 start=time.monotonic();l,h=atan(5);u,v=atan(239);lo,hi=16*l-4*v,16*h-4*u
 records=[];rows=zero=lll_failures=0;best=F(10**100);seen=set();per_c={};growth={}
 for c in range(4,13):
  S=samples(c)
  for d in (1,2,3):
   for a in sorted(set(((c-d)//2,max(1,(c-d)//3)))):
    b=c-d-a
    if a<1 or b<1:continue
    M,B=data(c,a,b)
    for budget in range(0,4*c+1,8):
     cost=F(1,2**budget)
     feats=S+[[cost if j==c else F(0) for j in range(c+1)], [cost if j<c else F(0) for j in range(c+1)]]
     D=lcm(*(v.denominator for row in feats for v in row))
     feature=s.Matrix([[int(v*D) for v in row] for row in feats]).T
     # identity tail recovers exact coordinates; no numerical relation to pi.
     extended=(B*feature).row_join(B)
     try:
      reduced=extended.lll(delta=F(3,4))
     except AssertionError:
      lll_failures+=1
      continue
     for row in reduced.tolist():
      rows+=1;R=list(map(int,row[-(c+1):]));alpha=R[-1];r1=sum(R[:-1])
      assert all(v==0 for v in (s.Matrix([list(M[j,:]) for j in range(a)]+[[sum(comb(j,k)*M[j,v] for j in range(k,c+1)) for v in range(c+1)] for k in range(b)])*s.Matrix(R)))
      if not alpha:zero+=1;continue
      rat=F(-4*r1,alpha*2**(c-1));p,q=rat.numerator,rat.denominator
      if p<=0 or (p,q) in seen:continue
      seen.add((p,q))
      if rat<lo:L,H=(lo-rat)*q*q,(hi-rat)*q*q
      elif rat>hi:L,H=(rat-hi)*q*q,(rat-lo)*q*q
      else:raise AssertionError('Need narrower pi enclosure')
      if c not in per_c or H<per_c[c][0]:per_c[c]=(H,{'c':c,'a':a,'b':b,'d':d,'p':str(p),'q':str(q),'Q_display':float(H)})
      if H<1 and (c not in growth or q>int(growth[c]['q'])):growth[c]={'c':c,'a':a,'b':b,'d':d,'p':str(p),'q':str(q),'Q_display':float(H)}
      if H<best:
       best=H
       rec={'c':c,'a':a,'b':b,'remaining_degree':d,'cost_weight_power_of_two':-budget,'alpha':str(alpha),'R_coefficients':[str(v) for v in R[:-1]],'R_at_one':str(r1),'p':str(p),'q':str(q),'quality_display':float(H),'max_R_coefficient_bits':max(abs(v).bit_length() for v in R[:-1])}
       scale=10**25;rec['quality_rational_bounds']=[str(F(L.numerator*scale//L.denominator,scale)),str(F(-((-H.numerator*scale)//H.denominator),scale))]
       records.append(rec)
 result={'status':'PASS','selection':'full saturated integer nullspace by Smith decomposition; exact LLL with rational sampling features and alpha/boundary cost; pi used only after selection to evaluate quality','coverage':{'c':[4,12],'remaining_degree':[1,2,3],'sample_x':['1/10','2/10','3/10','4/10','5/10','6/10','7/10','8/10','9/10']},'lll_failures_skipped':lll_failures,'lll_rows_examined':rows,'zero_alpha_rows':zero,'distinct_positive_rationals':len(seen),'elapsed_seconds':time.monotonic()-start,'records':records,'per_c_best':[per_c[k][1] for k in sorted(per_c)],'per_c_largest_q_below_one':[growth[k] for k in sorted(growth)],'scope':'Finite construction experiment; no new rational approximation record or asymptotic theorem. Smith right transform is unimodular, so the exact full integer kernel is used.'}
 Path(__file__).with_name('small_saturated_results.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':run()
