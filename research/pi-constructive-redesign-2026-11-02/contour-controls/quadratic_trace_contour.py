"""Critical algebraic endpoints: normalize the residue, then take Q(sqrt2) trace.
This is a bounded constructive experiment, not an infinite approximation claim.
"""
from fractions import Fraction as F
from math import comb,gcd,lcm,log,ceil
from pathlib import Path
import json
import sympy as sy
from endpoint_design import pibounds
K=sy.QQ.algebraic_field(sy.sqrt(2));zero=K.zero;one=K.one
v=K.from_sympy(sy.sqrt(2)-1)

def mul(a,b):
 c=[zero]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

def power(p,n):
 out=[one]
 for _ in range(n):out=mul(out,p)
 return out

def quotient(n,d):
 n=n[:];q=[zero]*max(1,len(n)-len(d)+1)
 while len(n)>=len(d):
  k=len(n)-len(d);value=n[-1]/d[-1];q[k]=value
  for j,x in enumerate(d):n[k+j]-=value*x
  while n and not n[-1]:n.pop()
 return q

def im_power(n):
 r,s=one,zero
 for _ in range(n):r,s=r-s*v,r*v+s
 return s

def rational(x):return F(int(x.numerator),int(x.denominator))
def trace(x):
 c=x.to_list()
 return 2*rational(c[-1]) if c else F(0)

def normalized_trace(a,b,c):
 N=mul(power([one+v*v,-2*one,one],2*a),power([-one,one],2*b))
 return from_numerator(N,c)

def from_numerator(N,c):
 D=[zero]*(c+1)+power([2*one,-one],c+1)
 P=quotient(N,D)
 residues=[]
 for j in range(c+1):
  degree=c-j
  residues.append(sum((N[k]*K.convert(F(comb(c+degree-k,degree-k),2**(c+1+degree-k))) for k in range(min(degree,len(N)-1)+1)),zero))
 r0=residues[0]
 if not r0:return None
 B=sum((coef*2*im_power(k+1)/(k+1) for k,coef in enumerate(P)),zero)
 for j in range(1,c+1):B+=4*residues[j]/j*im_power(j)/(one+v*v)**j
 T=trace(B/r0)
 return T,r0,B

def main():
 pl,ph=pibounds(300)
 rows=[];records=[];best=None
 for budget in range(1,15):
  for a in range(1,5):
   for b in range(5):
    c=budget-a-b
    if not 0<=c<=10:continue
    result=normalized_trace(a,b,c)
    if result is None:continue
    r,r0,B=result
    if r<=0:continue
    p,q=r.numerator,r.denominator
    if r<pl:lo,hi=(pl-r)*q*q,(ph-r)*q*q
    elif r>ph:lo,hi=(r-ph)*q*q,(r-pl)*q*q
    else:raise AssertionError('pi interval too wide')
    scale=10**18
    row={'a':a,'b':b,'c':c,'p':str(p),'q':str(q),'residue':str(K.to_sympy(r0)),'rational_part':str(K.to_sympy(B)),
         'quality_display':float((lo+hi)/2),'quality_strict_outward_bounds':[str(F(lo.numerator*scale//lo.denominator,scale)),str(F(-((-hi.numerator*scale)//hi.denominator),scale))]}
    rows.append(row)
    if best is None or hi<best:best=hi;records.append(row)
 result={'status':'PASS','candidate_count':len(rows),'field':'Q(sqrt2)','geometry':'poles0,2; endpoints1±i(sqrt2−1); conjugate endpoints use signed v\'=-sqrt2−1', 'identity':'r=Tr(B/r0); (I/r0)+sigma(I/r0)=-pi+r','scope':'Finite normalized trace experiment; no new pi record or infinite quality-zero theorem.', 'record_improvements':records,'best_10':sorted(rows,key=lambda r:r['quality_display'])[:10]}
 Path(__file__).with_name('quadratic_trace_checks.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'status':'PASS','count':len(rows),'records':records},indent=2))
if __name__=='__main__':main()
