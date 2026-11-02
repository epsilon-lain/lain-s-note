"""Exact constructive control: Salikhov-type paired poles and complex endpoints.
No floating pi, no continued fractions, and no target convergent supplied.
"""
from fractions import Fraction as F
from math import comb,gcd,lcm
from pathlib import Path
import json

def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

def power(p,n):
 out=[1]
 for _ in range(n):out=mul(out,p)
 return out

def quotient(n,d):
 n=n[:];q=[0]*max(1,len(n)-len(d)+1)
 while len(n)>=len(d):
  k=len(n)-len(d);assert n[-1]%d[-1]==0
  v=n[-1]//d[-1];q[k]=v
  for j,x in enumerate(d):n[k+j]-=v*x
  while n and not n[-1]:n.pop()
 return q

def imag(a,b,n):
 r,s=1,0
 for _ in range(n):r,s=r*a-s*b,r*b+s*a
 return s

def integral(a,b,c):
 # Q0=((x-4)^2+4)((x-6)^2+4)
 N=mul(power(mul([20,-8,1],[40,-12,1]),a),power([-5,1],2*b))
 D=[0]*(c+1)+power([10,-1],c+1)
 P=quotient(N,D)
 residues=[]
 for j in range(c+1):
  k=c-j
  residues.append(sum((F(N[t]*comb(c+k-t,k-t),10**(c+1+k-t)) for t in range(min(k,len(N)-1)+1)),F(0)))
 A=residues[0]/2
 B=sum((F(v*2*imag(4,2,k+1),k+1) for k,v in enumerate(P)),F(0))
 for j in range(1,c+1):
  B+=F(2,j)*residues[j]*(F(imag(6,2,j),40**j)-F(imag(4,-2,j),20**j))
 M=max(c,len(P))
 L=1
 for k in range(1,M+1):L=lcm(L,k)
 universal=2*10**(2*c+1)*4**c*L
 assert (universal*A).denominator==(universal*B).denominator==1
 return A,B,N,D,P,residues,universal

def atan(a):
 s=sum((F((-1)**k,(2*k+1)*a**(2*k+1)) for k in range(120)),F(0))
 return s,s+F(1,241*a**241)

def main():
 l,h=atan(5);u,v=atan(239);pl,ph=16*l-4*v,16*h-4*u
 rows=[];best=None;records=[]
 for cost in range(1,17):
  for a in range(1,5):
   for b in range(5):
    c=cost-a-b
    if not 0<=c<=8:continue
    A,B,*_,universal=integral(a,b,c)
    if not A:continue
    r=-B/A
    if r<=0:continue
    p,q=r.numerator,r.denominator
    if r<pl:lo,hi=(pl-r)*q*q,(ph-r)*q*q
    elif r>ph:lo,hi=(r-ph)*q*q,(r-pl)*q*q
    else:raise AssertionError('pi interval too wide')
    scale=10**18
    row={'a':a,'b':b,'c':c,'A':str(A),'B':str(B),'p':str(p),'q':str(q),
         'quality_strict_outward_bounds':[str(F(lo.numerator*scale//lo.denominator,scale)),str(F(-((-hi.numerator*scale)//hi.denominator),scale))],
         'quality_display':float((lo+hi)/2),'universal_denominator_bits':universal.bit_length(),
         'actual_joint_denominator_bits':lcm(A.denominator,B.denominator).bit_length()}
    rows.append(row)
    if best is None or hi<best:best=hi;records.append(row)
 result={'status':'PASS','candidate_count':len(rows),'scope':'Finite constructive experiment only, not a new approximation record or infinite quality-zero family.', 'parameter_bounds':{'a':[1,4],'b':[0,4],'c':[0,8],'a_b_c_sum':[1,16]},'record_improvements':records,'best_10':sorted(rows,key=lambda r:r['quality_display'])[:10]}
 Path(__file__).with_name('symmetric_contour_checks.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
