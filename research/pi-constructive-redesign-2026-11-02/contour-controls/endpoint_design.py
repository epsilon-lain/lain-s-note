"""Arithmetic endpoint design for a pi-only paired-pole contour.

s,t>0 generate integer poles 0,L and endpoint u+iv with argument
of z/(L-z) exactly pi/4. Uses no target convergent or continued fractions.
"""
from symmetric_contour import mul,power,quotient,imag
from fractions import Fraction as F
from math import comb,gcd,lcm,log,ceil
from pathlib import Path
import json

def contour(a,b,c,s,t):
 L=2*s*s+2*s*t+t*t;u=s*(2*s+t);v=s*t
 n1=u*u+v*v;n2=(L-u)**2+v*v
 assert n1==2*s*s*L and n2==t*t*L
 # Re(z/(L-z))=Im(z/(L-z)) and both positive.
 assert u*(L-u)-v*v==v*L
 Q0=mul([n1,-2*u,1],[n2,-2*(L-u),1])
 N=mul(power(Q0,a),power([-L,2],2*b))
 den=[0]*(c+1)+power([L,-1],c+1)
 P=quotient(N,den)
 r=[]
 for j in range(c+1):
  k=c-j
  r.append(sum((F(N[i]*comb(c+k-i,k-i),L**(c+1+k-i)) for i in range(min(k,len(N)-1)+1)),F(0)))
 A=r[0]/2
 B=sum((F(p*2*imag(u,v,k+1),k+1) for k,p in enumerate(P)),F(0))
 for j in range(1,c+1):
  B+=F(2,j)*r[j]*(F(imag(L-u,v,j),n2**j)-F(imag(u,-v,j),n1**j))
 M=max(c,len(P));ell=1
 for k in range(1,M+1):ell=lcm(ell,k)
 universal=2*L**(2*c+1)*lcm(2*s*s,t*t)**c*ell
 assert (universal*A).denominator==(universal*B).denominator==1
 return A,B,universal,(L,u,v),N,den

def atan(a,N):
 p=sum((F((-1)**k,(2*k+1)*a**(2*k+1)) for k in range(N)),F(0));tail=F(1,(2*N+1)*a**(2*N+1))
 return (p,p+tail) if N%2==0 else (p-tail,p)
def pibounds(N):
 l,h=atan(5,N);u,v=atan(239,N);return 16*l-4*v,16*h-4*u

def main():
 pl,ph=pibounds(240)
 records=[];rows=[];best=None
 for budget in range(1,11):
  for s in range(1,5):
   for t in range(1,5):
    if gcd(s,t)!=1:continue
    for a in range(1,4):
     for b in range(4):
      c=budget-a-b
      if not 0<=c<=6:continue
      A,B,univ,coords,*_=contour(a,b,c,s,t)
      if not A:continue
      r=-B/A
      if r<=0:continue
      p,q=r.numerator,r.denominator
      if r<pl:lo,hi=(pl-r)*q*q,(ph-r)*q*q
      elif r>ph:lo,hi=(r-ph)*q*q,(r-pl)*q*q
      else:raise AssertionError('pi interval too wide')
      assert lo>0
      D=lcm(A.denominator,B.denominator);g=gcd(int(D*A),int(D*B))
      scale=10**18
      row={'s':s,'t':t,'a':a,'b':b,'c':c,'p':str(p),'q':str(q),'A':str(A),'B':str(B),'geometry_L_u_v':coords,'actual_D_bits':D.bit_length(),'universal_D_bits':univ.bit_length(),'g_bits':g.bit_length(),'quality_display':float((lo+hi)/2),'quality_strict_outward_bounds':[str(F(lo.numerator*scale//lo.denominator,scale)),str(F(-((-hi.numerator*scale)//hi.denominator),scale))]}
      rows.append(row)
      if best is None or hi<best:best=hi;records.append(row)
 result={'status':'PASS','scope':'Finite geometry-and-exponent redesign; known rational approximants are benchmarks, not new pi records. No quality-zero sequence theorem.', 'candidate_count':len(rows),'bounds':{'s_t':[1,4],'gcd_s_t':1,'a':[1,3],'b':[0,3],'c':[0,6],'a_b_c_sum':[1,10]},'record_improvements':records,'best_10':sorted(rows,key=lambda r:r['quality_display'])[:10]}
 Path(__file__).with_name('endpoint_checks.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'status':'PASS','candidate_count':len(rows),'records':records},indent=2))
if __name__=='__main__':main()
