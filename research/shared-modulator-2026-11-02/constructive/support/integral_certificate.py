"""Rational integral error certificates: no value or digits of pi needed."""
from fractions import Fraction as F
from math import comb,gcd
from pathlib import Path
import json

def numerator(c,r,alpha):
 n=[0]*(c+1)
 for j in range(c+1):
  if j<c:n[j]+=(2*j+1)*r[j]
  if j:n[j]+=(2*j-2*c+1)*r[j-1]
  if j<c:n[j]+=alpha*comb(c-1,j)
 return n

def integral_interval(N,c,K=80):
 if c<1 or K<c:raise ValueError('require c>=1 and K>=c')
 poly=list(N);total=F(0)
 for k in range(K+1):
  moment=sum((F(t,2*j+1) for j,t in enumerate(poly)),F(0))
  total+=F((-1)**k*comb(c+k-1,k),3**k)*moment
  nxt=[0]*(len(poly)+1)
  for j,t in enumerate(poly):nxt[j]-=t;nxt[j+1]+=2*t
  poly=nxt
 scale=F(2,3)**c
 total*=scale
 # Every further absolute-series ratio is at most r.
 ratio=F(c+K+1,3*(K+2))
 assert ratio<1
 tail=scale*F(comb(c+K,K+1),3**(K+1))/(1-ratio)
 # Integral |N(x)| <= sum |N_j|/(2j+1).
 mass=sum((F(abs(t),2*j+1) for j,t in enumerate(N)),F(0))
 radius=tail*mass
 return total-radius,total+radius

def certificate(c,r,alpha,K=80):
 if c<3 or not alpha:raise ValueError('require c>=3, alpha nonzero')
 N=numerator(c,r,alpha);fraction=F(-sum(r),2**(c-3)*alpha)
 lo,hi=integral_interval(N,c,K)
 abslo=F(0) if lo<=0<=hi else min(abs(lo),abs(hi))
 abshi=max(abs(lo),abs(hi))
 factor=F(4*fraction.denominator**2,abs(alpha))
 return {'c':c,'R_coefficients':r,'alpha':alpha,'N_coefficients':N,'K':K,
  'p':str(fraction.numerator),'q':str(fraction.denominator),
  'quality_lower':str(factor*abslo),'quality_upper':str(factor*abshi),
  'quality_upper_display':float(factor*abshi),
  'certified_quality_below_one':factor*abshi<1}

def main():
 r=[339,1921,4473,5455,3560,1200,92]
 result=certificate(7,r,-339)
 assert (result['p'],result['q'])==('355','113')
 assert F('0.00340631193013')<F(result['quality_lower'])<F(result['quality_upper'])<F('0.00340631193015')
 out={'status':'PASS','method':'Uniform binomial expansion around x^2=1/2 with a geometric rational tail bound','pi_or_continued_fraction_input':False,'scope':'A known rational is a verification baseline; no new record or infinite quality-zero theorem.','baseline':result}
 Path(__file__).with_name('integral_certificate_checks.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k!='baseline'},indent=2));print({k:result[k] for k in ('p','q','quality_upper_display','certified_quality_below_one')})
if __name__=='__main__':main()
