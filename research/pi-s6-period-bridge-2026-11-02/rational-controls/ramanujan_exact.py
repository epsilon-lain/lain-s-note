"""A fully rational elliptic-period bridge, with reduced-denominator accounting.

The identity 16/pi = sum (42n+5)*binom(2n,n)^3/4096^n is the only
external analytic input (Ramanujan; Guillera 2018, formula (2)).
No decimal value of pi is used for any error certificate.
"""
from fractions import Fraction as F
from math import comb,gcd,log
from pathlib import Path
import json
import sympy as sy

def term(n):return F((42*n+5)*comb(2*n,n)**3,4096**n)
def pi_interval(N=150):
 s=sum((term(n) for n in range(N+1)),F(0));a=term(N+1)
 return F(16)/(s+F(32,31)*a),F(16)/(s+a)
def quality(r,lo,hi):
 assert not lo<=r<=hi
 p,q=r.numerator,r.denominator
 lower,upper=(lo-r,hi-r) if r<lo else (r-hi,r-lo)
 return lower*q*q,upper*q*q
def decimal_interval(lo,hi,digits=20):
 d=10**digits
 return [str(F(lo.numerator*d//lo.denominator,d)),str(F(-((-hi.numerator*d)//hi.denominator),d))]
def display(x):
 return float(x) if x<10**300 else "10^"+str((log(x.numerator)-log(x.denominator))/log(10))
def main():
 n=sy.symbols('n');difference=sy.expand(16*(42*n+5)*(n+1)**3-(42*n+47)*(2*n+1)**3)
 assert all(c>0 for c in sy.Poly(difference,n).all_coeffs())
 trunc=[];S=F(0);C=0
 for N in range(0,51):
  S+=term(N);C=4096*C+(42*N+5)*comb(2*N,N)**3
  assert S==F(C,4096**N)
  r=F(16)/S;p,q=r.numerator,r.denominator
  s2=N.bit_count();power=3*s2
  assert (C&-C)==2**power
  assert gcd(16*4096**N,C)==2**power
  assert p==2**(12*N+4-power) and q==C//2**power
  a=term(N+1)
  lo=16*a/(S*(S+a))*q*q
  upper=F(32,31)*a
  hi=16*upper/(S*(S+upper))*q*q
  assert lo>0
  if N in [0,1,2,3,5,10,20,50]:
   trunc.append({'N':N,'p':str(p),'q':str(q),'q_digits':len(str(q)),
    'exact_gcd_power_of_two':power,'quality_strict_bounds':decimal_interval(lo,hi),
    'quality_display':display((lo+hi)/2)})
 # Direct diagonal Pade of the integer-coefficient period series H(t).
 # Parameters and coefficient equations do not include pi or a target fraction.
 pl,ph=pi_interval(220);pade=[]
 coeff=[sy.Integer((42*k+5)*comb(2*k,k)**3) for k in range(26)]
 for n in range(1,11):
  matrix=sy.Matrix([[coeff[k-j] for j in range(1,n+1)] for k in range(n+1,2*n+1)])
  rhs=sy.Matrix([-coeff[k] for k in range(n+1,2*n+1)])
  qpoly=[sy.Integer(1)]+list(matrix.inv()*rhs)
  ppoly=[sum(qpoly[j]*coeff[k-j] for j in range(k+1)) for k in range(n+1)]
  for k in range(n+1,2*n+1):assert sum(qpoly[j]*coeff[k-j] for j in range(n+1))==0
  t=sy.Rational(1,4096)
  pval=sum(a*t**j for j,a in enumerate(ppoly));qval=sum(a*t**j for j,a in enumerate(qpoly))
  assert pval and qval
  r0=sy.cancel(16*qval/pval);p,q=map(int,r0.as_numer_denom());r=F(p,q)
  lo,hi=quality(r,pl,ph)
  pade.append({'order':[n,n],'series_coefficients_used':2*n+1,
   'p':str(p),'q':str(q),'q_digits':len(str(q)),
   'quality_strict_bounds':decimal_interval(lo,hi),'quality_display':display((lo+hi)/2),
   'source_polynomial_degrees':[len(ppoly)-1,len(qpoly)-1]})
  # Fair comparison: both consume exactly coefficients 0,...,2n.
  rawS=sum((term(k) for k in range(2*n+1)),F(0));rawr=16/rawS;rawa=term(2*n+1);rawq=rawr.denominator
  rawlo=16*rawa/(rawS*(rawS+rawa))*rawq*rawq
  rawtail=F(32,31)*rawa;rawhi=16*rawtail/(rawS*(rawS+rawtail))*rawq*rawq
  pade[-1].update({'same_coefficients_truncation':{'N':2*n,'p':str(rawr.numerator),'q':str(rawq),'q_digits':len(str(rawq)),'quality_display':display((rawlo+rawhi)/2)},'denominator_reduction_factor':str(F(rawq,q)),'quality_improvement_factor_lower':str(F((rawlo/hi*10**12).numerator//(rawlo/hi*10**12).denominator,10**12)),'quality_improvement_display':display(rawlo/hi)})
 out={'status':'PASS','external_identity':'16/pi=sum_(n>=0) (42n+5)*binom(2n,n)^3/4096^n',
  'primary_source':'https://arxiv.org/pdf/1807.07394, formula (2)',
  'error_certification':'Positive rational series tail; no decimal pi input',
  'uniform_term_ratio_bound':'1/32','ratio_certificate_positive_polynomial':str(difference),
  'truncation_gcd_exact':'2^(3*popcount(N))',
  'truncation_scope':[0,50],'truncation_samples':trunc,
  'pade_scope':[1,10],'pade_candidates':pade,
  'best_pade':min(pade,key=lambda r:F(r['quality_strict_bounds'][1])),
  'scope':'Constructive period-to-pi bridge and finite arithmetic controls. Known identity, not a new pi record or quality-zero theorem.'}
 Path(__file__).with_name('ramanujan_checks.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k!='pade_candidates'},indent=2))
if __name__=='__main__':main()
