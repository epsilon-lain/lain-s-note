"""Exact audits of the actual-initial differential Bezout witness."""
from fractions import Fraction as F
from math import gcd, log, sqrt
from pathlib import Path
import argparse
import json

def add(a,b):
 c=[F(0)]*max(len(a),len(b))
 for j,x in enumerate(a):c[j]+=x
 for j,x in enumerate(b):c[j]+=x
 while len(c)>1 and not c[-1]:c.pop()
 return c

def scale(a,s):return [s*x for x in a]
def mul(a,b):
 c=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 while len(c)>1 and not c[-1]:c.pop()
 return c
def derivative(a):return [j*x for j,x in enumerate(a)][1:] or [F(0)]

def vp(x,p):
 if not x:raise ValueError('valuation of zero is not used')
 e=0
 while x%p==0:e+=1;x//=p
 return e

def floorlog(n,p):
 e=0
 while n>=p:e+=1;n//=p
 return e

def main(N=1000):
 P=[[F(1)],[F(0),F(1)]];W=[[F(0)],[F(1)]]
 for n in range(1,25):
  for seq in (P,W):
   seq.append(scale(add(scale([F(0)]+seq[-1],F(2*n+1)),scale(seq[-2],F(-n))),F(1,n+1)))
 for n in range(26):
  cross=add(mul(derivative(P[n]),W[n]),scale(mul(P[n],derivative(W[n])),F(-1)))
  assert add(mul(P[n],P[n]),mul([F(1),F(0),F(-1)],cross))==[F(1)]
 primes=[p for p in range(3,N+1,2) if all(p%d for d in range(3,int(p**.5)+1,2))]
 A0,A=1,1;B0,B=0,1;fac=1;strict=[];tests=0;excluded=0
 forced_negative=0;prime_power_boundaries=[];removed_events=[]
 for n in range(1,N+1):
  if n>1:
   k=n-1;A0,A=A,(2*k+1)*A+k*k*A0;B0,B=B,(2*k+1)*B+k*k*B0;fac*=n
  assert (A<<n)%fac==0
  U=(A<<n)//fac;C=U//(U&-U);g=gcd(A,4*B)
  previous_factorial=fac//n
  assert A*B0-A0*B==(-1)**n*previous_factorial**2
  excess=g//gcd(g,fac);excess//=excess&-excess
  old=new=1;removed=[]
  for p in primes:
   if p>n:break
   if n==1:continue
   k=floorlog(n-1,p)-vp(n,p)
   m=vp(C,p);vfac=0;x=n
   while x:x//=p;vfac+=x
   if m>k:
    tests+=1
    previous_a_valuation=vp(A0,p)-vp(previous_factorial,p)
    expected=-vp(n,p)-previous_a_valuation
    assert vp(B,p)-vfac==expected
    assert vp(g,p)-vfac==expected<=0
    if expected<0:forced_negative+=1
    if k<0:
     prime_power_boundaries.append({'n':n,'p':p,'m':m,
                                    'k':k,'exact_surplus':expected})
   if vp(g,p)>vfac:assert 1<=m<=k
   cap=max(0,k)
   contribution=p**min(m,cap);old*=contribution
   if 1<=m<=k:new*=p**m
   elif contribution>1:
    removed.append(p);excluded+=1
    removed_events.append({'n':n,'p':p,'m':m,'k':k,
                           'removed_power':str(contribution)})
  assert excess and new%excess==0
  if old>new:
   strict.append({'n':n,'old_over_new':str(old//new),'removed_primes':removed})
 result={'status':'PASS','polynomial_identity_degrees':[0,25],
         'sequence_indices':[1,N],'forced_nonpositive_prime_cases':tests,
         'forced_strict_deficit_cases':forced_negative,
         'prime_power_boundary_cases':len(prime_power_boundaries),
         'first_prime_power_boundaries':prime_power_boundaries[:12],
         'strict_gain_sieve_improvements':len(strict),'removed_prime_events':excluded,
         'first_improvements':strict[:20],
         'first_removed_events':removed_events[:20],
         'simple_prime_mass_threshold':(2+log(2)-log(2+sqrt(2)))/3,
         'scope':'Exact finite audits of general proof; no prime-clustering asymptotic asserted.'}
 Path(__file__).with_name('differential_witness_checks.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))

if __name__=='__main__':
 ap=argparse.ArgumentParser()
 ap.add_argument('--N',type=int,default=1000)
 main(ap.parse_args().N)
