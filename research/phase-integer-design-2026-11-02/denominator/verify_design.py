#!/usr/bin/env python3
"""Independent replay of saved primitive fractions and exact error certificates."""
from fractions import Fraction as F
from math import gcd,lcm,comb,isqrt
from pathlib import Path
import json
import sys

sys.set_int_max_str_digits(0)
root=Path(__file__).parent
saved=json.loads((root/'results.json').read_text())

def partial(m):
    return 4*sum((F(1 if j%2==0 else -1,2*j+1) for j in range(m)),F())

def beta_by_integrating_polynomial(a,b):
    return sum((F((-1)**j*comb(b-1,j),a+j) for j in range(b)),F())

def alternate(k,n):
    value=F()
    for j in range(n):value+=F((-1)**j,(2*j+1)*k**(2*j+1))
    nxt=value+F((-1)**n,(2*n+1)*k**(2*n+1))
    return sorted([value,nxt])

a,b=alternate(5,70),alternate(239,20)
plo,phi=16*a[0]-4*b[1],16*a[1]-4*b[0]
assert [str(plo),str(phi)]==saved['pi_interval']

def err_bounds(value):
    assert value<plo or value>phi
    return min(abs(value-plo),abs(value-phi)),max(abs(value-plo),abs(value-phi))

tower_count=0
for row in saved['prime_power_rows']:
    m,p,r=row['m'],row['prime'],row['exponent']
    assert p>=3 and all(p%d for d in range(2,isqrt(p)+1)) and p**r==m
    R=partial(m);T=R-F(1,m)
    assert str(R.numerator)==row['old_numerator'] and str(R.denominator)==row['old_denominator']
    assert str(T.numerator)==row['new_numerator'] and str(T.denominator)==row['new_denominator']
    assert gcd(T.numerator,T.denominator)==1 and R.denominator%m==0
    factor=3 if p==3 and r%2==0 else 1
    assert R.denominator==factor*T.denominator
    A0,A1,A2=partial(m-1),R,partial(m+1)
    D=lcm(A0.denominator,A1.denominator,A2.denominator)
    aa,bb,cc=[x.numerator*(D//x.denominator) for x in (A0,A1,A2)]
    ratio=(A2-A1)/(A1-A0)
    h=(bb-aa)//ratio.denominator
    assert abs(h).bit_length()==row['predicted_common_factor_bits']
    N=aa*cc-bb*bb;Delta=aa-2*bb+cc
    g=gcd(abs(N),abs(D*Delta))
    assert g%abs(h)==0 and g.bit_length()==row['actual_common_factor_bits']
    assert F(N,D*Delta)==T
    elo,ehi=err_bounds(T);olo,ohi=err_bounds(R)
    lower=F(1,2*m*(m+1)*(2*m+1));upper=F(1,4*m**3)
    assert lower<elo<=ehi<upper
    bound=F(1,factor*factor*(4*m*m-1))
    assert bound==F(row['quality_ratio_strict_upper'])
    assert F(T.denominator**2,R.denominator**2)*ehi/olo<bound
    tower_count+=1

high_count=0
improved=worse=0
for row in saved['high_order_rows']:
    m,k=row['m'],row['k']
    value=partial(m)-F(1,m)
    base=value
    for j in range(1,k+1):
        value+=F((-1)**(j-1),2**(j-1))*beta_by_integrating_polynomial(2*m-j,2*j+1)
    assert value==F(row['fraction']) and str(value.denominator)==row['q']
    assert gcd(int(row['p']),int(row['q']))==1
    lower=F(1,2**k)*beta_by_integrating_polynomial(2*m-k,2*k+3)
    upper=2*lower
    assert [str(lower),str(upper)]==[row['error_lower'],row['error_upper']]
    elo,ehi=err_bounds(value);blo,bhi=err_bounds(base)
    assert lower<elo<=ehi<upper
    ratio0=F(value.denominator**2,base.denominator**2)*elo/bhi
    ratio1=F(value.denominator**2,base.denominator**2)*ehi/blo
    assert [str(ratio0),str(ratio1)]==row['quality_ratio_interval']
    cls='improved' if ratio1<1 else ('worse' if ratio0>1 else 'unresolved')
    assert cls==row['classification']
    improved+=cls=='improved';worse+=cls=='worse'
    high_count+=1

universal_count=0
for m in range(1,1001):
    R=partial(m);T=R-F((-1)**(m-1),m)
    assert lcm(R.denominator,m)%T.denominator==0
    upper=F(1,6) if m==1 else beta_by_integrating_polynomial(2*m-1,3)-beta_by_integrating_polynomial(2*m-1,5)/2
    bound=F(T.denominator**2,R.denominator**2)*upper/(F(1,m)-upper)
    assert bound<F(1,4)
    if m&(m-1)==0:
        row=next(x for x in saved['dyadic_sharpness_rows'] if x['m']==m)
        assert T.denominator==m*R.denominator and bound==F(row['certified_ratio_upper'])
    universal_count+=1

fast=json.loads((root/'fast_channels.json').read_text())
fast_count=0
def atan_partial(k,n):
    return sum((F((-1)**j,(2*j+1)*k**(2*j+1)) for j in range(n)),F())
def direct_aitken(values):
    x,y,z=values
    return (x*z-y*y)/(x-2*y+z)
for row in fast['rows']:
    n=row['n']
    separate=16*direct_aitken([atan_partial(5,n+j) for j in range(3)])-4*direct_aitken([atan_partial(239,n+j) for j in range(3)])
    mixed=direct_aitken([16*atan_partial(5,n+j)-4*atan_partial(239,n+j) for j in range(3)])
    assert str(separate)==row['separate_fraction'] and str(mixed)==row['mixed_fraction']
    se,me=err_bounds(separate),err_bounds(mixed)
    lower=F(mixed.denominator**2,separate.denominator**2)*me[0]/se[1]
    upper=F(mixed.denominator**2,separate.denominator**2)*me[1]/se[0]
    assert [str(lower),str(upper)]==row['quality_gain_interval'] and lower>1
    fast_count+=1

out={'status':'PASS','independent_prime_power_cases':tower_count,
     'independent_fast_channel_cases':fast_count,
     'universal_quarter_bound_checks':universal_count,
     'independent_high_order_cases':high_count,'actual_gcd_and_denominator_replays':tower_count+high_count,
     'high_order_improved':improved,'high_order_worse':worse,
     'scope':'Finite exact replay; the infinite theorem is proved in NOTE.md.'}
(root/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
