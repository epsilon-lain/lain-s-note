#!/usr/bin/env python3
"""Certified finite interfaces to the specified infinite-coordinate kernel.

Parameter denominators label original rational candidates; an infinite
rational-coordinate vector is not assigned a free finite common denominator.
"""
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
from math import sqrt,gcd
import json,sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'rational-scale-geometry-2026-11-01'))
from denominator_envelope_lab import pi_interval


def determinant(A):
    n=len(A);total=F(0)
    for p in permutations(range(n)):
        term=F((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n)))
        for i,j in enumerate(p):term*=A[i][j]
        total+=term
    return total


def distance_ratio(s,t):
    """Full Hilbert squared distance / scalar squared distance; continuous diagonal extension."""
    return (1+s*t/4)/(4*(1-s*s/4)*(1-t*t/4)*(1-s*t/4))


def prefix_ratio(s,t,m):
    return sum((sum(s**(k-1-j)*t**j for j in range(k))**2/F(4)**k
                for k in range(1,m+1)),F(0))


def pi_candidates():
    lo,hi=pi_interval()
    # Coarsen outward to a compact exact decimal interval before high powers.
    scale=10**40
    lo=F((lo*scale).numerator//(lo*scale).denominator,scale)
    u=hi*scale
    hi=F(-((-u.numerator)//u.denominator),scale)
    lo-=3;hi-=3
    rows=[]
    for p,q in ((22,7),(355,113)):
        a=F(p,q)-3
        assert 0<lo<hi<a<1
        kl,kh=q*q*(a-hi),q*q*(a-lo)
        rl,rh=distance_ratio(a,lo),distance_ratio(a,hi)
        assert F(1,4)<=rl<rh<=F(20,27)
        prefixes=[]
        for m in (1,2,4,8,12):
            pl,ph=prefix_ratio(a,lo,m),prefix_ratio(a,hi,m)
            assert F(1,4)<=pl<=ph
            assert pl<rl and ph<rh
            prefixes.append(dict(m=m,ratio_squared_lower=pl,ratio_squared_upper=ph,
                                  common_coordinate_denominator=(a/2).denominator**m))
        rows.append(dict(p=p,q=q,scalar_kappa_lower=kl,scalar_kappa_upper=kh,
                         hilbert_kappa_squared_lower=kl*kl*rl,
                         hilbert_kappa_squared_upper=kh*kh*rh,
                         ratio_squared_lower=rl,ratio_squared_upper=rh,
                         scalar_kappa_display=float((kl+kh)/2),
                         hilbert_kappa_display=sqrt(float((kl*kl*rl+kh*kh*rh)/2)),
                         ratio_display=sqrt(float((rl+rh)/2)),prefixes=prefixes,
                         q_is_original_parameter_denominator=True))
    return rows


def exact_small_volume_cost():
    rows=[]
    for Q in (2,3,5,10,100,1000):
        eps=F(1,Q);nodes=[F(0),eps,2*eps]
        G=[[1/(1-s*t/4) for t in nodes] for s in nodes]
        exact=determinant(G)
        den=(4*Q*Q-1)*(Q*Q-1)*(2*Q*Q-1)**2
        assert exact==F(Q*Q,den) and gcd(Q*Q,den)==1
        assert exact*den==Q*Q
        first_three=eps**6/16
        assert first_three<exact
        rows.append(dict(Q=Q,infinite_gram_determinant=exact,
                         first_three_coordinate_determinant=first_three,
                         true_reduced_denominator=den,
                         cleared_integer=Q*Q,volume_display=sqrt(float(exact))))
    return rows


def free_target_control():
    rows=[]
    for m in range(2,6):
        from math import factorial
        q=10**factorial(m)
        prefix=sum((F(1,10**factorial(j)) for j in range(1,m+1)),F(0))
        assert prefix.denominator==q
        upper=F(2,10**factorial(m+1))
        assert upper==2*F(1,q)**(m+1)
        rows.append(dict(m=m,denominator_decimal_digits=factorial(m)+1,
                         guaranteed_tail_bound='2*q^(-'+str(m+1)+')'))
    return dict(target='L=sum_{j>=1}10^(-j!)',prefixes=rows,
                scope='Constructing a chosen target; not a construction of pi approximants.')


def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x


def main():
    data=dict(status='passed',pi_candidates=pi_candidates(),
              exact_small_volume_cost=exact_small_volume_cost(),
              freely_chosen_target=free_target_control(),
              pi_bound_method='Machin rational interval, coarsened outward to 40 decimal places',
              scope=__doc__)
    (HERE/'pi_transfer_results.json').write_text(json.dumps(serial(data),indent=2)+'\n')
    print(json.dumps(dict(status='passed',pi_candidates=2,volume_cost_examples=6,
                          free_target_controls=4)))


if __name__=='__main__':main()
