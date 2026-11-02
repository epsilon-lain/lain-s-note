"""Certified necessary-condition sieve for Q_n <= 1/5 in one fixed family.

The odd-prime Wronskian theorem and U_n=2^n ell_n give a lower bound
using only U_n and an lcm. Exact A,B,g are computed only up to the last
surviving index. No decimal pi and no floating-point comparison are used.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path
import time


def v2(x):
    return (x & -x).bit_length()-1


def odd(x):
    return x >> v2(x)


def update_hash(chain, values):
    for value in values:
        payload=value.to_bytes(max(1,(value.bit_length()+7)//8),'big')
        chain.update(len(payload).to_bytes(8,'big'))
        chain.update(payload)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--limit',type=int,default=20000)
    parser.add_argument('--max-exact-index',type=int,default=5000)
    args=parser.parse_args()
    if not __debug__ or args.limit<2:
        parser.error('Assertions must be enabled and limit must be at least 2.')
    started=time.monotonic()
    up,u=1,2
    lcm=1
    survivors=[]
    selections=[]
    chain=hashlib.sha256()
    for n in range(1,args.limit+1):
        if n>1:
            k=n-1
            numerator=2*(2*n-1)*u+4*(n-1)*up
            assert numerator%n==0
            up,u=u,numerator//n
            lcm=lcm*k//gcd(lcm,k)
        e=lcm//gcd(lcm,n)
        s=v2(u)
        big_g=gcd(odd(u),odd(e))
        # Strict floor Q > lower_num/lower_den. Treat negative exponents
        # exactly rather than silently truncating a rational power of two.
        exponent=2*(n-s)
        lower_num=2 << max(exponent,0)
        lower_den=(2*n+1)*big_g*big_g << max(-exponent,0)
        is_survivor=5*lower_num<lower_den
        if is_survivor:
            survivors.append(n)
        update_hash(chain,(n,u,e,big_g))
        if n in (2,3,4,6,7,72,100,1000,5000,10000,args.limit):
            selections.append({'n':n,'U_bits':u.bit_length(),'v2_U':s,
                               'odd_sieve_gcd_bits':big_g.bit_length(),
                               'lower_bound_log2_floor_at_least':
                                   lower_num.bit_length()-lower_den.bit_length()-1,
                               'survives_threshold_one_fifth':is_survivor})
    # A survivor only means the inexpensive bound did not exclude it.
    # All are checked with actual reduced denominators below.
    exact_limit=max(survivors,default=0)
    exact_complete=exact_limit<=args.max_exact_index
    rows=[]
    if exact_complete:
        ap,a,bp,b=1,1,0,1
        fact=1
        survivor_set=set(survivors)
        for n in range(1,exact_limit+1):
            if n>1:
                k=n-1
                ap,a=a,(2*k+1)*a+k*k*ap
                bp,b=b,(2*k+1)*b+k*k*bp
                fact*=n
            if n not in survivor_set: continue
            g=gcd(a,4*b)
            floor=Fraction(2*fact*fact,(2*n+1)*g*g)
            ceil=2*floor
            rows.append({'n':n,'p':str(4*b//g),'q':str(a//g),
                         'strict_quality_lower':str(floor),
                         'strict_quality_upper':str(ceil),
                         'excluded_by_exact_lower_bound':floor>=Fraction(1,5)})
        unresolved=[row['n'] for row in rows if not row['excluded_by_exact_lower_bound']]
    else:
        unresolved=survivors[:]
    unique_best=(exact_complete and unresolved==[2] and
                 next(row for row in rows if row['n']==2)['strict_quality_upper']=='1/5')
    result={'status':'passed','limit':args.limit,
            'threshold':'1/5','sieve_survivors':survivors,
            'sieve_excluded_count':args.limit-len(survivors),
            'last_exact_index_required':exact_limit,
            'exact_stage_completed':exact_complete,
            'remaining_indices_after_exact_stage':unresolved,
            'unique_best_quality_index_in_checked_family':2 if unique_best else None,
            'unique_best_fraction':'3/1' if unique_best else None,
            'sieve_rows':selections,'exact_survivor_rows':rows,
            'U_lcm_odd_gcd_chain_sha256':chain.hexdigest(),
            'runtime_seconds':round(time.monotonic()-started,3),
            'scope':'Finite complete exclusion in the specified Legendre/arctangent family only. Surviving the sieve is not a good-approximation certificate; exact stage and analytic bounds are required. No all-rational or all-n conclusion.'}
    Path(__file__).with_name('gcd_sieve_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('sieve_rows','exact_survivor_rows')},indent=2))


if __name__=='__main__': main()
