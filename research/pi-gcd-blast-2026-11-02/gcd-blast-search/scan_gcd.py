"""Exact finite attack on the Legendre/arctangent gcd threshold.

Python standard library only. Quality bounds use the proved positive-error
integral, so the main scan requires no decimal expansion of pi.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path
import time


def encode_integer(value):
    return value.to_bytes(max(1, (value.bit_length()+7)//8), 'big')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--limit',type=int,default=5000)
    args = parser.parse_args()
    if not __debug__ or args.limit < 150:
        parser.error('Assertions must be enabled and limit must be at least 150.')
    started=time.monotonic()
    ap,a,bp,b=1,1,0,1
    fact=1
    hbest=Fraction(0)
    best_n=None
    after_two=Fraction(0)
    after_two_n=None
    min_quality_floor=None
    min_floor_n=None
    odd_cases=[]
    largest_odd_numerator=(0,1)
    selections=[]
    chain=hashlib.sha256()
    for n in range(1,args.limit+1):
        if n>1:
            k=n-1
            ap,a=a,(2*k+1)*a+k*k*ap
            bp,b=b,(2*k+1)*b+k*k*bp
            fact*=n
        g=gcd(a,4*b)
        common=gcd(g,fact)
        num,den=g//common,fact//common
        assert gcd(num,den)==1
        hsq=Fraction(n*num*num,den*den)
        qfloor=Fraction(2*den*den,(2*n+1)*num*num)
        qceil=2*qfloor
        if hsq>hbest:
            hbest,best_n=hsq,n
        if n>=3 and hsq>after_two:
            after_two,after_two_n=hsq,n
        if n!=2 and (min_quality_floor is None or qfloor<min_quality_floor):
            min_quality_floor,min_floor_n=qfloor,n
        odd=num >> ((num & -num).bit_length()-1)
        if odd>1:
            if len(odd_cases)<30:
                odd_cases.append({'n':n,'g_over_factorial_reduced_numerator':str(num),
                                  'odd_part':str(odd)})
            if odd>largest_odd_numerator[1]: largest_odd_numerator=(n,odd)
        for integer in (n,a,b,g,fact):
            payload=encode_integer(integer)
            chain.update(len(payload).to_bytes(8,'big'))
            chain.update(payload)
        if n in (1,2,3,4,6,10,72,100,107,1000,2000,5000,args.limit):
            selections.append({'n':n,'g_over_factorial_numerator':str(num),
                               'g_over_factorial_denominator_bits':den.bit_length(),
                               'A_bits':a.bit_length(),'g_bits':g.bit_length(),
                               'reduced_q_bits':(a//g).bit_length(),
                               'h_squared_numerator':str(hsq.numerator),
                               'h_squared_denominator_bits':hsq.denominator.bit_length(),
                               'quality_lower_bound_numerator_bits':qfloor.numerator.bit_length(),
                               'quality_lower_bound_denominator':str(qfloor.denominator)})
    assert best_n==2 and hbest==8
    assert after_two_n==3 and after_two==Fraction(4,3)
    assert min_floor_n==3 and min_quality_floor==Fraction(9,14)
    # At n=2 the analytic bound is Q_2<1/5, strictly below all other floors.
    assert Fraction(1,5)<min_quality_floor
    result={'status':'passed','limit':args.limit,
            'arithmetic':'exact integers; no numerical pi',
            'quality_bound_formula':'2/(2n+1)*(n!/g_n)^2 < Q_n < 4/(2n+1)*(n!/g_n)^2',
            'max_h_squared':{'n':best_n,'value':str(hbest)},
            'max_h_squared_for_n_at_least_3':{'n':after_two_n,'value':str(after_two)},
            'minimum_quality_lower_bound_excluding_n2':{'n':min_floor_n,'value':str(min_quality_floor)},
            'unique_best_quality_index_in_checked_family':2,
            'best_fraction':'3/1',
            'first_odd_excess_cases':odd_cases,
            'largest_odd_reduced_numerator':{'n':largest_odd_numerator[0],
                                             'value':str(largest_odd_numerator[1])},
            'recurrence_gcd_factorial_chain_sha256':chain.hexdigest(),
            'selected_rows':selections,
            'runtime_seconds':round(time.monotonic()-started,3),
            'scope':'Complete scan of this specified recurrence for 1<=n<=limit. It is not a scan of all rational numbers, an all-n bound, or a pi bad-approximability result.'}
    Path(__file__).with_name('gcd_scan_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('selected_rows','first_odd_excess_cases')},indent=2))


if __name__=='__main__': main()
