"""Greedy real-base wave cells and certified flat-prefix boundaries.

All decisions use integer/Fraction arithmetic; decimal strings are display only.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import random

from verify import pi_bounds, rounded_enclosure, decimal_string

ROOT=Path(__file__).parent


def greedy_digits(base, value, count):
    assert base>1 and 0<=value<1
    remainder=value
    out=[]
    for _ in range(count):
        scaled=base*remainder
        digit=scaled.__floor__()
        remainder=scaled-digit
        assert 0<=remainder<1
        out.append(digit)
    return out


def flat_prefix_bounds(base, digit, length):
    assert digit>=1 and base>digit+1 and length>=1
    tail=base**(-length)
    lower=digit*(1-tail)/(base-1)
    return lower,lower+tail


def root_bracket(value_lo,value_hi,digit,length,kind,bits=240):
    """Enclose the root for either prefix-cell boundary by rational bisection."""
    left=F(digit+1)
    right=F(2)+2*F(digit)/value_lo
    if kind=='lower':
        def boundary(b):
            return digit*(1-b**(-length))/(b-1)
    elif kind=='upper':
        def boundary(b):
            return digit*(1-b**(-length))/(b-1)+b**(-length)
    else:
        raise ValueError(kind)
    assert boundary(left)>value_hi and boundary(right)<value_lo
    for _ in range(bits):
        mid=(left+right)/2
        if boundary(mid)>value_hi:
            left=mid
        elif boundary(mid)<value_lo:
            right=mid
        else:
            raise ArithmeticError('Increase source enclosure precision')
    assert boundary(left)>value_hi and boundary(right)<value_lo
    return left,right


def roots(lo,hi,digit,length):
    ll,lu=root_bracket(lo,hi,digit,length,'lower')
    ul,uu=root_bracket(lo,hi,digit,length,'upper')
    center_lo,center_hi=1+F(digit)/hi,1+F(digit)/lo
    dl=(center_lo-lu,center_hi-ll)
    dr=(ul-center_hi,uu-center_lo)
    assert ll<=lu<center_lo<=center_hi<ul<=uu and dl[0]>0 and dr[0]>0
    width_lo=ul-lu; width_hi=uu-ll
    center=(center_lo+center_hi)/2
    ratio_lo=dr[0]/dl[1];ratio_hi=dr[1]/dl[0]
    scaled_lo=center_lo**length*dl[0]
    scaled_hi=center_hi**length*dl[1]
    scaled_right_lo=center_lo**length*dr[0]
    scaled_right_hi=center_hi**length*dr[1]
    return {'digit':digit,'prefix_length':length,
            'left_root_interval_exact':[str(ll),str(lu)],
            'right_root_interval_exact':[str(ul),str(uu)],
            'left_base_display':decimal_string((ll+lu)/2,25),
            'center_base_display':decimal_string(center,25),
            'right_base_display':decimal_string((ul+uu)/2,25),
            'left_radius_display':decimal_string(sum(dl)/2,20),
            'right_radius_display':decimal_string(sum(dr)/2,20),
            'width_display':decimal_string((width_lo+width_hi)/2,20),
            'right_to_left_radius_lower_exact':str(ratio_lo),
            'right_to_left_radius_upper_exact':str(ratio_hi),
            'right_to_left_radius_display':decimal_string((ratio_lo+ratio_hi)/2,18),
            'scaled_left_radius_display':decimal_string((scaled_lo+scaled_hi)/2,18),
            'scaled_right_radius_display':decimal_string((scaled_right_lo+scaled_right_hi)/2,18)}


def certified_integer_base(lo,hi,base,length):
    low,high=lo,hi
    digits=[]
    for _ in range(length):
        low*=base;high*=base
        digit=low.__floor__()
        assert digit==high.__floor__()
        digits.append(digit);low-=digit;high-=digit
    return digits


def run_checks():
    rng=random.Random(27102026)
    count=0
    for _ in range(7000):
        n=rng.randrange(1,6)
        base=F(n+1)+F(rng.randrange(1,40),rng.randrange(1,15))
        y=F(rng.randrange(0,10000),10000)
        length=rng.randrange(1,9)
        low,high=flat_prefix_bounds(base,n,length)
        assert (greedy_digits(base,y,length)==[n]*length)==(low<=y<high)
        # Every candidate starting inside a prefix cell satisfies recurrence.
        for start in (low, (low+high)/2):
            if start<1:
                assert greedy_digits(base,start,length)==[n]*length
        if high<1:
            assert greedy_digits(base,high,length)!=[n]*length
        count+=1
    return count


def main():
    low,high=pi_bounds()
    a,b,scale=rounded_enclosure(low,high)
    saved=json.loads((ROOT/'pi_interval.json').read_text())
    assert (a,b,scale)==tuple(map(int,(saved['lower'],saved['upper'],saved['scale'])))
    xlo,xhi=F(a-3*scale,scale),F(b-3*scale,scale)
    result={'pi_enclosure_regenerated':True,'exact_cell_tests':run_checks(),
            'target':'pi-3',
            'bases_8_and_114_digits':{
                '8':certified_integer_base(xlo,xhi,8,8),
                '114':certified_integer_base(xlo,xhi,114,8)},
            'rows':[],'n_growth_rows':[],
            'status':'finite certificates and elementary proofs only'}
    assert result['bases_8_and_114_digits']['8'][:3]==[1,1,0]
    assert result['bases_8_and_114_digits']['114'][:3]==[16,16,15]
    for n,lengths in [(1,[1,2,3,4,5,10,20,50]),(16,[1,2,3,4,8,12,20])]:
        for length in lengths:
            row=roots(xlo,xhi,n,length)
            result['rows'].append(row)
            print('n',n,'K',length,'left',row['left_base_display'],
                  'right',row['right_base_display'],'width',row['width_display'],
                  'right/left',row['right_to_left_radius_display'])
    for length in (2,3):
        for n in (1,2,4,8,16,32,64,128,1024):
            row=roots(xlo,xhi,n,length)
            left_lo,left_hi=map(F,row['left_root_interval_exact'])
            right_lo,right_hi=map(F,row['right_root_interval_exact'])
            ratio_lo=(right_lo-left_hi)*n**(length-1)/xhi**(length-2)
            ratio_hi=(right_hi-left_lo)*n**(length-1)/xlo**(length-2)
            result['n_growth_rows'].append({'n':n,'prefix_length':length,
                'ratio_lower_exact':str(ratio_lo),'ratio_upper_exact':str(ratio_hi),
                'ratio_display':decimal_string((ratio_lo+ratio_hi)/2,18)})
        print('n-growth K',length,'n=1024 ratio',result['n_growth_rows'][-1]['ratio_display'])
    by_key={(row['digit'],row['prefix_length']):row for row in result['rows']}
    for n,integer_base in ((1,8),(16,114)):
        two=by_key[n,2]
        three=by_key[n,3]
        assert F(two['left_root_interval_exact'][1]) < integer_base
        assert integer_base < F(two['right_root_interval_exact'][0])
        assert integer_base < F(three['left_root_interval_exact'][0])
    result['integer_base_cell_membership_certified']=True
    (ROOT/'base_space_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Exact cell checks:',result['exact_cell_tests'])
    print('Digits:',result['bases_8_and_114_digits'])


if __name__=='__main__':
    main()
