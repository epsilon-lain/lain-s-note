"""Optimal prefix-cell bounds and an oscillatory Liouville control.

Only exact integer/Fraction decisions; display decimals are not certificates.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
import random

from verify import atan_bounds, decimal_string

ROOT=Path(__file__).parent


def norm(x):
    t=x-x.__floor__()
    return min(t,1-t)


def norm_range(lo,hi):
    """Exact extrema of distance to Z on a closed rational interval."""
    assert lo<=hi
    integer_inside=-((-lo).__floor__())<=hi.__floor__()
    half_inside=-((-(lo-F(1,2))).__floor__())<=(hi-F(1,2)).__floor__()
    low=F() if integer_inside else min(norm(lo),norm(hi))
    high=F(1,2) if half_inside else max(norm(lo),norm(hi))
    return low,high


def cell_certificate(lower,upper,p,length=3):
    """Best possible bounds from the prefix alone, using cell closure."""
    factor=p**length
    a=(factor*lower).__floor__()
    assert a==(factor*upper).__floor__(), 'Insufficient source precision'
    q=p-1
    low,high=norm_range(F(q*a,factor),F(q*(a+1),factor))
    actual_lo,actual_hi=norm_range(q*lower,q*upper)
    assert q*low<=q*actual_lo<=q*actual_hi<=q*high
    digits=[];v=a
    for _ in range(length):
        digits.append(v%p);v//=p
    return {'base':str(p),'q':str(q),'prefix_digits':list(reversed(digits)),
            'prefix_length':length,'lower_exact':str(q*low),
            'upper_exact':str(q*high),
            'lower_display':decimal_string(q*low,22),
            'upper_display':decimal_string(q*high,22),
            'source_interval_lower_display':decimal_string(q*actual_lo,22),
            'source_interval_upper_display':decimal_string(q*actual_hi,22)}


def pi_records():
    al,au=atan_bounds(5,800);bl,bu=atan_bounds(239,250)
    raw_lo,raw_hi=16*al-4*bu,16*au-4*bl
    scale=10**1000
    lo=F((raw_lo*scale).__floor__(),scale)
    hi=F(-((-raw_hi*scale).__floor__()),scale)
    assert hi-lo==F(1,scale)
    cflo,cfhi=lo,hi
    a2,a1,b2,b1=0,1,1,0
    best_hi=None;records=[]
    for index in range(501):
        coef=cflo.__floor__()
        assert coef==cfhi.__floor__()
        numerator,q=coef*a1+a2,coef*b1+b2
        if q>=2:
            nearlo=(q*lo+F(1,2)).__floor__()
            nearhi=(q*hi+F(1,2)).__floor__()
            assert nearlo==nearhi==numerator
            dlow,dhigh=norm_range(q*lo,q*hi)
            el,eh=q*dlow,q*dhigh
            if best_hi is None or eh<best_hi:
                cert=cell_certificate(lo-3,hi-3,q+1)
                cert.update({'cf_index':index,'pi_numerator':str(numerator),
                             'denominator_decimal_digits':len(str(q)),
                             'E_interval_exact':[str(el),str(eh)]})
                records.append(cert);best_hi=eh
        cflo,cfhi=1/(cfhi-coef),1/(cflo-coef)
        a2,a1,b2,b1=a1,numerator,b1,q
    # A finite witness, independently certified from only three digits.
    assert F(records[-1]['upper_exact'])<F(1,20000)
    return records


def oscillatory_control():
    # y=0.overline(46)+sum_{j>=1}10^{-j!}; no carries in the digits.
    partial=F(46,99)+sum((F(1,10**factorial(j)) for j in range(1,6)),F())
    lower=partial+F(1,10**720)
    upper=partial+F(2,10**720)
    rows=[]
    for k in range(2,6):
        r=F(46,99)+sum((F(1,10**factorial(j)) for j in range(1,k+1)),F())
        q=99*10**factorial(k)
        assert (q*r).denominator==1
        bound=2*99**2*F(1,10**((k-1)*factorial(k)))
        assert 0<q*q*(lower-r)<q*q*(upper-r)<=bound
        row=cell_certificate(lower,upper,q+1)
        row.update({'k':k,'chosen_rational':str(r),
                    'chosen_representation_q':str(q),
                    'analytic_q_squared_error_upper':str(bound),
                    'analytic_upper_display':decimal_string(bound,22)})
        rows.append(row)
    n=10000
    special={factorial(j) for j in range(1,8)}
    digits=[(4 if i%2 else 6)+int(i in special) for i in range(1,n+1)]
    gaps=[abs(b-a) for a,b in zip(digits,digits[1:])]
    assert min(gaps)>=2
    energy=sum(g*g for g in gaps)
    return {'prefix_decimal_digits':''.join(map(str,digits[:80])),
            'checked_digits':n,'minimum_adjacent_digit_gap':min(gaps),
            'mean_squared_gap_exact':str(F(energy,n-1)),
            'mean_squared_gap_display':decimal_string(F(energy,n-1),22),
            'analytic_limit_mean_squared_gap':'4',
            'approximation_certificates':rows}


def range_checks():
    rng=random.Random(20261028)
    # Endpoints, integers, and half integers exhaust piecewise-linear extrema.
    for _ in range(3000):
        lo=F(rng.randrange(-100,100),rng.randrange(1,30))
        hi=lo+F(rng.randrange(0,80),rng.randrange(1,30))
        low,high=norm_range(lo,hi)
        points=[lo,hi]
        for k in range(2*lo.__floor__()-2,2*hi.__floor__()+4):
            if lo<=F(k,2)<=hi:points.append(F(k,2))
        assert low==min(map(norm,points)) and high==max(map(norm,points))
    return 3000


if __name__=='__main__':
    result={'interval_extrema_checks':range_checks(),
            'pi_records_among_convergents_0_to_500':pi_records(),
            'oscillatory_liouville_control':oscillatory_control(),
            'claim':'finite pi witnesses and proved control construction; no pi asymptotic verdict'}
    (ROOT/'wave_certificate_results.json').write_text(json.dumps(result,indent=2)+'\n')
    for row in result['pi_records_among_convergents_0_to_500']:
        print('pi convergent',row['cf_index'],'q digits',row['denominator_decimal_digits'],
              'E in',row['lower_display'],row['upper_display'])
    control=result['oscillatory_liouville_control']
    print('Control prefix:',control['prefix_decimal_digits'])
    print('Control minimum gap:',control['minimum_adjacent_digit_gap'],
          'energy:',control['mean_squared_gap_display'])
    for row in control['approximation_certificates']:
        print('Control k',row['k'],'analytic E upper',row['analytic_upper_display'])
