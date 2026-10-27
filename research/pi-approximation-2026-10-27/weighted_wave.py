"""Certified weighted-slope and period-return certificates, stdlib only."""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json
import random

from verify import certified_digits, pi_bounds, rounded_enclosure

ROOT = Path(__file__).parent


def norm(a):
    t = a - a.__floor__()
    return min(t, 1-t)


def analyze(label, lower, upper, scale, q, length=12):
    p = q+1
    digits = certified_digits(lower, upper, scale, p, length+1)
    # Ensure the nearest numerator is constant on the whole interval.
    n = (F(q*lower, scale)+F(1,2)).__floor__()
    assert n == (F(q*upper, scale)+F(1,2)).__floor__()
    signed_lo, signed_hi = F(q*lower,scale)-n, F(q*upper,scale)-n
    assert signed_lo*signed_hi > 0
    elo = q*min(abs(signed_lo), abs(signed_hi))
    ehi = q*max(abs(signed_lo), abs(signed_hi))
    out=[]
    for j in [2,4,length]:
        slopes = sum((F(digits[k]-digits[k-1],p**k)
                      for k in range(1,j+1)), F())
        tail=F(1,p**j)
        lo=q*max(F(), norm(slopes)-tail)
        hi=q*(norm(slopes)+tail)
        assert lo <= elo <= ehi <= hi
        out.append({'weighted_differences':j,
                    'truncated_sum_exact':str(slopes),
                    'normalized_error_lower_exact':str(lo),
                    'normalized_error_upper_exact':str(hi),
                    'lower_display':float(lo), 'upper_display':float(hi)})
    return {'target':label,'base':p,'q':q,'nearest_n':n,'digits':digits,
            'normalized_error_interval_exact':[str(elo),str(ehi)],
            'normalized_error_display':float((elo+ehi)/2),
            'prefix_certificates':out}


def identity_checks():
    rng=random.Random(27102026)
    count=0
    for _ in range(3000):
        p=rng.randrange(2,100)
        x=F(rng.randrange(1,100000),100001)
        for k in [1,2,3]:
            a=(p**k*x).__floor__()
            t=p**k*x-a
            assert (p**k-1)*x == a+t-x
            n=(((p**k-1)*x)+F(1,2)).__floor__()
            assert norm(t-x) == abs((p**k-1)*x-n)
            # Check the finite slope sum against its exact tail identity.
            total=F(); prev=(p*x).__floor__()
            for j in range(1,7):
                curr=(p**(j+1)*x).__floor__()-p*(p**j*x).__floor__()
                total+=F(curr-prev,p**j)
                prev=curr
            rem=(p**6*x)-(p**6*x).__floor__()
            next_rem=p*rem-(p*rem).__floor__()
            assert (p-1)*x-(p*x).__floor__() == total+F(next_rem-rem,p**6)
            count+=1
    return count


def main():
    lo,hi=pi_bounds()
    lower,upper,scale=rounded_enclosure(lo,hi)
    lower-=3*scale;upper-=3*scale
    r=isqrt(2*scale*scale)
    rows=[]
    for q in [7,113,1725033,6701487259]:
        rows.append(analyze('pi-3',lower,upper,scale,q))
    for q in [5,12,29,70,169,408,985]:
        rows.append(analyze('sqrt2-1',r-scale,r+1-scale,scale,q))
    out={'exact_identity_checks':identity_checks(),'rows':rows,
         'claim':'finite certified weighted-slope intervals; no pi asymptotic conclusion'}
    (ROOT/'weighted_wave_results.json').write_text(json.dumps(out,indent=2)+'\n')
    for row in rows:
        print(row['target'],'base',row['base'],'E',row['normalized_error_display'],
              'four-slope bounds',row['prefix_certificates'][1]['lower_display'],
              row['prefix_certificates'][1]['upper_display'])
    print('Exact identity checks:',out['exact_identity_checks'])


if __name__=='__main__':
    main()
