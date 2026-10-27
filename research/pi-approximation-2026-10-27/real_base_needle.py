"""Rational certificates for the real-base / rotating-needle model.

Uses exact rational arithmetic and an explicit alternating arctangent bound.
The reported decimal strings are only for display.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import random

from verify import pi_bounds, rounded_enclosure, decimal_string

ROOT = Path(__file__).parent


def atan_lower(z):
    """z-z^3/3 <= atan(z) <= z for 0<=z<=1."""
    assert 0 <= z <= 1
    return z-z**3/3


def case(pi_lo, pi_hi, q, n):
    x_lo, x_hi = pi_lo-3, pi_hi-3
    assert 0 < x_lo < x_hi < 1
    assert q >= 2 and 1 <= n <= q-1
    p=q+1
    beta_lo=1+F(n,1)/x_hi
    beta_hi=1+F(n,1)/x_lo
    assert not beta_lo<=p<=beta_hi
    bdist_lo=min(abs(beta_lo-p),abs(beta_hi-p))
    bdist_hi=max(abs(beta_lo-p),abs(beta_hi-p))
    r=3+F(n,q)
    assert not pi_lo<=r<=pi_hi
    dif_lo=min(abs(pi_lo-r),abs(pi_hi-r))
    dif_hi=max(abs(pi_lo-r),abs(pi_hi-r))
    norm_lo=q*q*dif_lo
    norm_hi=q*q*dif_hi
    # Monotonicity on either side of r ensures the extrema occur at endpoints.
    z1=abs(pi_lo-r)/(1+pi_lo*r)
    z2=abs(pi_hi-r)/(1+pi_hi*r)
    z_lo,z_hi=min(z1,z2),max(z1,z2)
    area_factor=F(q*q+(3*q+n)**2,2)
    area_lo=area_factor*atan_lower(z_lo)
    area_hi=area_factor*z_hi
    assert F(5,17)*norm_lo <= area_lo <= area_hi <= F(17,20)*norm_hi
    # Interval version of q^2*|x-n/q| = x*q*|p-beta_n|.
    assert x_lo*q*bdist_lo <= norm_hi
    assert x_hi*q*bdist_hi >= norm_lo
    return {'q':q,'a':3*q+n,'n':n,'integer_base_p':p,
            'beta_lower_exact':str(beta_lo),'beta_upper_exact':str(beta_hi),
            'beta_display':decimal_string((beta_lo+beta_hi)/2,25),
            'beta_gap_display':decimal_string((bdist_lo+bdist_hi)/2,22),
            'q2_error_lower_exact':str(norm_lo),'q2_error_upper_exact':str(norm_hi),
            'q2_error_display':decimal_string((norm_lo+norm_hi)/2,22),
            'needle_area_lower_exact':str(area_lo),'needle_area_upper_exact':str(area_hi),
            'needle_area_lower_display':decimal_string(area_lo,22),
            'needle_area_upper_display':decimal_string(area_hi,22)}


def identity_checks():
    rng=random.Random(20261027)
    for _ in range(12000):
        x=F(rng.randrange(1,10000),10001)
        q=rng.randrange(2,10000)
        n=rng.randrange(1,q)
        beta=1+F(n,1)/x
        p=q+1
        assert beta*x == n+x
        assert beta > n
        assert x-F(n,q) == x*(p-beta)/q
        assert q*q*abs(x-F(n,q)) == x*q*abs(p-beta)
    return 12000


def main():
    lo,hi=pi_bounds()
    lower,upper,scale=rounded_enclosure(lo,hi)
    supplied=json.loads((ROOT/'pi_interval.json').read_text())
    assert (lower,upper,scale)==tuple(map(int,
        (supplied['lower'],supplied['upper'],supplied['scale'])))
    pilo,pihi=F(lower,scale),F(upper,scale)
    rows=[case(pilo,pihi,q,n) for q,n in
          [(7,1),(113,16),(33102,4687),(1725033,244252)]]
    assert F(rows[1]['needle_area_upper_exact']) < F(1,85)
    out={'pi_enclosure_regenerated':True,'exact_identity_checks':identity_checks(),
         'arctan_certification':'z-z^3/3 <= atan(z) <= z; sector bounds exact',
         'pi_355_over_113_area_strictly_below_1_over_85':True,
         'results':rows,
         'claim':'finite examples and elementary equivalence; no pi asymptotic verdict'}
    (ROOT/'real_base_needle_results.json').write_text(json.dumps(out,indent=2)+'\n')
    for row in rows:
        print('q',row['q'],'a',row['a'],'beta',row['beta_display'],
              'base gap',row['beta_gap_display'],'E',row['q2_error_display'],
              'area in',row['needle_area_lower_display'],row['needle_area_upper_display'])
    print('Exact identity checks:',out['exact_identity_checks'])


if __name__=='__main__':
    main()
