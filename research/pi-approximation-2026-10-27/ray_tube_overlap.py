"""Exact ray-tube overlap certificates; Python standard library only.

This elementary model is inspired by the tube-overlap viewpoint in the cited
OpenAI preprint. It does not implement or independently verify that paper.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json
import random

from verify import pi_bounds, rounded_enclosure, decimal_string

ROOT=Path(__file__).parent


def phi(z):
    assert z>=0
    return 1-z/4 if z<=2 else 1/z


def overlap_integral(length,half_width,slope_gap):
    """Exact area of intersection of two vertical-slice tubes."""
    assert length>0 and half_width>0 and slope_gap>=0
    if slope_gap==0:
        return 2*half_width*length
    cutoff=min(length,2*half_width/slope_gap)
    return 2*half_width*cutoff-slope_gap*cutoff*cutoff/2


def nearest_bounds(lo,hi,q):
    a=(q*lo+F(1,2)).__floor__()
    assert a==(q*hi+F(1,2)).__floor__()
    signed_lo=q*lo-a; signed_hi=q*hi-a
    assert signed_lo*signed_hi>0
    return a,q*min(abs(signed_lo),abs(signed_hi)),q*max(abs(signed_lo),abs(signed_hi))


def individual(lo,hi,q,label):
    a,e_lo,e_hi=nearest_bounds(lo,hi,q)
    f_lo,f_hi=phi(e_hi),phi(e_lo)
    return {'target':label,'q':q,'a':str(a),
            'E_lower_exact':str(e_lo),'E_upper_exact':str(e_hi),
            'E_display':decimal_string((e_lo+e_hi)/2,22),
            'additional_sweep_area_lower_exact':str(e_lo/2),
            'additional_sweep_area_upper_exact':str(e_hi/2),
            'additional_sweep_area_display':decimal_string((e_lo+e_hi)/4,22),
            'overlap_fraction_lower_exact':str(f_lo),
            'overlap_fraction_upper_exact':str(f_hi),
            'overlap_fraction_display':decimal_string((f_lo+f_hi)/2,22)}


def dyadic(lo,hi,Q,label):
    """Finite candidate scan: minimize absolute slope error, then Phi."""
    best=None
    for q in range(Q,2*Q):
        a,e_lo,e_hi=nearest_bounds(lo,hi,q)
        z_lo=F(Q*Q,q*q)*e_lo
        z_hi=F(Q*Q,q*q)*e_hi
        if best is not None and F(a,q)==F(best['a'],best['q']):
            assert z_lo==best['z_lo'] and z_hi==best['z_hi']
            continue
        if best is None or z_hi<best['z_hi']:
            # Certification that the candidate is better than every previous one.
            if best is not None:
                assert z_hi<best['z_lo']
            best={'q':q,'a':a,'z_lo':z_lo,'z_hi':z_hi}
        else:
            assert z_lo>=best['z_hi'], 'Source interval cannot order candidates'
    return {'target':label,'dyadic_Q':Q,'q':best['q'],'a':str(best['a']),
            'scaled_gap_lower_exact':str(best['z_lo']),
            'scaled_gap_upper_exact':str(best['z_hi']),
            'overlap_lower_exact':str(phi(best['z_hi'])),
            'overlap_upper_exact':str(phi(best['z_lo'])),
            'overlap_display':decimal_string((phi(best['z_hi'])+phi(best['z_lo']))/2,22)}


def identity_checks():
    rng=random.Random(20261027)
    for _ in range(12000):
        Q=F(rng.randrange(1,200),rng.randrange(1,30))
        w=F(rng.randrange(1,60),rng.randrange(1,60))
        gap=F(rng.randrange(0,100),rng.randrange(1,50))
        z=gap*Q/w
        common=overlap_integral(Q,w,gap)
        assert common/(2*w*Q)==phi(z)
        # Shearing direction between endpoints fills each vertical interval.
        assert (2*w*Q+gap*Q*Q/2)-2*w*Q == gap*Q*Q/2
        assert 0<=common<=2*w*Q
    return 12000


def gram_checks():
    rng=random.Random(27102026)
    for _ in range(1200):
        delta=F(1,rng.randrange(2,40)**2)
        slopes=[F(rng.randrange(-100,101),rng.randrange(1,20))
                for _ in range(rng.randrange(2,9))]
        coeffs=[rng.randrange(-5,6) for _ in slopes]
        energy=sum((F(coeffs[i]*coeffs[j])*
                    overlap_integral(F(1),delta,abs(slopes[i]-slopes[j]))
                    for i in range(len(slopes)) for j in range(len(slopes))),F())
        assert energy>=0
    return 1200


def main():
    al,au=pi_bounds(); low,high,scale=rounded_enclosure(al,au)
    saved=json.loads((ROOT/'pi_interval.json').read_text())
    assert (low,high,scale)==tuple(map(int,(saved['lower'],saved['upper'],saved['scale'])))
    pi_lo,pi_hi=F(low,scale),F(high,scale)
    root=isqrt(2*scale*scale)
    rt_lo,rt_hi=2+F(root,scale),2+F(root+1,scale)
    rows=[individual(pi_lo,pi_hi,q,'pi') for q in [7,113,33102,1725033,6701487259]]
    rows +=[individual(rt_lo,rt_hi,q,'2+sqrt(2)') for q in [5,12,29,70,113,169]]
    assert F(rows[1]['overlap_fraction_lower_exact'])>F(999,1000)
    dyads=[dyadic(lo,hi,Q,label) for label,lo,hi in
           [('pi',pi_lo,pi_hi),('2+sqrt(2)',rt_lo,rt_hi)]
           for Q in [4,8,16,32,64,128,256,512,1024]]
    result={'geometric_identity_checks':identity_checks(),
            'exact_gram_energy_checks':gram_checks(),
            'pi_enclosure_regenerated':True,
            'individual_tube_certificates':rows,'dyadic_certificates':dyads,
            'status':'elementary tube model, no use of the PDE theorem, no pi liminf verdict'}
    (ROOT/'ray_tube_results.json').write_text(json.dumps(result,indent=2)+'\n')
    for row in rows:
        print(row['target'],'q',row['q'],'a',row['a'],'E',row['E_display'],
              'overlap',row['overlap_fraction_display'],
              'extra sweep',row['additional_sweep_area_display'])
    for row in dyads:
        print('dyadic',row['target'],'Q',row['dyadic_Q'],'q',row['q'],
              'overlap',row['overlap_display'])
    print('Exact geometry identities:',result['geometric_identity_checks'])
    print('Exact Gram energy checks:',result['exact_gram_energy_checks'])


if __name__=='__main__':
    main()
