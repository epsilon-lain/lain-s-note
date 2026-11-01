#!/usr/bin/env python3
"""Separate exact small determinants from quantization-induced apparent zeros.

This is a finite test of a specified six-column minor and measurement rules,
not a pi approximation theorem or a physical field model.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys
import sympy as sp

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'pi-construction-mechanisms-2026-11-01'))
from period_lab import original_matrix
from source_support import audit_source_files
sys.path.insert(0,str(HERE.parent/'rational-scale-geometry-2026-11-01'))
from denominator_envelope_lab import pi_interval

COLUMNS=[(0,(0,0)),(1,(0,0)),(0,(0,1)),(0,(0,2)),(0,(0,3)),(0,(1,0))]
ROWS=[(j,s,(0,0)) for j in range(2) for s in range(3)]


def quantize(x,h):
    z=F(x)/F(h)+F(1,2)
    return F(h)*(z.numerator//z.denominator)


def matrix(r1,r2):return original_matrix(COLUMNS,ROWS,(r1,r2),(3,3))


def gaussian(z):return bool(sp.re(z).is_Integer and sp.im(z).is_Integer)


def intrinsic_and_observed():
    delta=sp.Symbol('delta',real=True)
    A=matrix(sp.I*(6+delta),6*sp.I)
    assert sp.simplify(A.det()-108*sp.I*delta)==0
    assert A[:,5]-A[:,2]==sp.Matrix([0,0,0,sp.I*delta,0,0])
    assert A.subs(delta,0).rank()==5
    rows=[]
    h=F(1,100)
    for Q in (10,100,1000,10**4,10**6):
        d=F(1,Q)
        actual=A.subs(delta,sp.Rational(1,Q))
        observed=matrix(sp.I*sp.Rational(quantize(6+d,h)),6*sp.I)
        assert actual.rank()==6
        # q1=2Q represents r1=2i*(6Q+1)/(2Q), and q2=1.
        # Two s=2 rows need a factor 2; all other factors are integral.
        cleared=sp.diag(1,1,2,1,1,2)*actual*sp.diag(1,1,1,1,1,2*Q)
        assert all(gaussian(x) for x in cleared)
        assert cleared.det()==864*sp.I
        normalized=actual.copy()
        normalized[:,5]=-sp.I*Q*(actual[:,5]-actual[:,2])
        assert normalized[:,5]==sp.Matrix([0,0,0,1,0,0])
        assert normalized.det()==108
        rows.append(dict(Q=Q,center_gap=d,exact_rank=6,
                         exact_determinant=str(sp.factor(actual.det())),
                         true_determinant_modulus=F(108,Q),
                         observation_h=h,quantized_center=quantize(6+d,h),
                         observed_rank=observed.rank(),
                         observed_determinant=str(sp.factor(observed.det())),
                         fixed_normalization_condition_number_lower_bound=2*Q,
                         cleared_integer_determinant='864*I',
                         difference_column_magnification=Q,
                         determinant_after_magnification='108'))
    assert rows[2]['observed_rank']==5
    return rows


def actual_rational_centers():
    c1,c2=F(44,7),F(710,113)
    A=matrix(sp.I*sp.Rational(c1),sp.I*sp.Rational(c2))
    assert A.rank()==6
    result=[]
    for h in (F(1,10),F(1,50),F(1,100),F(1,1000)):
        q1,q2=quantize(c1,h),quantize(c2,h)
        observed=matrix(sp.I*sp.Rational(q1),sp.I*sp.Rational(q2))
        result.append(dict(h=h,observed_centers=[q1,q2],observed_rank=observed.rank(),
                           observed_determinant=str(sp.factor(observed.det()))))
    assert result[1]['observed_rank']==5 and result[2]['observed_rank']==6
    return dict(actual_rank=6,exact_determinant=str(sp.factor(A.det())),observations=result)


def pi_observation():
    lo,hi=pi_interval()
    records=[]
    for p,q in ((22,7),(355,113)):
        a=F(p,q)
        assert a>hi
        true_lo,true_hi=q*q*(a-hi),q*q*(a-lo)
        assert 0<true_lo<true_hi<F(1,2)
        samples=[]
        for power in (1,2,3,4,6):
            h=F(1,q**power)
            observed=quantize(lo,h)
            assert observed==quantize(hi,h)
            measured=q*q*abs(observed-a)
            budget=q*q*h/2
            assert max(abs(measured-true_lo),abs(measured-true_hi))<=budget
            samples.append(dict(h=h,power=power,observed_pi=observed,
                                measured_kappa=measured,error_budget=budget,
                                measured_kappa_display=float(measured)))
        assert samples[1]['measured_kappa']==0
        assert samples[3]['measured_kappa']>0
        records.append(dict(p=p,q=q,true_kappa_lower=true_lo,true_kappa_upper=true_hi,
                             true_kappa_display=float((true_lo+true_hi)/2),observations=samples))
    return records


def candidate_resolution():
    """Exact target, changing legal candidate sets; no measurement quantization."""
    lo,hi=pi_interval()
    out=[]
    for n in (7,8,14,28):
        r=quantize(lo,F(1,n))
        assert r==quantize(hi,F(1,n))
        low,high=((lo-r,hi-r) if r<lo else (r-hi,r-lo))
        assert low>0
        out.append(dict(n=n,nearest_fraction=r,exact_distance_lower=low,
                        exact_distance_upper=high,display_distance=float((low+high)/2)))
    assert out[0]['exact_distance_upper']<out[1]['exact_distance_lower']
    assert out[0]['nearest_fraction']==out[2]['nearest_fraction']==out[3]['nearest_fraction']
    assert (out[0]['nearest_fraction']*8).denominator!=1
    return out


def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,sp.Integer):return int(x)
    if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x


def main():
    payload=dict(status='passed',scope=__doc__,
                 symbolic_determinant='108*I*delta',
                 intrinsic_and_quantized=intrinsic_and_observed(),
                 actual_pi_approximant_centers=actual_rational_centers(),
                 pi_observations=pi_observation(),
                 exact_target_changing_candidate_resolution=candidate_resolution(),
                 source_audit=audit_source_files(['build/main.tex','build/sections/determinant.tex',
                                                  'build/sections/conclusion.tex']))
    (HERE/'mechanism_results.json').write_text(json.dumps(serial(payload),indent=2)+'\n')
    print(json.dumps(dict(status='passed',symbolic_identity=True,
                          rational_matrix_examples=5,pi_approximant_observations=10)))


if __name__=='__main__':main()
