"""Exact checks for the integer identities found in Zuuv v2 source.
This does not execute Zuuv or reproduce the 2025 computation.
"""
from fractions import Fraction as F
from random import Random
import json
from pathlib import Path

def matrix(coeffs):
    p0,p1,q0,q1=0,1,1,0
    for a in coeffs:
        p0,p1=p1,a*p1+p0
        q0,q1=q1,a*q1+q0
    return p1,p0,q1,q0

rng=Random(20261102)
count=0
for _ in range(3000):
    coeffs=[rng.randrange(0,7)]+[rng.randrange(1,30) for _ in range(rng.randrange(1,18))]
    p,pp,q,qp=matrix(coeffs)
    assert abs(p*qp-pp*q)==1
    t=F(rng.randrange(1,10000),rng.randrange(1,1000))
    x=(p*t+pp)/(q*t+qp)
    N,D=x.numerator,x.denominator
    cnum=D*pp-N*qp
    cden=N*q-D*p
    assert F(cnum,cden)==t
    count+=1

# Actual sign-only acceptance predicate from the disk-mode source.
p,pp,q,qp=matrix([1,1]);N,D=7,5
cnum,cden=D*pp-N*qp,N*q-D*p
assert (cnum,cden)==(-2,-3)
if cnum<0 and cden<0:cnum,cden=-cnum,-cden
assert cnum>0 and cden>0 and cnum<cden
assert F(cnum,cden)==F(2,3)
assert F(N,D)==F(7,5)

out={'matrix_correction_identity_checks':count,
     'sign_only_does_not_certify_prefix':{
         'input_fraction':[7,5], 'guessed_coefficients':[1,1],
         'raw_correction':[-2,-3], 'after_same_sign_normalization':[2,3],
         'residual_tail':'2/3 < 1',
         'true_regular_cf':[1,2,2],
         'scope':'Predicate counterexample only; not a demonstrated execution counterexample after the program truncates ten terms.'},
     'full_record_recomputed':False}
Path(__file__).with_name('static_formula_checks.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
