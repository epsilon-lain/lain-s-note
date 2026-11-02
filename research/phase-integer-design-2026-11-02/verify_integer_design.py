"""Exact checks of quotient denominators and bounded congruence designs."""
from fractions import Fraction as F
from math import gcd
from itertools import product
from pathlib import Path
import json

def egcd(a,b):
    if b==0:return abs(a),1 if a>=0 else -1,0
    g,x,y=egcd(b,a%b)
    return g,y,x-(a//b)*y

def bezout_all(values):
    g=values[0];coeff=[1]
    for value in values[1:]:
        gg,x,y=egcd(g,value)
        coeff=[x*c for c in coeff]+[y]
        g=gg
    return g,coeff

def main():
    fee_cases=0
    targets=[F(0),F(22,7),F(-5,3)]
    for A,B,C in product(range(-5,6),repeat=3):
        Delta=A-2*B+C
        if not Delta:continue
        x=B-A;d=gcd(abs(Delta),x*x);delta=Delta//d;k=x*x//d
        assert gcd(abs(delta),k)==1
        for D in range(1,9):
            R=F(A*C-B*B,D*Delta)
            h=gcd(D,abs(A*delta-k))
            assert R.denominator==D*abs(delta)//h
            if x*x%Delta==0:assert D%R.denominator==0
            for alpha in targets:
                assert R.denominator**2*abs(R-alpha)==F(D*abs(delta),h*h)*abs(delta*(A-D*alpha)-k)
            fee_cases+=1
    congruence_cases=0
    for m in range(1,13):
        for B1,B2 in product(range(-3,4),repeat=2):
            d,co=bezout_all([m,B1,B2])
            period=m//d
            assert sum(a*b for a,b in zip(co,[m,B1,B2]))==d
            image={(B1*c1+B2*c2)%m for c1,c2 in product(range(m),repeat=2)}
            assert len(image)==period
            assert image==set(range(0,m,d))
            for A in range(-3,4):
                solvable=(-A)%m in image
                assert solvable==(A%d==0)
                if solvable:
                    raw=[-(A//d)*v for v in co[1:]]
                    centered=[(v+period//2)%period-period//2 for v in raw]
                    assert all(abs(v)<=period//2 for v in centered)
                    assert (A+B1*centered[0]+B2*centered[1])%m==0
                congruence_cases+=1
    report={'status':'PASS','quotient_fee_cases':fee_cases,
            'targets_per_quotient':len(targets),'bounded_congruence_cases':congruence_cases,
            'checks':['exact_reduced_denominator','error_and_denominator_joint_identity',
                      'necessary_and_sufficient_modular_solvability','lattice_index','bounded_coefficient_solution']}
    Path(__file__).with_name('integer_design_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
