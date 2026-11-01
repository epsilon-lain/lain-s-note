#!/usr/bin/env python3
"""Exact checks of rational labelled scale geometry, not a theorem about pi."""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json


def nearest_integer(x):
    x=F(x)
    return (2*x.numerator+x.denominator)//(2*x.denominator)


def observation(x,n):
    k=nearest_integer(n*x)
    residual=n*x-k
    return k,residual,abs(residual)/n


def main():
    count=0
    for x in [F(22,7),F(355,113),F(3,7),F(-5,12),F(0),F(5,2)]:
        p,q=x.numerator,x.denominator
        for n in range(1,301):
            k,r,d=observation(x,n)
            assert F(-1,2)<=r<F(1,2)
            assert r==observation(x,n+q)[1]
            assert (d==0)==(n%q==0)
            if d:
                assert d>=F(1,q*n)
            count+=1
    nonmonotone=[observation(F(3,7),n)[2] for n in (5,6,7)]
    assert nonmonotone==[F(1,35),F(1,14),F(0)]
    two_errors=[F(observation(F(3,7),n)[0],n)-F(3,7) for n in (6,10)]
    assert two_errors==[F(1,14),F(-1,35)]
    shared_checks=0
    for n1 in range(1,16):
        for n2 in range(1,16):
            g=gcd(n1,n2)
            for denominator in range(1,21):
                t=F(1,denominator)
                present=(n1*t).denominator==1 and (n2*t).denominator==1
                assert present==((g*t).denominator==1)
                shared_checks+=1
    rationals=sorted({F(p,q) for q in range(1,21) for p in range(-q,2*q+1)})
    circle_checks=0
    for c in (F(1,10),F(1,3),F(1,2)):
        for left,right in zip(rationals,rationals[1:]):
            q,b=left.denominator,right.denominator
            r,s=c/q**2,c/b**2
            dx=right-left
            center_gap=dx**2+(r-s)**2-(r+s)**2
            assert center_gap==dx**2-4*r*s
            assert center_gap>=F(1-4*c*c,q*q*b*b)
            assert center_gap>=0
            circle_checks+=1
    # Rational circle parametrization avoids floating-point roots.
    crossings=[]
    for c in (F(1,2),F(1,3)):
        for u in [F(k,20) for k in range(0,20)]:
            kappa=2*c*u/(1+u*u)
            radical=c*(1-u*u)/(1+u*u)
            H=c-radical
            assert radical>=0 and radical**2==c*c-kappa*kappa
            assert H==kappa*kappa/(c+radical)
            assert kappa*kappa/(2*c)<=H<=kappa*kappa/c
            crossings.append({'c':str(c),'kappa':str(kappa),'H':str(H)})
    result={'status':'passed','scope':__doc__,
            'rational_observation_checks':count,
            'shared_grid_checks':shared_checks,
            'adjacent_rational_disk_checks':circle_checks,
            'circle_crossing_checks':len(crossings),
            'nonmonotonic_example':{'target':'3/7','n':[5,6,7],
                                      'distances':list(map(str,nonmonotone))},
            'two_scale_errors':list(map(str,two_errors)),
            'crossing_records':crossings}
    output=Path(__file__).with_name('integer_geometry_results.json')
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='crossing_records'}))


if __name__=='__main__':
    main()
