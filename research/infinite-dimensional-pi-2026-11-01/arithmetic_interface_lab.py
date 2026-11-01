#!/usr/bin/env python3
"""Exact checks of coherent infinite embeddings and finite arithmetic costs.

No inference about pi's infinite Diophantine approximation tail is made.
"""
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


def chord_squared(s,t):
    x,y=s/2,t/2
    return x*x/(1-x*x)+y*y/(1-y*y)-2*x*y/(1-x*y)


def prefix_squared(s,t,m):
    return sum(((s/2)**k-(t/2)**k)**2 for k in range(1,m+1))


def derivative_tail_squared(r,m):
    u=r*r/4
    return u**m/4*((m+1)**2/(1-u)+2*(m+1)*u/(1-u)**2+u*(1+u)/(1-u)**3)


def nearest_integer_distance(x):
    k=x.numerator//x.denominator
    return min(x-k,k+1-x)


def observed_distance_squared(s,t,m):
    return sum(((nearest_integer_distance(q*s)-nearest_integer_distance(q*t))/4**q)**2
               for q in range(1,m+1))


def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x


def main():
    exact_chords=0;tail_checks=0;series_checks=0
    tail_records=[]
    for r in (F(1,2),F(1),F(3,2)):
        u=r*r/4
        full_derivative=(1+u)/(4*(1-u)**3)
        points=[r*F(j,4) for j in range(-4,5)]
        for s in points:
            for t in points:
                value=chord_squared(s,t)
                d=(s-t)**2
                assert d/4<=value<=full_derivative*d
                exact_chords+=1
                for m in (1,2,4,8):
                    prefix=prefix_squared(s,t,m)
                    assert 0<=value-prefix<=derivative_tail_squared(r,m)*d
                    tail_checks+=1
        for m in (0,1,2,4,8):
            direct=sum(F(k*k)*u**(k-1)/4 for k in range(1,m+1))
            assert direct+derivative_tail_squared(r,m)==full_derivative
            series_checks+=1
            tail_records.append(dict(r=r,m=m,C_squared=full_derivative,
                                     relative_tail_constant_squared=4*derivative_tail_squared(r,m)))
    denominator_checks=[]
    for candidate in (F(22,7),F(355,113),F(3),F(10,3)):
        base=(candidate-3)/2
        for m in (1,2,3,5,8):
            coords=[base**k for k in range(1,m+1)]
            actual=lcm(*(x.denominator for x in coords))
            assert actual==base.denominator**m
            assert all((actual*x).denominator==1 for x in coords)
            denominator_checks.append(dict(candidate=candidate,m=m,
                                           base=base,minimum_common_denominator=actual))
    reflection=0;observation_bounds=0
    lipschitz_squared=F(272,3375)
    for s in (F(1,3),F(2,3),F(2,7),F(1,100),F(-3,4)):
        for q in range(1,31):
            assert nearest_integer_distance(q*s)==nearest_integer_distance(q*(1-s))
            reflection+=1
        for t in (F(1,3),F(2,3),F(3,7),F(0),F(7,5)):
            for m in (1,3,8,16):
                extra=observed_distance_squared(s,t,m)
                assert 0<=extra<=lipschitz_squared*(s-t)**2
                observation_bounds+=1
    # One nonzero observation coordinate tends to zero uniformly for every x.
    artificial_smallness=[dict(q=q,universal_coordinate_upper=F(1,2*4**q),
                               q_squared_times_upper=F(q*q,2*4**q))
                          for q in (1,2,4,8,16,32)]
    payload=dict(status='passed',scope=__doc__,exact_chord_checks=exact_chords,
                 relative_tail_checks=tail_checks,exact_series_splits=series_checks,
                 tail_constants=tail_records,denominator_checks=denominator_checks,
                 reflection_checks=reflection,observation_lipschitz_checks=observation_bounds,
                 all_denominator_lipschitz_squared=lipschitz_squared,
                 weighting_alone_produces_small_coordinates=artificial_smallness,
                 pi_badly_approximable_conclusion='not_established')
    (HERE/'arithmetic_interface_results.json').write_text(json.dumps(serial(payload),indent=2)+'\n')
    print(json.dumps(dict(status='passed',exact_chords=exact_chords,
                          relative_tail_checks=tail_checks,denominator_cases=len(denominator_checks),
                          observation_bounds=observation_bounds)))


if __name__=='__main__':main()
