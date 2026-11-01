#!/usr/bin/env python3
"""Exact finite audits of a genuinely infinite-dimensional feature curve.

The infinite assertions are proved in kernel_geometry.md. These checks neither
construct new pi approximants nor establish the pi017 theorem.
"""
from fractions import Fraction as F
from itertools import combinations
from math import factorial
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


def determinant(matrix):
    a=[list(map(F,row)) for row in matrix]
    sign=1;out=F(1)
    for k in range(len(a)):
        pivot=next((i for i in range(k,len(a)) if a[i][k]),None)
        if pivot is None:return F(0)
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];sign=-sign
        v=a[k][k];out*=v
        for i in range(k+1,len(a)):
            c=a[i][k]/v
            for j in range(k+1,len(a)):a[i][j]-=c*a[k][j]
            a[i][k]=F(0)
    return sign*out


def kernel(x,y):return 1/(1-x*y)


def gram(xs,L=None):
    return [[kernel(x,y) if L is None else sum(((x*y)**k for k in range(L)),F(0))
             for y in xs] for x in xs]


def closed_determinant(xs):
    numerator=F(1);denominator=F(1)
    for x,y in combinations(xs,2):numerator*=(x-y)**2
    for x in xs:
        for y in xs:denominator*=1-x*y
    return numerator/denominator


def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [serial(v) for v in x]
    return x


def main():
    matrices=[];tails=[]
    for N in range(1,7):
        for centered in (False,True):
            xs=[F(2*k-N+1,2*N+3) if centered else F(k+1,2*N+3)
                for k in range(N)]
            actual=determinant(gram(xs));predicted=closed_determinant(xs)
            assert actual==predicted and actual>0
            matrices.append(dict(N=N,xs=xs,determinant=actual))
            r=max(map(abs,xs));B=1/(1-r*r)
            previous=F(0)
            for L in (max(1,N-1),N,N+1,2*N+3):
                part=determinant(gram(xs,L))
                assert previous<=part<=actual
                bound=factorial(N)*N*r**(2*L)*B**N
                assert actual-part<=bound
                if L<N:assert part==0
                tails.append(dict(N=N,L=L,determinant_gap=actual-part,
                                  explicit_absolute_error_bound=bound))
                previous=part
    repeat=[F(1,3),F(1,3),F(2,5)]
    assert determinant(gram(repeat))==closed_determinant(repeat)==0
    scalar=[]
    for a in range(-4,5):
        for b in range(-4,5):
            x,y=F(a,6),F(b,6)
            norm2=kernel(x,x)+kernel(y,y)-2*kernel(x,y)
            expression=(x-y)**2*(1+x*y)/((1-x*x)*(1-y*y)*(1-x*y))
            assert norm2==expression
            r=max(abs(x),abs(y))
            # x=t/2; norm2 >= (s-t)^2/4 = (x-y)^2.
            lower=(x-y)**2
            upper=(x-y)**2*(1+r*r)/(1-r*r)**3
            assert lower<=norm2<=upper
            scalar.append(dict(x=x,y=y,norm_squared=norm2,lower=lower,upper=upper))
    families=[]
    for Q in (2,3,4,5,10,100,1000):
        eps=F(1,Q);xs=[F(0),eps/2,eps]
        D=(4*Q*Q-1)*(Q*Q-1)*(2*Q*Q-1)**2
        expected=F(Q*Q,D)
        actual=determinant(gram(xs))
        assert actual==expected
        assert actual.denominator==D and actual.numerator==Q*Q
        ratio=16*Q**6*actual
        assert 1<=ratio<=1/(1-F(9,4*Q*Q))
        families.append(dict(Q=Q,epsilon=eps,determinant=actual,
                             determinant_times_clearing_denominator=actual*D,
                             exact_reduced_denominator=D,
                             ratio_to_epsilon6_over16=ratio,
                             normalized_volume_squared_limit=1))
    result=dict(status='passed',scope=__doc__,
                kernel_identity_checks=len(matrices),repeated_node_check=True,
                determinant_tail_checks=len(tails),scalar_bilipschitz_checks=len(scalar),
                matrices=matrices,tail_bounds=tails,
                scalar_distance_examples=scalar,
                rational_node_families=families,
                pi_claim='None. Infinite feature dimension does not by itself produce stronger scalar rational approximation.')
    (HERE/'kernel_results.json').write_text(json.dumps(serial(result),indent=2)+'\n')
    print(json.dumps(dict(status='passed',kernel_identity_checks=len(matrices),
                          determinant_tail_checks=len(tails),scalar_distance_checks=len(scalar))))


if __name__=='__main__':main()
