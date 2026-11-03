#!/usr/bin/env python3
"""Exact checks of a bounded real address space and its rational grid layers."""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


def phi(u):
    assert -1<u<1
    return u/(1-abs(u))


def address(x): return x/(1+abs(x))


def main():
    N=100
    grid={};gcd_checks=0
    for n in range(1,N+1):
        for a in range(1-n,n):
            y=phi(F(a,n));g=gcd(a,n)
            assert y==F(a,n-abs(a))
            assert gcd(a,n-abs(a))==g
            assert y.denominator==(n-abs(a))//g
            assert address(y)==F(a,n)
            grid.setdefault(y,n)
            gcd_checks+=1
    target=set()
    for q in range(1,N+1):
        for p in range(q-N,N-q+1):
            if gcd(p,q)==1 and abs(p)+q<=N:
                y=F(p,q);target.add(y)
                assert grid[y]==q+abs(p)
    assert set(grid)==target
    transfer_checks=0
    for X in [F(22,7),F(355,113),F(-22,7),F(-1,3)]:
        v=address(X)
        for n in range(2,61):
            for a in range(1-n,n):
                if a==0 or (a>0)!=(v>0): continue
                u=F(a,n);P=phi(u);q=P.denominator;B=u.denominator
                assert abs(P-X)==abs(u-v)/((1-abs(u))*(1-abs(v)))
                assert q*q*abs(P-X)==B*B*abs(u-v)*(1-abs(u))/(1-abs(v))
                transfer_checks+=1
    examples=[]
    for X in [F(22,7),F(355,113),F(-22,7)]:
        u=address(X)
        examples.append({'value':str(X),'address':str(u),'minimal_layer':u.denominator})
    out={'status':'passed','complete_layers_through':N,'grid_gcd_checks':gcd_checks,
         'distinct_rational_outputs':len(grid),'quality_transfer_checks':transfer_checks,
         'examples':examples,'scope':'Exact finite checks; full real coverage uses infinite-precision addresses.'}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__': main()
