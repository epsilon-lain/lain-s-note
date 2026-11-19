#!/usr/bin/env python3
"""Fixed weights; legal polynomial columns Y^h X2^a (X1-X2)^b; raw full-log jets."""
from fractions import Fraction as F
from math import comb

from verify_two_nodes import determinant, multiply_series


def raw_fixed(height, r1, r2):
    assert height >= 1
    delta = r1-r2
    orders = [(height-2*s-1)//4 for s in range((height+1)//2)]
    cols = [(h,a,b) for h,B in enumerate(orders) for a in range(B+1)
            for b in range(2*(B-a+1))]
    assert all(h+a+b <= height for h,a,b in cols)
    orders_t = len(orders)
    log = [F(0)] + [F((-1)**(s+1),s) for s in range(1,orders_t)]
    rows = []
    for j in (0,1):
        x2 = [j*r2]+log[1:]
        powers = [[F(1)] + [F(0)]*(orders_t-1)]
        for _ in range(orders[0]):
            powers.append(multiply_series(powers[-1],x2))
        for s,B in enumerate(orders):
            for i in range(B+1):
                for k in range(B+1-i):
                    def term(h,a,b):
                        out = F(0)
                        for u in range(min(k,a)+1):
                            if i+k-u > b: continue
                            dcoeff = (comb(b,i)*comb(b-i,k-u)*(-1)**(k-u)
                                      *(j*delta)**(b-i-k+u))
                            tx = sum(comb(h,v)*powers[a-u][s-v]
                                     for v in range(min(h,s)+1))
                            out += comb(a,u) * dcoeff * tx
                        return out
                    rows.append([term(*col) for col in cols])
    exponent = sum(sum(i*i for i in range(1,B+2)) for B in orders)
    assert len(rows)==len(cols)
    return len(rows), exponent, abs(determinant(rows))

if __name__=='__main__':
    for h in range(1,12):
        for r1,r2 in [(F(22,7),F(355,113)),(F(2,5),F(9,11))]:
            size,e,d=raw_fixed(h,r1,r2)
            assert d == abs(r1-r2)**e
        print(f'H={h} M={size} E={e}')
