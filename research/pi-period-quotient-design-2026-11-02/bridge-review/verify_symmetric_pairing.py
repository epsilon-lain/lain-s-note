"""Independent exact audit of squared-period coordinates and defects."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json


def add(*rows):
    out=[F(0)]*max(map(len,rows))
    for row in rows:
        for i,x in enumerate(row):out[i]+=x
    while len(out)>1 and out[-1]==0:out.pop()
    return out


def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    while len(out)>1 and out[-1]==0:out.pop()
    return out


def scale(a,c):return [c*x for x in a]
def derivative(a):return [j*a[j] for j in range(1,len(a))] or [F(0)]


def invariant(y,ab):
    dy=derivative(y);ddy=derivative(dy)
    return add(mul([0,2,-2],mul(y,ddy)),
               mul([0,-1,1],mul(dy,dy)),
               mul([2,-3],mul(y,dy)),scale(mul(y,y),-4*ab))


def coefficients(a,b,n):
    out=[F(1)]
    for k in range(n):out.append(out[-1]*(k+a)*(k+b)/(k+1)**2)
    return out


def valuation(x,p):
    count=0
    while x%p==0:x//=p;count+=1
    return count


def digit_sum(n,b):
    total=0
    while n:total+=n%b;n//=b
    return total


def run():
    highest=30;integerchecks=clausen=defects=0
    aa=coefficients(F(1,12),F(5,12),highest)
    uu=[v*1728**n for n,v in enumerate(aa)]
    for n,u in enumerate(uu):
        assert u.denominator==1
        assert valuation(u.numerator,2)==2*digit_sum(n,2)
        assert valuation(u.numerator,3)==digit_sum(n,3)
        integerchecks+=1
    ff=coefficients(F(1,4),F(1,4),highest)
    ww=[v*64**n for n,v in enumerate(ff)]
    assert all(x.denominator==1 for x in ww)
    for n in range(highest+1):
        assert sum(ww[i]*ww[n-i] for i in range(n+1))==comb(2*n,n)**3
        clausen+=1
    for a,b in [(F(1,4),F(1,4)),(F(1,12),F(5,12))]:
        cc=coefficients(a,b,15)
        for n in range(16):
            f=cc[:n+1];y=mul(f,f)
            left=invariant(y,a*b)
            source=[F(0)]*n+[-4*(n+a)*(n+b)*cc[n]]
            right=mul(source,mul(mul(f,f),f))
            assert add(left,scale(right,-1))==[F(0)]
            defects+=1
    out=dict(status='passed',integer_level_one_terms=integerchecks,
             integer_clausen_convolutions=clausen,
             exact_squared_truncation_defects=defects,
             period_normalization='Canonical log companion has W=1/(z sqrt(1-z)); physical companion introduces1/(2pi i).',
             scope='Exact bounded algebra audits; structural invariants alone supply no denominator or pi-quality theorem.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':run()
