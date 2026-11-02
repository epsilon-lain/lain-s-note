"""Independent exact audit of moment cancellation and period equivalence."""
from fractions import Fraction as F
from math import comb, factorial, gcd
from pathlib import Path
import json


def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return out


def endpoint(a,b):
    return [F(0)]*a+[F((-1)**j*comb(b,j)) for j in range(b+1)]


def rising(a,k):
    x=F(1)
    for j in range(k): x*=a+j
    return x


def rod(a,b,k):
    out=[F(0)]*(k+1)
    factor=F(4**k,factorial(k))
    for j in range(k+1):
        coef=factor*(-1)**j*comb(k,j)*rising(F(2*a+1,2)+j,k-j)*rising(b+k-j+1,j)
        for l in range(k-j+1): out[j+l]+=coef*(-1)**l*comb(k-j,l)
    assert all(x.denominator==1 for x in out)
    return [int(x) for x in out]


def moment(n,j=0):
    return sum((F(2*x,2*(l+j)+1) for l,x in enumerate(n)),F(0))


def period(n,c):
    assert len(n)<=c+1
    d=[sum((n[j]*comb(j,l)*(-1)**(j-l) for j in range(l,len(n))),F(0)) for l in range(c+1)]
    u,v=F(0),F(1)
    beta,alpha=d[c],F(0)
    for s in range(1,c+1):
        if s>=2:
            u=F(2*s-3,2*s-2)*u+F(1,(s-1)*2**s)
            v=F(2*s-3,2*s-2)*v
        beta+=d[c-s]*u
        alpha+=d[c-s]*v
    return beta,alpha


def run():
    rows=[];checks=0;zero_alpha=[]
    for m in range(1,9):
        a,b,c=m,2*m,4*m+1
        for k in range(m+2):
            h=rod(a,b,k)
            n=mul(endpoint(a,b),h)
            for j in range(k): assert moment(n,j)==0
            old=period(n,c)
            new=period(endpoint(a+k,b+k),c+k)
            scalar=4**k*comb(c+k-1,k)
            assert old==tuple(scalar*x for x in new)
            alpha=old[1]
            if not alpha:
                zero_alpha.append([m,k])
                continue
            r=-4*old[0]/old[1]
            rr=-4*new[0]/new[1]
            assert r==rr
            checks+=1
            if k==m+1:
                rows.append(dict(m=m,k=k,h=h,alpha=str(alpha),p=str(r.numerator),q=str(r.denominator),scalar=str(scalar)))
    out=dict(status='passed',maximum_m=8,period_and_reduced_fraction_checks=checks,zero_alpha_cases=zero_alpha,maximal_cancellation_samples=rows,scope='Exact finite algebra audit; no quality-zero limit claimed.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='maximal_cancellation_samples'},indent=2))


if __name__=='__main__':run()
