"""Independent exact checks for LOCAL_POLES_AND_EXCESS_SIEVE.md."""
from fractions import Fraction as F
from math import gcd, lcm, factorial
from pathlib import Path
import json


def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:
        p.pop()
    return p


def add(*polys):
    out=[F(0)]*max(map(len,polys))
    for p in polys:
        for i,c in enumerate(p):
            out[i]+=c
    return trim(out)


def scale(p,c):
    return trim([x*c for x in p])


def mul(p,q):
    out=[F(0)]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        for j,y in enumerate(q):
            out[i+j]+=x*y
    return trim(out)


def deriv(p):
    return trim([i*p[i] for i in range(1,len(p))] or [F(0)])


def polynomial_checks(maximum=16):
    ell,w=[[F(1)],[F(0),F(1)]],[[F(0)],[F(1)]]
    for n in range(1,maximum):
        for seq in (ell,w):
            seq.append(scale(add(scale([F(0)]+seq[n],2*n+1),
                                 scale(seq[n-1],n)),F(1,n+1)))
    count=0
    for n in range(1,maximum+1):
        assert mul([1,0,1],deriv(ell[n]))==scale(add([F(0)]+ell[n],ell[n-1]),n)
        assert mul([1,0,1],deriv(w[n]))==add(scale(add([F(0)]+w[n],w[n-1]),n),scale(ell[n],-1))
        left=mul([1,0,1],add(mul(ell[n],deriv(w[n])),scale(mul(deriv(ell[n]),w[n]),-1)))
        right=add([F((-1)**n)],scale(mul(ell[n],ell[n]),-1))
        assert left==right
        count+=3
    return count


def v2(x):
    assert x>0
    return (x & -x).bit_length()-1


def odd(x):
    return x>>v2(x)


def two_power(k):
    return F(2**k) if k>=0 else F(1,2**(-k))


def sieve_checks(limit=1000):
    u0,u1=1,2
    a0,a1=1,1
    b0,b1=0,1
    fact=1
    lprev=1
    exclusions=[]
    survivors=[]
    sharpened=[]
    counterexample=None
    for n in range(2,limit+1):
        numerator=2*(2*n-1)*u1+4*(n-1)*u0
        assert numerator%n==0
        u0,u1=u1,numerator//n
        a0,a1=a1,(2*n-1)*a1+(n-1)**2*a0
        b0,b1=b1,(2*n-1)*b1+(n-1)**2*b0
        fact*=n
        assert F(a1,fact)==F(u1,2**n)
        s=v2(u1)
        oo=odd(u1)
        lnext=lcm(lprev,n,n+1)
        ee=gcd(lprev//gcd(lprev,n),lnext//(n+1))
        gg=gcd(oo,odd(ee))
        hh=gcd(oo//gg,gg)
        iterations=0
        while True:
            nxt=gcd(gg,hh*hh)
            if nxt==hh:
                break
            hh=nxt
            iterations+=1
        refined=gg//hh
        actual_g=gcd(a1,4*b1)
        ratio=F(actual_g,fact)
        assert refined%odd(ratio.numerator)==0
        assert ratio<=two_power(s-n)*refined
        floor=2*two_power(2*(n-s))/((2*n+1)*refined**2)
        actual_analytic_floor=F(2,2*n+1)/ratio**2
        assert floor<=actual_analytic_floor
        if floor>=1:
            exclusions.append(n)
        else:
            survivors.append(n)
        if gg!=refined:
            sharpened.append({'n':n,'removed_divisor':str(hh),
                              'saturation_steps':iterations})
        if n==72:
            def vp(x,p):
                power=0
                while x%p==0:
                    x//=p
                    power+=1
                return power
            counterexample={'n':n,'p':17,'v_A':vp(a1,17),
                            'v_B':vp(b1,17),'v_factorial':vp(fact,17)}
            assert counterexample=={'n':72,'p':17,'v_A':5,'v_B':5,'v_factorial':4}
        lprev=lcm(lprev,n)
    return {'maximum_n':limit,'actual_gcd_comparisons':limit-1,
            'quality_threshold':'1',
            'indices_excluded_by_strict_sieve':len(exclusions),
            'surviving_indices':survivors,
            'whole_prime_removal_improvements':sharpened,
            'odd_excess_counterexample':counterexample}


def main():
    result={'status':'passed',
            'exact_polynomial_differential_identities':polynomial_checks(),
            **sieve_checks(),
            'scope':'The finite comparisons certify the sieve against actual '
                    'gcds. They do not establish an all-n exclusion of good '
                    'subsequences or a global statement about pi.'}
    Path(__file__).with_name('local_poles_sieve_verification.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
