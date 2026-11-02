"""Independent exact audits of the joint Legendre/arctangent arithmetic notes."""
from fractions import Fraction as F
from math import factorial, gcd
from pathlib import Path
import json


def trim(a):
    a = list(a)
    while len(a)>1 and not a[-1]:
        a.pop()
    return a


def add(*args):
    a = [0]*max(map(len,args))
    for b in args:
        for j,x in enumerate(b):
            a[j] += x
    return trim(a)


def scale(a,c):
    return trim([c*x for x in a])


def multiply(a,b):
    result = [0]*(len(a)+len(b)-1)
    for j,x in enumerate(a):
        for k,y in enumerate(b):
            result[j+k] += x*y
    return trim(result)


def derivative(a):
    return trim([j*a[j] for j in range(1,len(a))] or [0])


def telescoping_identity():
    # r=N/d, v=V/d. Separately compare the coefficients of n and 1
    # in d(v*r^n)/dx after extracting r^(n-1)/d and clearing d^2.
    d,N,V = [1,0,1],[0,1,-1],[-1,2,1]
    r_derivative_numerator = add(multiply(derivative(N),d),
                                scale(multiply(N,derivative(d)),-1))
    v_derivative_numerator = add(multiply(derivative(V),d),
                                scale(multiply(V,derivative(d)),-1))
    nn,nd,dd = multiply(N,N),multiply(N,d),multiply(d,d)
    assert multiply(V,r_derivative_numerator) == add(scale(nn,4),scale(nd,4),scale(dd,-1))
    assert multiply(v_derivative_numerator,N) == add(scale(nn,4),scale(nd,2))
    return 2


def atan_interval(x, terms=600):
    assert 0 < x < 1
    result, power = F(0),x
    for k in range(terms):
        result += (-1)**k*power/(2*k+1)
        power *= x*x
    next_value = result+(-1)**terms*power/(2*terms+1)
    return min(result,next_value),max(result,next_value)


def pi_interval():
    a,b=atan_interval(F(1,5))
    c,d=atan_interval(F(1,239))
    return 16*a-4*d,16*b-4*c


def target_interval(a,c,pi_bounds):
    lo,hi=pi_bounds
    if a == 1:
        return c*lo/4,c*hi/4
    if a > 1:
        al,ah=atan_interval(1/a)
    else:
        al,ah=atan_interval(a)
        al,ah=lo/2-ah,hi/2-al
    return c*al,c*ah


def family(a,c,limit,pi_bounds):
    ell,w=[F(1),a],[F(0),F(1)]
    for n in range(1,limit):
        ell.append(((2*n+1)*a*ell[n]+n*ell[n-1])/(n+1))
        w.append(((2*n+1)*a*w[n]+n*w[n-1])/(n+1))
    tl,th=target_interval(a,c,pi_bounds)
    fractions,quality_bounds=[],[]
    checks=0
    for n in range(limit+1):
        value=c*w[n]/ell[n]
        p,q=value.numerator,value.denominator
        assert q>0 and gcd(p,q)==1
        fractions.append(value)
        if n%2==0:
            low,high=tl-value,th-value
        else:
            low,high=value-th,value-tl
        assert 0<low<high
        strict_lower=c*a/((a*a+1)*(2*n+1)*ell[n]**2)
        strict_upper=c/(a*(2*n+1)*ell[n]**2)
        assert strict_lower<low<high<strict_upper
        quality_bounds.append((q*q*low,q*q*high))
        checks+=1
        if n:
            cross=ell[n]*w[n-1]-ell[n-1]*w[n]
            assert cross==F((-1)**n,n)
            prior=fractions[n-1]
            assert abs(value-prior)==c/(n*ell[n]*ell[n-1])
            assert q*prior.denominator>=n*ell[n]*ell[n-1]/c
            product_lower=a*a*n*n/((a*a+1)**2*(4*n*n-1))
            assert quality_bounds[n][0]*quality_bounds[n-1][0]>product_lower
            checks+=4
        if n>=2:
            prior=fractions[n-2]
            separation=c*(2*n-1)*a/(n*(n-1)*ell[n]*ell[n-2])
            assert abs(value-prior)==separation
            assert q*prior.denominator>=1/separation
            product_lower=F(n*n*(n-1)**2,(2*n-1)**2*(2*n+1)*(2*n-3))/(a*a+1)**2
            assert quality_bounds[n][0]*quality_bounds[n-2][0]>product_lower
            checks+=3
    return checks


def integer_checks(limit,pi_bounds):
    aa,bb=[1,1],[0,1]
    for n in range(1,limit):
        aa.append((2*n+1)*aa[n]+n*n*aa[n-1])
        bb.append((2*n+1)*bb[n]+n*n*bb[n-1])
    gg=[gcd(a,4*b) for a,b in zip(aa,bb)]
    checks=0
    lo,hi=pi_bounds
    for n in range(limit+1):
        r=F(4*bb[n],aa[n])
        assert r.denominator==aa[n]//gg[n]
        assert r.numerator==4*bb[n]//gg[n]
        ell=F(aa[n],factorial(n))
        assert r.denominator/ell==F(factorial(n),gg[n])
        # Exact interval equality for the two ways of computing joint cost.
        error=(lo-r,hi-r) if n%2==0 else (r-hi,r-lo)
        assert error[0]>0
        for e in error:
            assert r.denominator**2*e==F(factorial(n),gg[n])**2*ell**2*e
        checks+=4
        if n:
            assert aa[n]*bb[n-1]-aa[n-1]*bb[n]==(-1)**n*factorial(n-1)**2
            assert (4*factorial(n-1)**2)%(gg[n]*gg[n-1])==0
            if n>=2:
                # Square of h_n*h_(n-1)<=4 sqrt((n-1)/n).
                hsquare=F((gg[n]*gg[n-1])**2*n*(n-1),
                          (factorial(n)*factorial(n-1))**2)
                assert hsquare<=F(16*(n-1),n)
            checks+=3
    return checks


def main():
    bounds=pi_interval()
    assert 3<bounds[0]<bounds[1]<F(22,7)
    assert bounds[1]-bounds[0]<F(1,10**800)
    cases=[(F(1),F(4),100),(F(5),F(1),100),
           (F(239),F(1),100),(F(1,2),F(7,3),40)]
    result={
        'status':'passed',
        'all_n_telescoping_polynomial_identities':telescoping_identity(),
        'families':[{'a':str(a),'c':str(c),'max_n':n,
                     'exact_checks':family(a,c,n,bounds)} for a,c,n in cases],
        'integer_and_gcd_checks_through_n_100':integer_checks(100,bounds),
        'pi_interval_width_less_than':'10^(-800)',
        'scope':'Exact algebra and finite inequalities only. Laplace constants '
                'and the general all-n claims are independently derived in '
                'INDEPENDENT_JOINT_ARITHMETIC_REVIEW.md. No subsequence or '
                'global pi bad-approximability claim.',
    }
    Path(__file__).with_name('independent_joint_arithmetic_checks.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
