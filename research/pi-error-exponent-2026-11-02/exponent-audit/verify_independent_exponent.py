"""Minimal independent exact checks for the fixed-family exponent audit.

No floating-point powers or logarithms. All-index implications are established
in the written proof audit, not extrapolated from this finite run.
"""
from fractions import Fraction as F
from math import factorial,gcd,comb,isqrt
from pathlib import Path
import json


def v2(x):
    if isinstance(x,F):
        return v2(x.numerator)-v2(x.denominator)
    assert x
    x=abs(x)
    return (x&-x).bit_length()-1


def valuation_formula(n):
    m,r=divmod(n,4)
    return (2*m,2*m,2*m+2,2*m+3+v2(m+1))[r]


def multiply(x,y):
    a,b=x;c,d=y
    return a*c+2*b*d,a*d+b*c


def scale(x,s):
    return x[0]*s,x[1]*s


def sign(x):
    a,b=x
    if not b:
        return (a>0)-(a<0)
    if not a:
        return (b>0)-(b<0)
    if a>0 and b>0:
        return 1
    if a<0 and b<0:
        return -1
    difference=a*a-2*b*b
    out=(difference>0)-(difference<0)
    return out if a>0 else -out


def lambda_power(k):
    base=(F(1),F(1)) if k>=0 else (F(-1),F(1))
    k=abs(k)
    out=(F(1),F(0))
    while k:
        if k&1:
            out=multiply(out,base)
        base=multiply(base,base)
        k//=2
    return out


def main():
    maximum=500
    A=[1,1];B=[0,1]
    for n in range(1,maximum):
        A.append((2*n+1)*A[-1]+n*n*A[-2])
        B.append((2*n+1)*B[-1]+n*n*B[-2])
    negative=[]
    lcm=1
    lcms=[1]
    uniform_checks=0
    exponents=tuple(map(F,('1/2','1','4/3','3/2','2','5/2','3','4')))
    for n in range(maximum+1):
        g=gcd(A[n],4*B[n])
        assert v2(g)==valuation_formula(n)-(n==3)
        q=A[n]//g
        assert (q%2==0)==(n==3)
        if not n:
            continue
        r=F(4*B[n],A[n])
        previous=F(4*B[n-1],A[n-1])
        actual=v2(r-previous)
        expected=2+2*v2(factorial(n-1))-valuation_formula(n)-valuation_formula(n-1)
        m,residue=divmod(n,4)
        piecewise=(4*m-3-2*m.bit_count()-3*v2(m) if residue==0 else
                   4*m+2-2*m.bit_count() if residue==1 else
                   4*m-2*m.bit_count() if residue==2 else
                   4*m-1-2*m.bit_count()-v2(m+1))
        assert actual==expected==piecewise
        if actual<0:
            negative.append(n)
        if n>=5:
            assert actual>=0
        nf=factorial(n)
        ell=F(A[n],nf);R=F(g,nf)
        d=v2(R)
        E=R.numerator//2**v2(R.numerator)
        D=R.denominator//2**v2(R.denominator)
        assert R==F(2)**d*F(E,D)
        assert d==valuation_formula(n)-v2(nf)-(n==3)
        beta_twice=2*(valuation_formula(n)-v2(nf))+n
        beta_pieces=(2*n.bit_count() if residue==0 else
                     2*n.bit_count()-1 if residue==1 else
                     2*n.bit_count()+2 if residue==2 else
                     2*n.bit_count()+2*v2(n+1)-1)
        assert beta_twice==beta_pieces
        assert 2*d+n>=0
        lcm=lcm*n//gcd(lcm,n)
        lcms.append(lcm)
        odd_lcm=lcm//2**v2(lcm)
        assert odd_lcm%D==0
        assert R*R*2**n*lcm*lcm>=1
        if n>40:
            continue
        explicit=sum(F(nf,2**j*factorial(j)**2*factorial(n-2*j))
                     for j in range(n//2+1))
        assert ell==explicit
        lam=lambda_power(n)
        assert sign((lam[0]-ell,lam[1]))>=0
        assert sign(((2*n+1)*ell-lam[0],-lam[1]))>=0
        for c in exponents:
            k,l=c.numerator,c.denominator
            power=lambda_power((k-2*l)*n)
            # Raise each uniform bound to l, making all powers integral.
            lower=scale(power,F(2)**l/(R**k*(2*n+1)**(l+max(k-2*l,0))))
            upper=scale(power,F(4)**l*F(2*n+1)**(max(2*l-k,0)-l)/R**k)
            analytic_lower=F(2)**l*ell**(k-2*l)/((2*n+1)**l*R**k)
            analytic_upper=F(4)**l*ell**(k-2*l)/((2*n+1)**l*R**k)
            assert sign((analytic_lower-lower[0],-lower[1]))>=0
            assert sign((upper[0]-analytic_upper,upper[1]))>=0
            uniform_checks+=1
    assert negative==[3,4]
    halving_checks=0
    for m in range(1,maximum//2+1):
        assert comb(2*m,m)%(lcms[2*m]//lcms[m])==0
        halving_checks+=1
        if 2*m+1<=maximum:
            assert ((2*m+1)*comb(2*m,m))%(lcms[2*m+1]//lcms[m])==0
            halving_checks+=1
    # Independent construction of the n=50000 necessary-count certificate.
    # No values of A or B at this index are needed or computed.
    index=50000
    def valuation(x,p):
        count=0
        while x%p==0:
            x//=p
            count+=1
        return count
    def floor_log(x,p):
        power=p
        count=0
        while power<=x:
            count+=1
            power*=p
        return count
    small=[p for p in range(3,isqrt(index-1)+1,2)
           if all(p%d for d in range(2,isqrt(p)+1))]
    small_budget=1
    for p in small:
        cap=max(0,min(floor_log(index-1,p)-valuation(index,p),
                      floor_log(index+1,p)-valuation(index+1,p)))
        small_budget*=p**cap
    dA=valuation_formula(index)-(index-index.bit_count())
    assert dA==-24994
    right=2**(1-2*dA)
    base=3*(2*index+1)
    assert base*(index-1)**3200<=right<base*(index-1)**3202
    assert base*small_budget**2*(index-1)**3120<=right<base*small_budget**2*(index-1)**3122
    assert len(small)==47 and len(str(small_budget))==188
    result={'status':'passed','maximum_index':maximum,
            'actual_gcd_saturation_checks':maximum+1,
            'rational_difference_valuation_checks':maximum,
            'negative_difference_indices':negative,
            'net_odd_and_binary_factorizations':maximum,
            'odd_deficit_lcm_divisibility_checks':maximum,
            'exact_weak_lower_gcd_bounds':maximum,
            'elementary_lcm_halving_divisibilities':halving_checks,
            'uniform_ell_bounds_in_Qsqrt2':40,
            'uniform_real_exponent_bound_checks':uniform_checks,
            'independent_n50000_certificate':{
                'small_odd_primes':len(small),'small_budget_decimal_digits':len(str(small_budget)),
                'dA':dA,'required_power_of_two':1-2*dA,
                'necessary_positive_odd_prime_count':1601,
                'necessary_large_positive_odd_prime_count':1561,
                'previous_integer_fails_each_necessary_test':True},
            'scope':'Finite exact checks only; the written audit supplies all-index proofs, asymptotic constants, and rate-limit boundary arguments.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
