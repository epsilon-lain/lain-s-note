"""Exact finite audits of NET_PRIME_BUDGET_AND_ERROR_EXPONENT.md.

Uses Python standard library and Q(sqrt(2)); no floating-point decision,
large-index A/B recurrence, or numerical pi/e constants are needed.
"""
from fractions import Fraction as F
from math import gcd, isqrt, lcm
from pathlib import Path
import json


def valuation(x,p):
    assert x
    count=0
    while x%p==0:
        x//=p
        count+=1
    return count


def odd(x):
    return x//(2**valuation(x,2))


def power2(k):
    return F(2**k) if k>=0 else F(1,2**(-k))


def primes_to(n):
    return [p for p in range(3,n+1,2)
            if all(p%d for d in range(2,isqrt(p)+1))]


def floor_log(n,p):
    assert n>=1
    count=0
    while n>=p:
        count+=1
        n//=p
    return count


def a_valuation_formula(n):
    k,r=divmod(n,4)
    return (2*k,2*k,2*k+2,2*k+3+valuation(k+1,2))[r]


def beta_formula(n):
    r=n%4
    s=n.bit_count()
    return (F(s),F(s)-F(1,2),F(s+1),
            F(s+valuation(n+1,2))-F(1,2))[r]


def sign_quad(u,v):
    """Sign of rational u+v sqrt(2), exactly."""
    if not u: return (v>0)-(v<0)
    if not v: return (u>0)-(u<0)
    if u>0 and v>0:return 1
    if u<0 and v<0:return -1
    diff=u*u-2*v*v
    assert diff
    sign=(diff>0)-(diff<0)
    return sign if u>0 else -sign


def lambda_power(k):
    u,v=1,0
    base=(1,1) if k>=0 else (-1,1)
    for _ in range(abs(k)):
        u,v=u*base[0]+2*v*base[1],u*base[1]+v*base[0]
    assert sign_quad(u,v)>0
    return u,v


def atan_interval(a,terms):
    x=F(1,a)
    power=x
    total=F(0)
    for k in range(terms):
        total+=(-1)**k*power/(2*k+1)
        power*=x*x
    next_sum=total+(-1)**terms*power/(2*terms+1)
    return min(total,next_sum),max(total,next_sum)


def error_interval(target,r):
    lo,hi=target
    if r<lo:return lo-r,hi-r
    if r>hi:return r-hi,r-lo
    raise AssertionError('Increase target interval precision')


def e_interval(terms):
    total=term=F(1)
    for k in range(1,terms+1):
        term/=k
        total+=term
    tail=term/F(terms+1)*F(terms+2,terms+1)
    return total,total+tail


def log_interval(x,terms=100):
    assert x>1
    y=(x-1)/(x+1)
    power=y
    total=F(0)
    for k in range(terms):
        total+=2*power/(2*k+1)
        power*=y*y
    tail=2*power/((2*terms+1)*(1-y*y))
    return total,total+tail


def decimal_interval(lo,hi,digits=12):
    scale=10**digits
    lower=(lo*scale).numerator//(lo*scale).denominator
    upper=-((-hi*scale).numerator//(-hi*scale).denominator)
    def write(k):
        sign='-' if k<0 else ''
        s=str(abs(k)).zfill(digits+1)
        return sign+s[:-digits]+'.'+s[-digits:]
    return {'lower':write(lower),'upper':write(upper),
            'endpoints':'exact decimal rationals, outward enclosure'}


def constants():
    scale=10**60
    root=isqrt(2*scale*scale)
    sqrt_lo,sqrt_hi=F(root,scale),F(root+1,scale)
    assert sqrt_lo**2<2<sqrt_hi**2
    a_lo=log_interval(1+sqrt_lo)[0]
    a_hi=log_interval(1+sqrt_hi)[1]
    log2=log_interval(F(2))
    b_lo,b_hi=log2[0]/2,log2[1]/2
    critical_lo=2*a_lo/(a_lo+b_hi)
    critical_hi=2*a_hi/(a_hi+b_lo)
    thresholds=[]
    for c in (F(1),F(3,2),F(2),F(5,2),F(3)):
        coef=1-2/c
        lo=coef*(a_lo if coef>=0 else a_hi)+b_lo
        hi=coef*(a_hi if coef>=0 else a_lo)+b_hi
        thresholds.append({'c':str(c),'theta':decimal_interval(lo,hi),
                           'strictly_positive_certified':lo>0})
    return {'critical_c_star':decimal_interval(critical_lo,critical_hi),
            'a_log_lambda':decimal_interval(a_lo,a_hi),
            'b_half_log_two':decimal_interval(b_lo,b_hi),
            'thresholds':thresholds}


def finite_candidate_audits():
    a5=atan_interval(5,500)
    a239=atan_interval(239,200)
    pi_bounds=(16*a5[0]-4*a239[1],16*a5[1]-4*a239[0])
    e_bounds=e_interval(80)
    assert e_bounds[1]<3
    A0,A=1,1
    B0,B=0,1
    factorial=1
    L=1
    counts={'exact_signed_prime_budget':0,'beta_and_saturation':0,
            'quadratic_field_uniform_ell_bound':0,
            'all_real_c_bounds_at_rational_samples':0,
            'deficit_lcm_divisibility':0,'prime_caps_and_count_bounds':0}
    exponents=(F(1),F(3,2),F(2),F(5,2),F(3))
    accepted={str(c):[] for c in exponents}
    rejected={str(c):[] for c in exponents}
    samples=[]
    ps=primes_to(100)
    for n in range(1,101):
        if n>1:
            A0,A=A,(2*n-1)*A+(n-1)**2*A0
            B0,B=B,(2*n-1)*B+(n-1)**2*B0
        factorial*=n
        L=lcm(L,n)
        g=gcd(A,4*B)
        ratio=F(g,factorial)
        E,D=odd(ratio.numerator),odd(ratio.denominator)
        d=valuation(g,2)-valuation(factorial,2)
        dA=valuation(A,2)-valuation(factorial,2)
        t=valuation(A,2)-valuation(g,2)
        beta=F(n,2)+dA
        assert ratio==power2(d)*E/D
        assert d==dA-t
        counts['exact_signed_prime_budget']+=1
        assert valuation(A,2)==a_valuation_formula(n)
        assert beta==beta_formula(n)
        assert t==(1 if n==3 else 0)
        if n>=4:assert (A//g)%2==1
        counts['beta_and_saturation']+=1
        assert odd(L)%D==0
        counts['deficit_lcm_divisibility']+=1
        ell=F(A,factorial)
        u,v=lambda_power(n)
        assert sign_quad((2*n+1)*ell-u,-v)>=0
        assert sign_quad(u-ell,v)>=0
        counts['quadratic_field_uniform_ell_bound']+=1
        r=F(4*B,A)
        q=r.denominator
        eps_lo,eps_hi=error_interval(pi_bounds,r)
        assert F(2,(2*n+1)*ell**2)<eps_lo
        assert eps_hi<F(4,(2*n+1)*ell**2)
        Z=F(factorial,g)
        for c in exponents:
            numerator,denominator=c.numerator,c.denominator
            k=(numerator-2*denominator)*n
            lu,lv=lambda_power(k)
            lower_coef=F(2,2*n+1)**denominator*Z**numerator
            lower_coef/=F(2*n+1)**max(numerator-2*denominator,0)
            upper_coef=F(4,2*n+1)**denominator*Z**numerator
            upper_coef*=F(2*n+1)**max(2*denominator-numerator,0)
            actual_lo=eps_lo**denominator*q**numerator
            actual_hi=eps_hi**denominator*q**numerator
            assert sign_quad(actual_lo-lower_coef*lu,-lower_coef*lv)>0
            assert sign_quad(upper_coef*lu-actual_hi,upper_coef*lv)>0
            counts['all_real_c_bounds_at_rational_samples']+=1
            if actual_hi<e_bounds[0]**denominator:
                accepted[str(c)].append(n)
            elif actual_lo>e_bounds[1]**denominator:
                rejected[str(c)].append(n)
            else:
                raise AssertionError(('Insufficient precision',n,c))
        factors=[]
        remainder=E
        for p in ps:
            if remainder%p:continue
            e=valuation(remainder,p)
            factors.append((p,e))
            remainder//=p**e
            assert n>=2
            cap=min(floor_log(n-1,p)-valuation(n,p),
                    floor_log(n+1,p)-valuation(n+1,p))
            assert e<=cap
            assert p**e<=n-1
            if p*p>n-1:assert e==1
        assert remainder==1
        if n>=2: assert E<=(n-1)**len(factors)
        counts['prime_caps_and_count_bounds']+=1
        if n in (10,20,50,72,100):
            samples.append({'n':n,'E_positive_odd':str(E),'D_negative_odd':str(D),
                            'd_binary':d,'beta':str(beta),'t':t,
                            'positive_odd_primes':[p for p,e in factors],
                            'reduced_q_digits':len(str(q))})
    return {'counts':counts,'target':'epsilon <= e/q^c',
            'strictly_certified_accepted_indices_through_100':accepted,
            'strictly_certified_rejected_count_through_100':{k:len(v) for k,v in rejected.items()},
            'selected_exact_signed_budgets':samples,
            'scope':'Finite candidate classification only; does not predict all larger indices.'}


def necessary_prime_count_50000():
    n=50000
    dA=a_valuation_formula(n)-(n-n.bit_count())
    exponent=1-2*dA
    assert dA==-24994 and exponent==49989
    rhs=2**exponent
    base=3*(2*n+1)
    def minimum(M):
        lo,hi=0,1
        while base*M*M*(n-1)**(2*hi)<=rhs:hi*=2
        while lo<hi:
            mid=(lo+hi)//2
            if base*M*M*(n-1)**(2*mid)>rhs:hi=mid
            else:lo=mid+1
        assert base*M*M*(n-1)**(2*lo)>rhs
        if lo:assert base*M*M*(n-1)**(2*(lo-1))<=rhs
        return lo
    M=1
    caps=[]
    for p in primes_to(isqrt(n-1)):
        cap=max(0,min(floor_log(n-1,p)-valuation(n,p),
                      floor_log(n+1,p)-valuation(n+1,p)))
        M*=p**cap
        caps.append({'p':p,'cap':cap})
    total=minimum(1)
    large=minimum(M)
    assert len(caps)==47 and total==1601 and large==1561
    return {'n':n,'target':'epsilon <= e/q^2','uses_e_less_than_3':True,
            'd_A':dA,'required_power_of_two':exponent,
            'necessary_positive_odd_prime_count_at_least':total,
            'small_prime_cutoff':isqrt(n-1),'small_odd_prime_count':len(caps),
            'small_prime_caps':caps,'maximum_small_prime_budget':str(M),
            'necessary_large_positive_odd_prime_count_at_least':large,
            'comparisons':'At each reported minimum N, the exact integer inequality passes; at N-1 it fails.',
            'scope':'Unconditional necessary counts only. No A_50000 or B_50000 recurrence or successful approximant is computed.'}


def main():
    result={'status':'passed','arithmetic':'exact rationals, integers and Q(sqrt(2)); Python standard library',
            'constants':constants(),
            'finite_candidates':finite_candidate_audits(),
            'integer_only_density_illustration':necessary_prime_count_50000(),
            'scope':'All-n budget and necessary-density theorems are proved analytically; no all-n subexponential excess bound or global pi irrationality exponent is claimed.'}
    Path(__file__).with_name('net_prime_budget_verification.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'constants':result['constants'],
                      'finite_counts':result['finite_candidates']['counts'],
                      'accepted_indices':result['finite_candidates']['strictly_certified_accepted_indices_through_100'],
                      'n_50000_necessary_counts':{
                          k:result['integer_only_density_illustration'][k]
                          for k in ('necessary_positive_odd_prime_count_at_least',
                                    'necessary_large_positive_odd_prime_count_at_least')}},indent=2))


if __name__=='__main__':
    main()
