"""Independent exact audits of the multiple-band and window results.

Python standard library only. No author verifier is imported. The finite
checks do not prove a global bound for S(n) or for the smooth prime part.
"""
from fractions import Fraction as Q
from math import gcd, isqrt
from pathlib import Path
import json


def odd(x):
    x=abs(x)
    assert x
    return x//(x&-x)


def valuation(x,p):
    assert x
    out=0
    while x%p==0:
        x//=p; out+=1
    return out


def primes(limit):
    return [p for p in range(3,limit+1,2)
            if all(p%d for d in range(2,isqrt(p)+1))]


def cube_floor(n):
    b=0
    while (b+1)**3<=n: b+=1
    return b


def floor_log(n,p):
    out=0
    while n>=p:
        n//=p; out+=1
    return out


def main():
    limit=500
    A=[1,1]; B=[0,1]; factorial=[1,1]
    for n in range(1,limit+6):
        A.append((2*n+1)*A[-1]+n*n*A[-2])
        B.append((2*n+1)*B[-1]+n*n*B[-2])
        factorial.append(factorial[-1]*(n+1))
    U=[]
    for n in range(len(A)):
        assert (2**n*A[n])%factorial[n]==0
        U.append(2**n*A[n]//factorial[n])
    ps=primes(limit)
    prefix=[1]
    for b in range(1,16): prefix.append(prefix[-1]*B[b])
    budgets=forced_deficits=quality_checks=high_gain_checks=0
    computable_floor_checks=0
    examples=[]
    for n in range(2,limit+1):
        g=gcd(A[n],4*B[n]); R=Q(g,factorial[n])
        E,D=odd(R.numerator),odd(R.denominator)
        C=odd(U[n]); d=valuation(R.numerator,2)-valuation(R.denominator,2)
        assert R==Q(2)**d*Q(E,D)
        # Independent prime-by-prime construction of the whole-prime sieve.
        # The author's implementation uses gcd saturation instead.
        Ghat=1
        for p in ps:
            alpha=valuation(C,p)
            cap=min(floor_log(n-1,p)-valuation(n,p),
                    floor_log(n+1,p)-valuation(n+1,p))
            if 0<alpha<=cap: Ghat*=p**alpha
        assert Ghat%E==0
        s=valuation(U[n],2)
        assert d<=s-n
        base_floor=Q(2,2*n+1)*Q(2)**(2*(n-s))/Ghat**2
        for cut in sorted({1,2,3,cube_floor(n),min(8,n)}):
            hp=[p for p in ps if p<=n and p*p>n and p*(cut+1)>n]
            H=1
            for p in hp: H*=p
            P=gcd(C,H); F=gcd(H,prefix[cut])
            M=H//gcd(H,P*F)
            assert D%M==0
            assert (D*P*F)%H==0
            budgets+=1
            eh=1
            for p in hp:
                b,r=divmod(n,p)
                assert 1<=b<=cut and b<p
                e=valuation(g,p)-valuation(factorial[n],p)
                if P%p and F%p:
                    assert e==-1
                    forced_deficits+=1
                if E%p==0:
                    assert E%(p*p)!=0 and P%p==0
                    eh*=p
            assert P%eh==0 and E%eh==0
            low_count=sum(p<=n and p not in hp for p in ps)
            assert E//eh <= (n-1)**low_count
            high_gain_checks+=1
            actual_floor=Q(2,2*n+1)/R**2
            proposed_floor=Q(2,2*n+1)*Q(2)**(-2*d)*Q(H,E*P*F)**2
            assert proposed_floor<=actual_floor
            quality_checks+=1
            new_floor=base_floor*M*M
            assert base_floor<=new_floor<=actual_floor
            computable_floor_checks+=1
            if n in (21,107,500) and cut==cube_floor(n):
                examples.append({'n':n,'bands':cut,'s':s,'C':str(C),
                                 'H':str(H),'P':str(P),'F':str(F),'M':str(M),
                                 'Ghat':str(Ghat),'E':str(E),'D':str(D),
                                 'gain_only_floor':str(base_floor),
                                 'improved_floor':str(new_floor),
                                 'exact_analytic_floor':str(actual_floor),
                                 'improvement_factor':str(M*M)})
    window_checks=0
    for n in range(1,61):
        for h in range(0,6):
            coeff={-1:(Q(0),Q(1)),0:(Q(1),Q(0))}
            den=1
            for j in range(h):
                coeff[j+1]=tuple((2*(2*n+2*j+1)*coeff[j][i]
                                    +4*(n+j)*coeff[j-1][i])/(n+j+1)
                                   for i in (0,1))
                den*=n+j+1
            cleared={j:tuple(den*t for t in v) for j,v in coeff.items()}
            assert all(t.denominator==1 for v in cleared.values() for t in v)
            assert all(v[0]*U[n]+v[1]*U[n-1]==U[n+j]
                       for j,v in coeff.items())
            for i in range(-1,h):
                j=i+1
                ai,bi=map(int,cleared[i]); aj,bj=map(int,cleared[j])
                witness=bj*U[n+i]-bi*U[n+j]
                assert den*witness==U[n]*(bj*ai-bi*aj)
                if not witness: continue
                assert witness%U[n]==0
                for degree,value in ((1,witness),(2,witness*witness)):
                    required=odd(U[n])//gcd(odd(U[n]),den**degree)
                    assert odd(value)%required==0
                    window_checks+=1
    midpoint_checks=0
    midpoint_zero=[]
    for p in ps:
        h=(p-1)//2
        character_sum=0
        for x in range(p):
            v=pow((x*(x*x+2*x+2))%p,h,p)
            character_sum+=(v if v in (0,1) else -1)
        trace=-character_sum
        assert (trace-U[h])%p==0
        assert trace*trace<=4*p
        if p>=5:
            assert (A[h]%p==0)==(trace==0)
            if trace==0:
                assert p%8==7
                midpoint_zero.append(p)
        midpoint_checks+=1
    edge_checks=0
    for n in range(2,301):
        for K in sorted({0,1,cube_floor(n),isqrt(n),n//4}):
            edge=[p for p in ps if 2*p>n and p<n and A[n-p]%p==0
                  and min(n-p,2*p-n-1)<=K]
            assert len(edge)<=2*(K+1)
            product=1
            for p in edge:
                k=min(n-p,2*p-n-1)
                assert p==n-k or 2*p==n+1+k
                assert p>=n-K or 2*p<=n+1+K
                product*=p
            assert product<=n**(2*(K+1))
            edge_checks+=1
    result={
        'status':'passed','arithmetic':'Independent exact integers and rational fractions',
        'maximum_budget_index':limit,'multiple_band_divisibility_checks':budgets,
        'forced_deficit_prime_checks':forced_deficits,
        'high_gain_and_low_capacity_checks':high_gain_checks,
        'finite_quality_floor_checks':quality_checks,
        'computable_overlap_certificate_checks':computable_floor_checks,
        'concrete_overlap_certificate_examples':examples,
        'formal_linear_and_quadratic_window_witness_checks':window_checks,
        'midpoint_point_count_checks':midpoint_checks,
        'midpoint_zero_primes_in_range':midpoint_zero,
        'edge_affine_count_checks':edge_checks,
        'scope':'Finite independent checks of stated identities. No uniform bound for S(n), the primorial gcd, or pi bad approximability.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
