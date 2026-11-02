"""Independent exact checks for the first-layer global-budget audit.

No author verification module is imported. The separate finite scan uses the
integral U recurrence modulo each prime, without the reflection shortcut.
All-index statements are supplied by the written audit, not extrapolation.
"""
from fractions import Fraction as F
from math import factorial, gcd, isqrt, prod
from hashlib import sha256
from pathlib import Path
import json


def primes_up_to(n):
    out=[]
    for x in range(2,n+1):
        if all(x%p for p in out if p*p<=x):
            out.append(x)
    return out


def v2(x):
    assert x>0
    return (x&-x).bit_length()-1


def odd(x):
    return x>>v2(x)


def prime_v(x,p):
    out=0
    while x%p==0:
        x//=p
        out+=1
    return out


def quadratic_power(exponent):
    # Exact pair (a,b) representing (2+sqrt(2))**exponent.
    a,b=1,0
    c,d=2,1
    while exponent:
        if exponent&1:
            a,b=a*c+2*b*d,a*d+b*c
        c,d=c*c+2*d*d,2*c*d
        exponent//=2
    return a,b


def below_quadratic(integer,pair):
    a,b=pair
    return integer<=a or (integer-a)**2<=2*b*b


def finite_budget(limit=500):
    A,B=[1,1],[0,1]
    for n in range(1,limit):
        A.append((2*n+1)*A[-1]+n*n*A[-2])
        B.append((2*n+1)*B[-1]+n*n*B[-2])
    fac=[factorial(n) for n in range(limit+1)]
    U=[A[n]*2**n//fac[n] for n in range(limit+1)]
    assert all(A[n]*2**n%fac[n]==0 for n in range(limit+1))
    V=[odd(u) for u in U]
    prefix=[1]
    for n in range(limit+1):
        assert 2*v2(U[n])>=n
        assert below_quadratic(V[n],quadratic_power(n))
        prefix.append(prefix[-1]*V[n])
    ps=[p for p in primes_up_to(limit) if p%2]
    L=1
    count=fold_count=edge_count=local_count=0
    for n in range(2,limit+1):
        # At entry L=lcm(1,...,n-1).
        L=L*(n-1)//gcd(L,n-1)
        Lo=odd(L)
        g=gcd(A[n],4*B[n])
        R=F(g,fac[n])
        E,D=odd(R.numerator),odd(R.denominator)
        d=v2(R.numerator)-v2(R.denominator)
        upper=[p for p in ps if n<2*p and p<n]
        exceptional=[p for p in upper if A[n-p]%p==0]
        lifted=[p for p in upper if g%(p*p)==0]
        H,P,T=prod(upper),prod(exceptional),prod(lifted)
        left=E*odd(n)*H*H
        assert D*P*P*Lo%left==0
        assert D*P*T*Lo%left==0
        for p in upper:
            actual=prime_v(g,p)-1
            assert actual==-1+(p in exceptional)+(p in lifted)
            local_count+=1
        previous_lower=F(2,2*n+1)/R**2
        new_lower=F(2)**(1-2*d)*odd(n)**2*H**4/F((2*n+1)*P**4*Lo**2)
        sharper_lower=F(2)**(1-2*d)*odd(n)**2*H**4/F((2*n+1)*P**2*T**2*Lo**2)
        assert previous_lower>=sharper_lower>=new_lower
        count+=1
        folds=[]
        for p in exceptional:
            k=min(n-p,2*p-n-1)
            assert 0<=k<=(n-1)//3
            assert V[k]%p==0
            folds.append((p,k))
            fold_count+=1
        # Several genuinely different cutoffs; no asymptotic inference.
        cube=0
        while (cube+1)**3<=n:
            cube+=1
        cuts={0,1,2,isqrt(n),cube,(n-1)//3}
        for K in cuts:
            K=min(K,(n-1)//3)
            Pedge=prod(p for p,k in folds if k<=K)
            capacity=prefix[K+1]
            assert capacity%Pedge==0
            assert below_quadratic(capacity,quadratic_power(K*(K+1)//2))
            edge_count+=1
    reflection=0
    for p in [p for p in ps if p<=199]:
        ell=[A[r]*pow(fac[r],-1,p)%p for r in range(p)]
        for r in range(p):
            assert ell[p-1-r]==(-1)**((p-1)//2+r)*ell[r]%p
            reflection+=1
    return {'maximum_index':limit,'exact_aggregate_divisibilities_each':count,
            'exact_finite_quality_bounds':count,
            'top_half_valuation_checks':local_count,
            'reflection_values_modulo_primes_through_199':reflection,
            'folded_anomaly_checks':fold_count,
            'odd_height_bounds':limit+1,'edge_capacity_checks':edge_count}


def independent_scan(limit=10000):
    factors=[[] for _ in range(limit+1)]
    steps=0
    for p in [p for p in primes_up_to(limit-1) if p%2]:
        end=min(p-1,limit-p)
        inverses=[0]*(end+1)
        if end>=1:
            inverses[1]=1
        before,current=1,2
        for r in range(2,end+1):
            inverses[r]=-(p//r)*inverses[p%r]%p
            before,current=current,(2*(2*r-1)*current+4*(r-1)*before)*inverses[r]%p
            if current==0:
                factors[p+r].append(p)
            steps+=1
    events=0
    nonempty=0
    largest=0
    digest=sha256()
    for n in range(1,limit+1):
        ps=factors[n]
        P=prod(ps)
        # An integer proof of P**6<e**n, with 8/3<e.
        value=P**6
        if value.bit_length()>n:
            assert value*3**n<8**n
        # Strict maximum at n=21; comparison uses only integers.
        if n==21:
            assert P==17
        else:
            assert P**21<17**n
        if ps:
            nonempty+=1
        largest=max(largest,len(ps))
        events+=len(ps)
        for p in ps:
            digest.update(f'{n}:{p}\n'.encode())
    assert events==957 and nonempty==911 and largest==3
    return {'limit':limit,'method':'Unreflected integral-U modular recurrence',
            'modular_steps':steps,'anomaly_events':events,
            'nonempty_indices':nonempty,'maximum_cluster_size':largest,
            'S_n_less_than_n_over_6_each_index':limit,
            'unique_maximum_ratio':'ln(17)/21 at n=21',
            'ascending_event_stream_sha256':digest.hexdigest()}


def triplet_lifts():
    # These six first-layer clusters are the author's finite observations.
    # Their gcd valuations are checked anew by direct recurrence mod p**2.
    rows=((2772,(2383,2671,2689)),(9559,(5741,6841,7459)),
          (29259,(21139,23027,25343)),(29158,(17359,23327,26497)),
          (40219,(28393,29723,37813)),(73219,(37813,44959,53231)))
    output=[]
    for n,ps in rows:
        checks=[]
        for p in ps:
            assert n<2*p and p<n
            assert all(p%d for d in range(2,isqrt(p)+1))
            m=p*p
            a0,a1,b0,b1=1,1,0,1
            for k in range(1,n):
                a0,a1=a1,((2*k+1)*a1+k*k*a0)%m
                b0,b1=b1,((2*k+1)*b1+k*k*b0)%m
            assert a1%p==0 and b1%p==0
            assert a1!=0 or b1!=0
            checks.append({'p':p,'v_p_g':1,'v_p_g_minus_v_p_factorial':0,
                           'A_n_mod_p_squared':a1,'B_n_mod_p_squared':b1})
        output.append({'n':n,'modular_checks':checks})
    return {'triples_checked':len(rows),'prime_pairs_checked':18,
            'method':'Direct A,B recurrence modulo p squared; no shifted-transfer verifier imported',
            'all_gcd_valuations_equal_one':True,
            'no_positive_excess_among_these_first_layer_clusters':True,
            'rows':output}


def main():
    data={'status':'passed','arithmetic':'Exact integers and rationals; no floating logarithms',
          'finite_aggregate_and_fold':finite_budget(),
          'independent_direct_scan':independent_scan(),
          'independent_triplet_lifts':triplet_lifts(),
          'scope':'Finite independent audits only. Written arguments prove the all-index budget, reflection and edge capacity. Eventual S(n)<=n/6 is not proved.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n')
    shown=dict(data)
    shown['independent_triplet_lifts']={k:v for k,v in data['independent_triplet_lifts'].items() if k!='rows'}
    print(json.dumps(shown,indent=2))


if __name__=='__main__':
    main()
