"""Exact certificates for the mean-degree collision relaxation.

Uses only the Python standard library. No period-rank or pi assertion is
inferred from the relaxed multiplicity profiles or the finite tests.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, ceil, gcd
from pathlib import Path
import json


def continuum_value(m,t):
    mu=F(m,m+1)
    variance=F(m,(m+1)**2*(m+2))
    boundary=F(m,m+2)
    if t<=boundary:
        return F(2*(m+1),m+2)*(boundary/t)**m
    if t<=mu:
        return 1+(mu-t)**2/variance
    return F(1)


def integrate(m,polynomial,upper=F(1)):
    return sum(F(m)*coefficient*upper**(m+j)/F(m+j)
               for j,coefficient in enumerate(polynomial))


def square(polynomial):
    answer=[F(0)]*(2*len(polynomial)-1)
    for i,a in enumerate(polynomial):
        for j,b in enumerate(polynomial):
            answer[i+j]+=a*b
    return answer


def continuum_checks():
    checked=0
    for m in range(1,9):
        mu=F(m,m+1)
        variance=F(m,(m+1)**2*(m+2))
        boundary=F(m,m+2)
        for A in (F(1,4),F(1,2),F(3,4),F(1)):
            multiplier=F(m+1)/A**(m+1)
            polynomial=[multiplier*A,-multiplier]
            t=boundary*A
            assert integrate(m,polynomial,A)==1
            assert integrate(m,[F(0)]+polynomial,A)==t
            assert integrate(m,square(polynomial),A)==continuum_value(m,t)
            checked+=1
        for z in (F(0),F(1,4),F(1,2),F(3,4),F(1)):
            t=boundary+z*(mu-boundary)
            slope=(mu-t)/variance
            polynomial=[1+slope*mu,-slope]
            assert min(polynomial[0],sum(polynomial))>=0
            assert integrate(m,polynomial)==1
            assert integrate(m,[F(0)]+polynomial)==t
            assert integrate(m,square(polynomial))==continuum_value(m,t)
            checked+=1
        value=F(2*(m+1),m+2)
        assert continuum_value(m,boundary)==value
        assert -m*value/boundary==-2*(m+1)
        assert m*(m+1)*value/boundary**2==2/variance
        assert continuum_value(m,mu)==1
    return checked


def types(weights,H):
    return [a for a in product(*(range(H//w+1) for w in weights))
            if sum(w*x for w,x in zip(weights,a))<=H]


def real_optimizer(costs,M,t):
    """Exact finite optimum with sum n=M and sum cost*n<=M*t.

    Return a nonnegative affine-cutoff profile and certify the KKT
    conditions exactly; these conditions suffice by strict convexity.
    """
    xs=sorted(costs)
    count=len(xs)
    if t>=sum(xs)/count:
        profile=[F(M,count)]*count
        return profile,F(M*M,count),sum(xs)/count
    if t==xs[0]:
        zeros=sum(x==xs[0] for x in xs)
        profile=[F(M,zeros) if x==xs[0] else F(0) for x in xs]
        return profile,F(M*M,zeros),xs[0]
    first=second=F(0)
    for number,x in enumerate(xs,1):
        first+=x
        second+=x*x
        denominator=number*second-first*first
        if not denominator:
            continue
        B=F(M)*(first-number*t)/denominator
        A=(F(M)+B*first)/number
        if B<=0 or A-B*x<0:
            continue
        if number<count and A-B*xs[number]>0:
            continue
        profile=[max(F(0),A-B*y) for y in xs]
        assert sum(profile)==M
        mean=sum(y*n for y,n in zip(xs,profile))/M
        assert mean==t
        value=sum(n*n for n in profile)
        assert value==M*(A-B*t)
        assert all((n==A-B*y and n>=0) if n else A-B*y<=0
                   for y,n in zip(xs,profile))
        return profile,value,mean
    raise AssertionError('No finite optimizer found')


def compositions(total,length):
    if length==1:
        yield (total,)
        return
    for first in range(total+1):
        for tail in compositions(total-first,length-1):
            yield (first,)+tail


def exhaustive_profiles():
    checked=cases=0
    for weights,H in (((1,),3),((1,1),2),((1,2),3)):
        costs=sorted(F(sum(w*a for w,a in zip(weights,index)),H)
                     for index in types(weights,H))
        for M in (2,3,4):
            for t in (F(1,4),F(1,2),F(3,4)):
                profile,relaxed,mean=real_optimizer(costs,M,t)
                observed=None
                for counts in compositions(M,len(costs)):
                    weighted=sum(x*n for x,n in zip(costs,counts))
                    if weighted>M*t:
                        continue
                    value=sum(n*n for n in counts)
                    assert value>=relaxed
                    # Projection inequality follows from the exact dual signs.
                    distance=sum((F(n)-r)**2 for n,r in zip(counts,profile))
                    assert value-relaxed>=distance
                    observed=value if observed is None else min(observed,value)
                    checked+=1
                assert observed is not None
                assert mean<=t
                cases+=1
    return {'finite_problems':cases,'feasible_integer_profiles':checked}


def packet_checks():
    records=[]
    for weights in ((2,),(2,3),(2,3,5)):
        m=len(weights)
        K=2;theta=F(1,2);w0=1
        v0=2*K*theta**m*w0
        source_weights=tuple(2*w for w in weights)
        for H in (24,48,96):
            all_types=types(weights,H)
            costs=[F(sum(w*a for w,a in zip(weights,index)),H)
                   for index in all_types]
            M=0
            for beta in types(source_weights,H-1):
                gap=H-sum(v*b for v,b in zip(source_weights,beta))
                M+=K*ceil(F(gap)/v0)
            raw_count=sum((H-sum(w*a for w,a in zip(weights,index)))//w0+1
                          for index in all_types)
            assert 0<M<=raw_count
            for t in (F(1,3),F(m,m+2),F(m,m+1),F(9,10)):
                _,value,mean=real_optimizer(costs,M,t)
                continuum=F(1,2*(m+1)*w0)*continuum_value(m,t)
                records.append({'m':m,'weights':weights,'H':H,'types':len(costs),
                                'source_rows':M,'raw_columns':raw_count,
                                'mean_upper_bound':str(t),'actual_relaxed_mean':str(mean),
                                'finite_optimal_collision_over_MH':str(value/(M*H)),
                                'proved_continuum_limit':str(continuum)})
    return records


def concentration_checks():
    records=[]
    theta=F(1,2);B=F(3,5);C=F(13,10)
    for m,t in ((40,F(497,1000)),(80,F(541,1000)),
                (120,F(558,1000)),(200,F(573,1000))):
        lower=F(1,3*(m+2))*(B*F(m,m+2)/t)**m
        assert lower>2
        K=(C**m).numerator//(C**m).denominator
        v0=2*K*(theta/B)**m
        # ln2<7/10 and the previously proved arithmetic upper bound.
        # Remaining transverse-weight errors can be made arbitrarily small.
        base_budget=14*K*B**m+1-F(m,m+2)*theta+F(28,5)/v0
        assert base_budget<1
        source_mean=F(m,m+2)*theta
        assert source_mean<t<F(m,m+2)
        records.append({'m':m,'target_mean_threshold':str(t),
                        'certified_collision_saving_lower':str(lower),
                        'arithmetic_plus_analytic_budget_upper_before_weight_errors':str(base_budget),
                        'forced_normalized_displacement_lower':str(t-source_mean)})
    return records


def lcm_checks():
    lcm=1
    count=0
    for N in range(1,201):
        lcm=lcm*N//gcd(lcm,N)
        r=N//2
        central=comb(2*r,r)
        assert lcm%central==0
        assert central*(2*r+1)>=2**(2*r)
        if N>=2:
            assert central*(N+1)>=2**(N-1)
            count+=1
    return count


def parameter_checks():
    checked=0
    for theta in (F(1,4),F(1,2),F(3,4),F(9,10)):
        for b in range(1,100):
            B=F(b,100)
            for c in range(101,201):
                C=F(c,100)
                if C*B<=1 and C*theta>B:
                    assert B*B<theta
                    assert 1+theta-2*B>0
                    checked+=1
    B=F(3,5);C=F(13,10);theta=F(1,2)
    assert C*B==F(39,50)
    assert C*theta/B==F(13,12)
    assert 1+theta-2*B==F(3,10)
    assert F(7071,5000)**2<2<F(14143,10000)**2
    return {'admissible_rational_parameter_triples':checked,
            'original_family_relaxed_gap':'3/10',
            'optimized_half_theta_gap_interval':['857/10000','429/5000'],
            'interval_meaning':'strictly between; exact value 3/2-sqrt(2)'}


if __name__=='__main__':
    result={'status':'passed','continuum_optimizer_integrals':continuum_checks(),
            'finite_profile_checks':exhaustive_profiles(),
            'fixed_weight_packet_checks':packet_checks(),
            'concentration_thresholds':concentration_checks(),
            'lcm_divisibility_checks':lcm_checks(),
            'parameter_checks':parameter_checks(),
            'scope':'Exact mass/mean relaxation and parameter arithmetic. '
                    'No rank certificate for the optimizing targets and no pi claim. '
                    'Asymptotic statements follow from the written proofs.'}
    Path(__file__).with_name('moment_collision_verification.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:len(v) if isinstance(v,list) else v for k,v in result.items()},indent=2))
