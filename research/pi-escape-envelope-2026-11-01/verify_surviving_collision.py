"""Exact checks of surviving-type minimization and collision constants.

Uses standard-library integer/rational arithmetic. Asymptotics and the actual
period inequalities are proved in the accompanying note, not by these samples.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, ceil
from pathlib import Path
import json


def integer_minimum(caps, retained):
    assert 0 <= retained <= sum(caps)
    if retained == sum(caps):
        return sum(x*x for x in caps)
    level = 0
    while sum(min(n, level+1) for n in caps) <= retained:
        level += 1
    base = [min(n, level) for n in caps]
    residual = retained-sum(base)
    assert residual < sum(n > level for n in caps)
    return sum(x*x for x in base)+residual*(2*level+1)


def continuous_minimum(caps, retained):
    if retained == sum(caps):
        return F(sum(x*x for x in caps)), F(max(caps))
    remaining = F(retained)
    fixed_cost = F(0)
    ordered = sorted(caps)
    for i, cap in enumerate(ordered):
        level = remaining/(len(ordered)-i)
        if level <= cap:
            return fixed_cost+(len(ordered)-i)*level*level, level
        remaining -= cap
        fixed_cost += cap*cap
    raise AssertionError('No capped minimizer found')


def exhaustive_checks():
    cases = assignments = 0
    for caps in ((1,), (2, 4), (1, 3, 4), (2, 2, 5), (1, 2, 3, 4)):
        best = {}
        for row in product(*(range(c+1) for c in caps)):
            count = sum(row)
            cost = sum(x*x for x in row)
            best[count] = min(best.get(count, cost), cost)
            assignments += 1
        for retained, cost in best.items():
            real_cost, level = continuous_minimum(caps, retained)
            assert integer_minimum(caps, retained) == cost
            assert real_cost <= cost <= real_cost+2*sum(caps)+len(caps)
            assert sum(min(F(n), level) for n in caps) == retained
            cases += 1
    return cases, assignments


def packet_caps(H, v0, vi, K):
    caps = []
    for beta in product(*(range((H-1)//v+1) for v in vi)):
        gap = H-sum(v*b for v, b in zip(vi, beta))
        if gap > 0:
            caps.append(K*ceil(F(gap, v0)))
    direct = 0
    for s in range(ceil(F(H, v0))):
        for beta in product(*(range((H-v0*s-1)//v+1) for v in vi)):
            if v0*s+sum(v*b for v, b in zip(vi, beta)) < H:
                direct += K
    assert sum(caps) == direct
    return caps


def packet_checks():
    records = []
    for vi in ((3,), (3, 7), (3, 7, 15)):
        m = len(vi)
        for H in (30, 60, 120):
            caps = packet_caps(H, 1, vi, 2)
            M = sum(caps)
            for t in (F(1, 4), F(1, 2), F(3, 4)):
                rho = t**(m+1)
                escaped = int(rho*M)
                cost = integer_minimum(caps, M-escaped)
                real_cost, level = continuous_minimum(caps, M-escaped)
                profile = 1-(m+2)*rho+(m+1)*t**(m+2)
                predicted = F(4, m+2)*profile
                assert real_cost <= cost
                assert profile >= F(m+2, 2*(m+1))*(1-rho)**2
                records.append({'m':m, 'vi':vi, 'H':H, 'K':2, 'rows':M,
                                'escaped':escaped, 'target_fraction':str(rho),
                                'real_cap_level':str(level),
                                'exact_minimum_survivor_collision':cost,
                                'normalized_finite_minimum':str(F(cost, M*H)),
                                'proved_limit_at_target_fraction':str(predicted)})
    return records


def ordinary_coefficient_checks():
    count = 0
    # Rational surrogate omega=6i, testing coefficient estimates only.
    # No assertion that e^(6i)=1 is made.
    for K in (1, 2, 3, 5):
        for j in range(K):
            for d in range(25):
                for s in range(d+1):
                    assert comb(d, s)*(6*j)**(d-s) <= (7*K)**d
                    count += 1
    return count


def dimension_checks():
    records = []
    w_min = 10**100
    # The positive atanh series gives 2/3<ln2<25/36.
    for m, epsilon_bound in ((40, F(58,100)), (80,F(22,1000)),
                             (120,F(7,10000)), (160,F(21,1000000))):
        K = (13**m)//(10**m)
        w0 = F(5,3)**m
        v0 = 2*K*F(1,2)**m*w0
        epsilon_sq = 12*(m+1)*F(5,6)**m
        assert epsilon_sq < epsilon_bound**2 < 1
        cost_upper = (1-F(m,2*(m+2))+F(50,9)/v0+14*K/w0
                      + F(28*K+2,w_min)+F(m,(m+2)*w_min))
        assert cost_upper < 1
        records.append({'m':m,'K':K,'v0':str(v0),
                        'retained_fraction_epsilon_squared':str(epsilon_sq),
                        'certified_epsilon_upper_bound':str(epsilon_bound),
                        'conditional_analytic_plus_arithmetic_cost_upper':str(cost_upper),
                        'minimum_transverse_weight_assumed':str(w_min),
                        'scope':'Conditional parameters only; arbitrary interpolation '
                                'thresholds must also hold. No approximation tuple is constructed.'})
    for m in (1,2,3,10,40):
        for numerator in range(21):
            t=F(numerator,20)
            rho=t**(m+1)
            profile=1-(m+2)*rho+(m+1)*t**(m+2)
            assert profile>=F(m+2,2*(m+1))*(1-rho)**2
    return records


if __name__ == '__main__':
    cases, assignments = exhaustive_checks()
    records=packet_checks()
    result={'status':'passed','exhaustive_survivor_allocations':assignments,
            'exact_optimization_cases':cases,'complete_packet_cases':records,
            'ordinary_coefficient_checks':ordinary_coefficient_checks(),
            'conditional_highdimensional_bounds':dimension_checks(),
            'scope':'Finite rational minimization, source counts, coefficient '
                    'bounds and parameter certificates. Infinite conclusions '
                    'use the written proofs. No pi approximation or endpoint '
                    'contradiction is certified.'}
    Path(__file__).with_name('surviving_collision_verification.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:(len(v) if isinstance(v,list) else v)
                      for k,v in result.items()},indent=2))
