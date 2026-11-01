"""Exact finite checks accompanying the positive-fraction transport proofs.

All arithmetic is rational or integral. The asymptotic statements are proved
in the note, not inferred from these finite samples.
"""
from fractions import Fraction as F
from itertools import product, combinations
from math import comb, factorial, prod
from pathlib import Path
import json


def lattice(costs,height):
    return [a for a in product(*(range((height-1)//c+1) for c in costs))
            if sum(c*x for c,x in zip(costs,a))<height]


def determinant(a):
    a=[[F(x) for x in row] for row in a];out=F(1)
    for k in range(len(a)):
        pivot=next((i for i in range(k,len(a)) if a[i][k]),None)
        if pivot is None:return F(0)
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];out=-out
        value=a[k][k];out*=value
        for i in range(k+1,len(a)):
            ratio=a[i][k]/value
            for j in range(k+1,len(a)):
                a[i][j]-=ratio*a[k][j]
    return out


def check_cost_minimum():
    checked=0
    for costs,height in [((1,),4),((2,),7),((1,2),3),((2,3),5),((1,2,3),3)]:
        inside=lattice(costs,height);b=len(inside)
        key=lambda a:(sum(c*x for c,x in zip(costs,a)),a)
        inside=sorted(inside,key=key)
        bigger=sorted(lattice(costs,4*height+sum(costs)),key=key)
        outside=[a for a in bigger if a not in set(inside)]
        pool=inside+outside[:min(3,b)]
        cost=lambda a:sum(c*x for c,x in zip(costs,a))
        initial=sum(map(cost,inside))
        minimum={}
        for target in combinations(pool,b):
            k=len(set(target)-set(inside))
            increase=sum(map(cost,target))-initial
            minimum[k]=min(minimum.get(k,increase),increase)
            checked+=1
        for k,value in minimum.items():
            expected=sum(cost(a) for a in outside[:k])-sum(cost(a) for a in inside[b-k:])
            assert value==expected
    return checked


def check_fiber_determinants():
    checked=0
    for costs,height in [((1,2),5),((2,3),7),((3,9),14),((2,3,5),8)]:
        source=lattice(costs,height)
        source.sort(key=lambda a:(a[1:],a[0]))
        ys=sorted({a[1:] for a in source})
        lengths={y:sum(a[1:]==y for a in source) for y in ys}
        for rho in (F(1,3),F(1,2),F(1)):
            target=[];expected=F(1);cost_total=0;replaced=0
            for y in ys:
                n=lengths[y];k=(rho*n).__floor__()
                js=list(range(n-k))+list(range(n,n+k))
                target.extend((x,)+y for x in js)
                numerator=prod(js[j]-js[i] for i in range(n) for j in range(i+1,n))
                expected*=F(numerator,prod(factorial(t) for t in range(n)))
                cost_total+=costs[0]*k*k;replaced+=k
            matrix=[[prod(comb(a,b) for a,b in zip(alpha,beta)) if all(a>=b for a,b in zip(alpha,beta)) else 0
                     for alpha in target] for beta in source]
            assert determinant(matrix)==expected>0
            assert len(set(target)-set(source))==replaced
            difference=sum(sum(c*x for c,x in zip(costs,a)) for a in target)-sum(sum(c*x for c,x in zip(costs,b)) for b in source)
            assert difference==cost_total
            checked+=1
    return checked


def packet_counts(costs,height,rho,K=2):
    v0=1;vi=tuple(c+1 for c in costs);theta=F(1,2)
    cmin=min(costs);c1=costs[0]
    cutoff=((1+rho)*height+rho*c1)/(cmin+(1+rho)*v0)
    total=moved=cost_total=0;fibers=0
    for N in range(height):
        h=height-N
        ys=lattice(costs[1:],h) if len(costs)>1 else [()]
        bulk=N>=cutoff
        for y in ys:
            if sum(y)>N:continue
            remaining=h-sum(c*x for c,x in zip(costs[1:],y))
            n=min((remaining+c1-1)//c1,N-sum(y)+1)
            if n<=0:continue
            total+=n
            if not bulk:continue
            assert n==(remaining+c1-1)//c1
            k=(rho*n).__floor__()
            if k:
                max_first=n+k-1
                a=(max_first,)+y
                assert sum(a)<=N
                assert sum(F(v,2)*x for v,x in zip(vi,a))<=height
                assert c1*n+sum(c*x for c,x in zip(costs[1:],y))>=h
            moved+=k;cost_total+=c1*k*k;fibers+=1
    # Independent enumeration in the original (s,beta) coordinates.
    original=sum(height-sum(v*x for v,x in zip(vi,beta))
                 for beta in lattice(vi,height))
    assert total==original
    m=len(costs);A=prod(F(v,c) for v,c in zip(vi,costs))
    t0=F(cmin,1)/(cmin+(1+rho))
    return {'costs':list(costs),'H':height,'rho':str(rho),'rows':K*total,
            'replaced_rows':K*moved,'weighted_displacement':K*cost_total,
            'actual_fraction':str(F(moved,total)),
            'actual_normalized_cost':str(F(cost_total,total*height)),
            'predicted_fraction_limit':str(rho*A*t0**(m+1)),
            'predicted_cost_limit':str(2*rho*rho*A*t0**(m+2)/F(m+2)),
            'bound_at_actual_fraction':str(F(moved,total)**2/(F(m+1)*A)),
            'checked_bulk_fibers':fibers}


def run():
    records=[packet_counts(costs,h,rho) for costs in ((3,9),(3,9,25))
             for h in (40,80,160) for rho in (F(1,2),F(1))]
    # Exact algebraic check of the dimension-dependent bracket constants.
    dimensional=[]
    for m in (1,2,3,10,40):
        for rho in (F(1,2),F(1)):
            lower=rho*rho/F(m+1);upper=2*rho*rho/F(m+2)
            assert lower<=upper
            dimensional.append({'m':m,'rho':str(rho),'universal_fraction_lower':str(lower),
                                'constructive_large_weight_upper':str(upper)})
    result={'status':'passed','exhaustive_target_set_checks':check_cost_minimum(),
            'nonzero_fiber_determinant_checks':check_fiber_determinants(),
            'complete_packet_support_checks':records,'large_weight_constants':dimensional,
            'scope':'Exact finite combinatorial minima, fiber Pascal determinants, full source count and target-budget checks. The written arguments prove asymptotics; no period-minor estimate, interpolation threshold, or pi approximation result is certified.'}
    Path(__file__).with_name('positive_fraction_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','exhaustive_target_set_checks':result['exhaustive_target_set_checks'],
                      'nonzero_fiber_determinant_checks':result['nonzero_fiber_determinant_checks'],
                      'complete_packet_cases':len(records)},indent=2))


if __name__=='__main__':run()
