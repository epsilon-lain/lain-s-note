"""Exact checks of grouped transport coefficients and integer-grid bases.

Uses only Python's standard library. The period-function determinant is not
estimated here, and no Diophantine conclusion is inferred from finite tests.
"""
from fractions import Fraction as F
from itertools import product, combinations
from math import comb, factorial, prod
from pathlib import Path
import json


def determinant(matrix):
    a=[list(map(F,row)) for row in matrix]
    ans=F(1)
    for k in range(len(a)):
        pivot=next((i for i in range(k,len(a)) if a[i][k]),None)
        if pivot is None:return F(0)
        if pivot!=k:
            a[k],a[pivot]=a[pivot],a[k]
            ans=-ans
        p=a[k][k]
        ans*=p
        for i in range(k+1,len(a)):
            q=a[i][k]/p
            for j in range(k+1,len(a)):
                a[i][j]-=q*a[k][j]
    return ans


def dot(a,b):return sum(x*y for x,y in zip(a,b))
def multi_binom(a,b):return prod(comb(x,y) if x>=y else 0 for x,y in zip(a,b))
def mpower(x,a):return prod(v**e for v,e in zip(x,a))
def lower_indices(m,d):
    return sorted((a for a in product(range(d+1),repeat=m) if sum(a)<=d),key=lambda a:(sum(a),a))


def test_one_dimensional():
    cases=0;sharp_cases=0;partial_cases=0
    for d in range(9):
        for b in range(1,d+2):
            source=tuple(range(b))
            for target in combinations(range(d+1),b):
                got=determinant([[comb(a,beta) if a>=beta else 0 for a in target] for beta in source])
                expected=F(prod(target[j]-target[i] for i in range(b) for j in range(i+1,b)),prod(factorial(k) for k in range(b)))
                assert got==expected and got>0
                shift=sum(target)-b*(b-1)//2
                k=sum(a>=b for a in target)
                assert shift>=k*k
                partial_cases+=1
                minimum=tuple(range(b-k))+tuple(range(b,b+k))
                assert (shift==k*k)==(target==minimum)
                if k==b:assert got==1 if shift==b*b else True
                if target==minimum:sharp_cases+=1
                for zeta in (F(1,7),F(-2,5)):
                    transport=[[F(comb(a,beta))*((-zeta)**(a-beta)) if a>=beta else F(0) for a in target] for beta in source]
                    assert determinant(transport)==expected*(-zeta)**shift
                cases+=1
    return {'Pascal_minor_cases':cases,'quadratic_partial_boundary_cases':partial_cases,'sharp_minimum_cases':sharp_cases}


def test_rounded_grids():
    counts=0;integer_affine=0;records=[]
    for costs in ((1,1),(1,2),(2,3),(1,2,3)):
        m=len(costs)
        for cut in range(1,6):
            b=sorted((a for a in product(range(cut),repeat=m) if dot(costs,a)<cut),key=lambda a:(sum(a),a))
            for tau in (F(1),F(5,4),F(3,2),F(2)):
                nodes=lambda i,n:(tau.numerator*n)//tau.denominator+i
                target=[tuple(nodes(i,a[i]) for i in range(m)) for a in b]
                got=determinant([[multi_binom(a,beta) for a in target] for beta in b])
                expected=F(1)
                for beta in b:
                    for i in range(m):
                        expected*=F(prod(nodes(i,beta[i])-nodes(i,k) for k in range(beta[i])),factorial(beta[i]))
                assert got==expected and got>0 and got.denominator==1
                gamma=tuple(sum(a[i]-beta[i] for a,beta in zip(target,b)) for i in range(m))
                assert all(g>=0 for g in gamma)
                zeta=tuple(F(i+1,11) for i in range(m))
                transport=[[multi_binom(a,beta)*mpower(tuple(-z for z in zeta),tuple(x-y for x,y in zip(a,beta))) if all(x>=y for x,y in zip(a,beta)) else F(0) for a in target] for beta in b]
                assert determinant(transport)==got*mpower(tuple(-z for z in zeta),gamma)
                if tau.denominator==1:
                    assert got==tau**sum(map(sum,b))
                    integer_affine+=1
                counts+=1
                records.append({'costs':list(costs),'cut':cut,'size':len(b),'dilation':str(tau),'coefficient':str(got),'total_error_exponent':list(gamma)})
    return {'rounded_grid_cases':counts,'integer_affine_cases':integer_affine,'records':records}


def test_grouping():
    # Full packet with K=2,v0=1,v=(3,4),H=5, in two coordinates.
    v0=1;v=(3,4);height=5;m=2;centers=2
    zeta=(F(2,7),F(-1,5))
    groups=[]
    sources=[]
    for j in range(centers):
        for d in range(height):
            b=[beta for beta in lower_indices(m,d) if v0*(d-sum(beta))+dot(v,beta)<height]
            if not b:continue
            groups.append((j,d,b))
            sources += [(j,d,beta) for beta in b]
    def entry(source,target):
        j,d,beta=source;k,e,a=target
        if (j,d)!=(k,e) or any(x<y for x,y in zip(a,beta)):return F(0)
        return multi_binom(a,beta)*mpower(tuple(-z for z in zeta),tuple(x-y for x,y in zip(a,beta)))
    passed=0
    for seed in range(24):
        targets=[];expected=F(1)
        for g,(j,d,b) in enumerate(groups):
            options=list(combinations(lower_indices(m,d),len(b)))
            chosen=options[(seed*(g+1)+g)%len(options)]
            targets += [(j,d,a) for a in chosen]
            minor=determinant([[multi_binom(a,beta) for a in chosen] for beta in b])
            gamma=tuple(sum(a[i] for a in chosen)-sum(beta[i] for beta in b) for i in range(m))
            if minor:
                assert all(x>=0 for x in gamma)
                expected*=minor*mpower(tuple(-z for z in zeta),gamma)
            else:expected=F(0)
        actual=determinant([[entry(s,t) for t in targets] for s in sources])
        assert actual==expected
        passed+=1
    # Replace the target in the first singleton group by an unused row in
    # the last group. The cardinality mismatch forces a zero determinant.
    targets=list(sources)
    candidates=[(1,height-1,a) for a in lower_indices(m,height-1) if (1,height-1,a) not in targets]
    targets[0]=candidates[0]
    assert len(set(targets))==len(targets)
    assert determinant([[entry(s,t) for t in targets] for s in sources])==0
    return {'grouped_packet_cases':passed,'cardinality_mismatch_cases':1,'source_rows':len(sources)}


def test_cancellation_and_shell():
    b=((0,0),(1,0),(0,1));a=((0,0),(1,1),(2,2))
    coefficient=determinant([[multi_binom(target,beta) for target in a] for beta in b])
    assert coefficient==0
    z1,z2=F(1,7),F(-2,9)
    boundary=[[-z2,-2*z1*z2*z2],[-z1,-2*z1*z1*z2]]
    assert all(x for row in boundary for x in row)
    assert determinant(boundary)==0
    shell=(2,3)
    assert len([x for x in shell if x>=1])==2
    assert len([x for x in shell if x-2!=0])==1
    return {'nonzero_assignments_cancel':True,'shell_shadow_counterexample':True}


if __name__=='__main__':
    result={'status':'passed','one_dimensional':test_one_dimensional(),'rounded_grids':test_rounded_grids(),'grouping':test_grouping(),'counterexamples':test_cancellation_and_shell(),'scope':'Exact transport coefficients, quadratic error-order laws and lattice interpolation determinants only. No nonvanishing or size bound for the actual period mixed minors is established.'}
    Path(__file__).with_name('exterior_transport_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    concise=dict(result);concise['rounded_grids']={k:v for k,v in result['rounded_grids'].items() if k!='records'}
    print(json.dumps(concise,indent=2))
