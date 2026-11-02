"""Exact checks of selected-column costs and residual-degree identities.

Finite checks use the Python standard library. The weighted-simplex limits
are proved in the note, not inferred from these finite samples.
"""
from fractions import Fraction as F
from itertools import product, combinations
from math import comb, factorial, ceil, prod
from pathlib import Path
import json


def determinant(matrix):
    a=[[F(x) for x in row] for row in matrix]
    out=F(1)
    for k in range(len(a)):
        pivot=next((i for i in range(k,len(a)) if a[i][k]),None)
        if pivot is None:return F(0)
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];out=-out
        d=a[k][k];out*=d
        for i in range(k+1,len(a)):
            q=a[i][k]/d
            for j in range(k+1,len(a)):a[i][j]-=q*a[k][j]
    return out


def types(weights,H):
    return [a for a in product(*(range(H//w+1) for w in weights))
            if sum(w*x for w,x in zip(weights,a))<=H]


def columns(w0,weights,H):
    return [(h,a) for a in types(weights,H)
            for h in range((H-sum(w*x for w,x in zip(weights,a)))//w0+1)]


def cdf(m,a):return a**m*((m+1)-m*a)


def truncated_mean(m,a):
    return F(m,m+2)*a*((m+2)-(m+1)*a)/((m+1)-m*a)


def quantile_interval(m,density=F(1,2),steps=64):
    lo,hi=F(0),F(1)
    for _ in range(steps):
        mid=(lo+hi)/2
        if cdf(m,mid)<density:lo=mid
        else:hi=mid
    assert cdf(m,lo)<=density<=cdf(m,hi)
    return lo,hi


def exhaustive_column_cost():
    choices=counts=0
    for w0,weights,H in ((1,(2,),3),(2,(2,3),4),(1,(2,3),3)):
        cols=columns(w0,weights,H)
        costs=[sum(w*x for w,x in zip(weights,a)) for _,a in cols]
        for number in range(1,len(cols)+1):
            smallest=sum(sorted(costs)[:number])
            observed=None
            for sub in combinations(range(len(cols)),number):
                value=sum(costs[i] for i in sub)
                observed=value if observed is None else min(observed,value)
                choices+=1
            assert observed==smallest
            counts+=1
    return choices,counts


def finite_packet_counts():
    records=[]
    for weights in ((3,),(3,7),(3,7,15)):
        m=len(weights);theta=F(1,2);K=2;w0=1
        v0=2*K*theta**m*w0
        vi=tuple(2*w for w in weights)
        lo,hi=quantile_interval(m)
        for H in (40,80,160):
            cols=columns(w0,weights,H)
            M=0;source_cost=0
            for beta in types(vi,H-1):
                gap=H-sum(v*b for v,b in zip(vi,beta))
                count=K*ceil(F(gap)/v0)
                M+=count
                source_cost+=count*sum(w*b for w,b in zip(weights,beta))
            assert 0<M<=len(cols)
            costs=sorted(sum(w*a for w,a in zip(weights,alpha)) for h,alpha in cols)
            selected=sum(costs[:M])
            records.append({'m':m,'weights':weights,'H':H,'source_rows':M,
                            'all_columns':len(cols),'finite_density':str(F(M,len(cols))),
                            'minimum_column_transverse_mean':str(F(selected,M*H)),
                            'source_transverse_mean':str(F(source_cost,M*H)),
                            'finite_gap':str(F(selected-source_cost,M*H)),
                            'proved_limit_mean_interval':[str(truncated_mean(m,lo)),str(truncated_mean(m,hi))]})
    return records


def period_entry(row,column,Z=F(2)):
    j,s,a=row;h,alpha=column
    if any(b>x for b,x in zip(a,alpha)):return F(0)
    pref=prod(comb(x,b) for x,b in zip(alpha,a))
    degree=sum(alpha)-sum(a)
    return pref*sum(F(comb(degree,k)*h**(s-k),factorial(s-k))
                    *(j*Z)**(degree-k) for k in range(min(s,degree)+1))


def coefficient(a,d,column,R):
    h,alpha=column
    if any(b>x for b,x in zip(a,alpha)):return F(0)
    residual=sum(alpha)-sum(a)
    if d<residual:return F(0)
    return prod(comb(x,b) for x,b in zip(alpha,a))*R**d*F(h**(d-residual),factorial(d-residual))


def degree_checks():
    source=[(0,0),(1,0),(0,1)]
    candidates=[a for a in product(range(3),repeat=2) if sum(a)<=2]
    cols=[(0,(0,0)),(0,(1,0)),(0,(0,1)),(1,(1,0)),
          (1,(0,1)),(0,(2,0)),(0,(1,1)),(1,(2,1))]
    checked=nonzero=coefficient_checks=0
    R=F(5);q=(F(3),F(7));zeta=(F(1,10),F(1,20))
    for target in combinations(candidates,3):
        P=[[prod(comb(a_i,b_i) if a_i>=b_i else 0 for a_i,b_i in zip(a,b))
            for a in target] for b in source]
        pascal=determinant(P)
        Gamma=tuple(sum(a[i] for a in target)-sum(b[i] for b in source) for i in range(2))
        if pascal:assert min(Gamma)>=0
        rows=[(1,2-sum(a),a) for a in target]
        for chosen in combinations(cols,3):
            value=determinant([[period_entry(row,c) for c in chosen] for row in rows])
            D=tuple(sum(c[1][i] for c in chosen)-sum(b[i] for b in source) for i in range(2))
            E=tuple(sum(c[1][i] for c in chosen)-sum(a[i] for a in target) for i in range(2))
            assert D==tuple(Gamma[i]+E[i] for i in range(2))
            if value:assert min(E)>=0
            if value and pascal:
                assert all(0<=Gamma[i]<=D[i] for i in range(2))
                left=prod(q[i]**D[i]*zeta[i]**Gamma[i]*R**E[i] for i in range(2))
                right=prod((q[i]*R)**E[i]*(q[i]*zeta[i])**Gamma[i] for i in range(2))
                assert left==right
                nonzero+=1
            # Exact row/column factorization used before the analytic bound.
            degree_labels=(2,3,4)
            C=[[coefficient(a,d,c,R) for c in chosen] for a,d in zip(target,degree_labels)]
            normalized=[[C[i][j]*R**sum(target[i])/R**sum(c[1])
                         for j,c in enumerate(chosen)] for i in range(3)]
            assert determinant(C)==R**sum(E)*determinant(normalized)
            coefficient_checks+=1;checked+=1
    assert nonzero>0
    return {'target_column_pairs':checked,'nonzero_joint_terms':nonzero,
            'residual_degree_factorizations':coefficient_checks}


def rational_dimension_bounds():
    records=[]
    for m,bound in ((40,F(451,1000)),(80,F(475,1000)),(120,F(483,1000))):
        lo,hi=quantile_interval(m)
        lower=truncated_mean(m,lo)-F(m,2*(m+2))
        upper=truncated_mean(m,hi)-F(m,2*(m+2))
        assert lower>bound
        records.append({'m':m,'cutoff_interval':[str(lo),str(hi)],
                        'gap_interval':[str(lower),str(upper)],
                        'certified_lower':str(bound)})
    return records


if __name__=='__main__':
    choices,cases=exhaustive_column_cost()
    result={'status':'passed','exhaustive_column_subsets':choices,
            'cheapest_column_cardinalities':cases,'packet_counts':finite_packet_counts(),
            'degree_checks':degree_checks(),'half_density_bounds':rational_dimension_bounds(),
            'scope':'Exact finite quantiles, degree identities, coefficient factorization '
                    'and rational enclosing constants. Asymptotics use the written proofs. '
                    'No actual approximation tuple or pi property is certified.'}
    Path(__file__).with_name('column_cost_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:len(v) if isinstance(v,list) else v for k,v in result.items()},indent=2))
