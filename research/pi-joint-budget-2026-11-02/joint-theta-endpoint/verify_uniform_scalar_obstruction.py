"""Exact checks for the uniform scalar mean/collision obstruction.

Only rational arithmetic is used. Positivity of the actual determinant,
time-grade matching, or existence of approximation tuples is not checked.
"""
from fractions import Fraction as F
from itertools import product, combinations
from pathlib import Path
import json


def moment(m,t):
    a=F(m,m+2);mu=F(m,m+1)
    variance=F(m,(m+1)**2*(m+2))
    if t<=a:
        return F(2*(m+1),m+2)*(a/t)**m
    if t<=mu:
        return 1+(mu-t)**2/variance
    return F(1)


def collision(m,w0,c,t):
    return c*moment(m,t)/F(2*(m+1)*w0)


def envelope(m,theta,K,w0,c,u,y,t):
    s=theta*F(m,m+2)
    v0=2*K*theta**m*w0
    kappa=1+2*c/(theta*v0)
    return 14*K*y/w0+kappa*(u-s)-2*(t-s)-collision(m,w0,c,t)


def low_witness(m,root_theta,K,w0,root_ratio,c):
    theta=root_theta**2
    a=F(m,m+2);s=a*theta
    u=a*root_ratio**2
    t=a*root_theta*root_ratio
    assert s<=t<=u<=a
    y=moment(m,u)/F(4*(m+1))
    value=envelope(m,theta,K,w0,c,u,y,t)
    z=collision(m,w0,c,u)
    x=F(m+2)*(u-s)/theta
    r=(root_ratio/root_theta)**m
    assert r*r==(1+x/m)**m
    formula=a*(root_ratio-root_theta)**2+z*(7*K/c+r*r*x/K-r)
    assert value==formula
    assert value>0
    if x>=F(1,4):
        assert x*r*r/K-r>=-F(K)/(4*x)>=-K
    else:
        assert r<2
    return value


def high_witness(m,root_theta,K,w0,u,c):
    theta=root_theta**2
    a=F(m,m+2);s=a*theta;t=a*root_theta
    assert s<=t<=a<=u<=1
    y=F(1,4*(m+1))
    value=envelope(m,theta,K,w0,c,u,y,t)
    p=K*root_theta**m
    rho=1/(w0*root_theta**m)
    x=F(m)*(1-theta)/theta
    lower=(u-a)+a*(1-root_theta)**2+rho/F(m+2)*(
        F(7,2)*F(m+2,m+1)*p+c*x/p-c)
    # The difference is the nonnegative excess LCM cost above u=a.
    v0=2*K*theta**m*w0
    assert value-lower==2*c/(theta*v0)*(u-a)
    assert value>=lower>0
    return value


def witness_checks():
    low=high=0
    # c=ln2/2 lies strictly between these rationals. The expression is
    # affine in c at fixed data, so both endpoint checks cover it.
    for m in (1,2,3,5,10,20,40,80):
        a=F(m,m+2)
        for root in (F(1,4),F(1,2),F(3,4),F(9,10),F(99,100),F(1)):
            for K in (1,2,7):
                for w0 in (F(1),F(5,3)**m):
                    for c in (F(1,3),F(7,20)):
                        for mixing in (F(0),F(1,2),F(1)):
                            ratio=root+(1-root)*mixing
                            low_witness(m,root,K,w0,ratio,c)
                            low+=1
                        for mixing in (F(0),F(1,2),F(1)):
                            u=a+(1-a)*mixing
                            high_witness(m,root,K,w0,u,c)
                            high+=1
    return {'low_mean_endpoint_checks':low,'high_mean_endpoint_checks':high}


def types(weights,H):
    return [a for a in product(*(range(H//w+1) for w in weights))
            if sum(w*x for w,x in zip(weights,a))<=H]


def columns(w0,weights,H):
    return [(h,a) for a in types(weights,H)
            for h in range((H-sum(w*x for w,x in zip(weights,a)))//w0+1)]


def finite_column_inequalities():
    count=duals=0
    for w0,weights,H in ((1,(2,),3),(2,(2,3),4),(1,(2,3),3)):
        all_types=types(weights,H)
        cols=columns(w0,weights,H)
        costs={a:F(sum(w*x for w,x in zip(weights,a)),H) for a in all_types}
        for M in range(1,len(cols)+1):
            for chosen in combinations(cols,M):
                multiplicity={a:0 for a in all_types}
                for h,a in chosen:multiplicity[a]+=1
                hsum=sum(h for h,_ in chosen)
                norm=sum(n*n for n in multiplicity.values())
                assert 2*hsum>=norm-M
                mean=sum(costs[a]*n for a,n in multiplicity.items())/M
                # Finite dual tests used in passing from column counts to
                # the continuum scalar collision inequality.
                for cutoff in (F(1,2),F(1),F(3,2)):
                    if cutoff<=mean:continue
                    denominator=sum(max(F(0),cutoff-x)**2 for x in costs.values())
                    dual=M*M*(cutoff-mean)**2/denominator
                    assert norm>=dual
                    assert 2*hsum>=dual-M
                    duals+=1
                count+=1
    return {'selected_column_subsets':count,'finite_column_moment_duals':duals}


def selected_column_moments():
    records=[]
    for weights in ((2,),(2,3),(2,3,5)):
        m=len(weights);w0=1
        for H in (24,48,96):
            cols=columns(w0,weights,H)
            M=len(cols)//2
            selectors={
                'cheapest_time':lambda col:(col[0],sum(col[1])),
                'cheapest_transverse':lambda col:(sum(w*x for w,x in zip(weights,col[1])),col[0]),
                'mixed_cost':lambda col:(3*col[0]+sum(w*x for w,x in zip(weights,col[1])),col[0])}
            for label,key in selectors.items():
                chosen=sorted(cols,key=key)[:M]
                counts={}
                for h,a in chosen:counts[a]=counts.get(a,0)+1
                hsum=sum(h for h,_ in chosen)
                square=sum(n*n for n in counts.values())
                assert 2*hsum>=square-M
                u=F(sum(sum(w*x for w,x in zip(weights,a)) for _,a in chosen),M*H)
                y=F(w0*hsum,M*H)
                records.append({'m':m,'H':H,'selection':label,'selected_columns':M,
                                'all_columns':len(cols),'transverse_mean':str(u),
                                'time_mean':str(y),
                                'exact_time_from_multiplicity_lower':str(F(w0*(square-M),2*M*H)),
                                'continuum_half_density_lower_at_observed_mean':str(moment(m,u)/F(4*(m+1)))})
    return records


def joint_sequence_checks():
    records=[]
    for m in (4,8,16,32,64,128):
        for exponent in (1,2):
            root=1-F(1,m**exponent)
            theta=root*root
            for K in (1,m):
                w0=F(m+1,m)**m
                ratio=(1+root)/2
                low_values=[low_witness(m,root,K,w0,ratio,c)
                            for c in (F(1,3),F(7,20))]
                high_values=[high_witness(m,root,K,w0,F(m,m+2),c)
                             for c in (F(1,3),F(7,20))]
                records.append({'m':m,'root_theta':str(root),'theta':str(theta),'K':K,
                                'w0':str(w0),'lower_enclosed_low_witness':str(min(low_values)),
                                'lower_enclosed_high_witness':str(min(high_values))})
    return records


if __name__=='__main__':
    result={'status':'passed','scalar_witness_checks':witness_checks(),
            'finite_column_checks':finite_column_inequalities(),
            'selected_column_moments':selected_column_moments(),
            'joint_theta_sequences':joint_sequence_checks(),
            'scope':'Same-column scalar mass/mean relaxation only. '
                    'No arbitrary-column coordinate matching, time-grade matching, '
                    'actual-period nonvanishing, or pi assertion is certified.'}
    Path(__file__).with_name('uniform_scalar_verification.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:len(v) if isinstance(v,list) else v for k,v in result.items()},indent=2))
