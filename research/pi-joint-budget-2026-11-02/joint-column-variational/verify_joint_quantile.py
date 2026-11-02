"""Exact checks of joint simplex costs and the degree-constrained frontier.

No finite moment-feasible selection is asserted to give a nonzero minor.
"""
from fractions import Fraction as F
from itertools import product, combinations
from math import comb, factorial, prod
from pathlib import Path
import json


def cdf_and_integral(knots,q):
    d=len(knots)-1
    assert len(set(knots))==len(knots)
    terms=[(max(F(0),q-k),prod(other-k for j,other in enumerate(knots) if j!=i))
           for i,k in enumerate(knots)]
    return (sum(x**d/den for x,den in terms),
            sum(x**(d+1)/den for x,den in terms)/F(d+1))


def clip_triangle(costs,q):
    polygon=[(F(0),F(0)),(F(1),F(0)),(F(0),F(1))]
    value=lambda p:sum(x*c for x,c in zip(p,costs))
    out=[]
    for p,r in zip(polygon,polygon[1:]+polygon[:1]):
        a,b=value(p),value(r);inside_a=a<=q;inside_b=b<=q
        if inside_a!=inside_b:
            t=(q-a)/(b-a)
            out.append(tuple(x+t*(y-x) for x,y in zip(p,r)))
        if inside_b:out.append(r)
    if len(out)<3:return F(0),F(0)
    area=moment=F(0);p=out[0]
    for r,s in zip(out[1:],out[2:]):
        twice=(r[0]-p[0])*(s[1]-p[1])-(r[1]-p[1])*(s[0]-p[0])
        triangle=twice/2
        assert triangle>=0
        area+=triangle;moment+=triangle*(value(p)+value(r)+value(s))/3
    return 2*area,2*moment


def geometric_cdf_checks():
    count=0
    for costs in ((F(1),F(3)),(F(2),F(5)),(F(3),F(-1)),(F(1),F(3,2))):
        knots=(F(0),)+costs
        for q in (F(-2),F(-1),F(-1,2),F(0),F(1,4),F(1,2),F(1),
                  F(3,2),F(2),F(3),F(4),F(5),F(6)):
            probability,moment=clip_triangle(costs,q)
            f,j=cdf_and_integral(knots,q)
            assert f==probability
            assert q*f-j==moment
            count+=1
    return count


def homogeneous_moment(knots,k):
    h=[F(1)]+[F(0)]*k
    for x in knots:
        new=[F(0)]*(k+1)
        for a in range(k+1):
            for b in range(k+1-a):new[a+b]+=h[a]*x**b
        h=new
    d=len(knots)-1
    return F(factorial(k)*factorial(d),factorial(d+k))*h[k]


def integrated_density_moment(knots,k):
    d=len(knots)-1;maximum=max(knots);out=F(0)
    for i,x in enumerate(knots):
        den=prod(other-x for j,other in enumerate(knots) if j!=i)
        length=maximum-x
        out+=F(d,1)/den*sum(F(comb(k,j),d+j)*x**(k-j)*length**(d+j)
                            for j in range(k+1))
    return out


def higher_moment_checks():
    count=0
    all_knots=(F(0),F(1),F(3),F(7),F(11),F(18))
    for d in range(1,6):
        for shift in (F(0),F(-2)):
            knots=tuple(x+shift for x in all_knots[:d+1])
            for k in range(5):
                assert integrated_density_moment(knots,k)==homogeneous_moment(knots,k)
                count+=1
    return count


def quantile_enclosure(knots,density=F(1,2),steps=40):
    lo,hi=min(knots),max(knots)
    for _ in range(steps):
        mid=(lo+hi)/2
        f,_=cdf_and_integral(knots,mid)
        if f<density:lo=mid
        else:hi=mid
    fl,jl=cdf_and_integral(knots,lo)
    fh,jh=cdf_and_integral(knots,hi)
    assert fl<=density<=fh and fl>0
    mean_lo=lo-jl/fl;mean_hi=hi-jh/fh
    assert mean_lo<=mean_hi
    return lo,hi,mean_lo,mean_hi


def columns(weights,H):
    return [a for a in product(*(range(H//w+1) for w in weights))
            if sum(w*x for w,x in zip(weights,a))<=H]


def finite_selection_checks():
    records=[];duals=0
    for weights,H,costs,theta in (
        ((1,2),4,(F(1,2),F(2)),F(1,2)),
        ((1,2,3),4,(F(1),F(3,2),F(2)),F(1,3)),
        ((2,2,3),4,(F(3),F(1),F(2)),F(1,2))):
        raw=columns(weights,H);number=len(raw)//2;m=len(weights)-1
        source=theta/F(m+2)
        means=[]
        for chosen in combinations(raw,number):
            mean=tuple(F(sum(a[i]*weights[i] for a in chosen),number*H)
                       for i in range(m+1))
            means.append(mean)
        feasible=[x for x in means if min(x[1:])>=source]
        assert feasible
        optimum=min(sum(c*x for c,x in zip(costs,mean)) for mean in feasible)
        for eta in product((F(0),F(1,2),F(1),F(2)),repeat=m):
            modified=(costs[0],)+tuple(costs[i+1]-eta[i] for i in range(m))
            dual=min(sum(c*x for c,x in zip(modified,mean)) for mean in means)+sum(eta)*source
            assert dual<=optimum
            duals+=1
        records.append({'weights':weights,'H':H,'columns':len(raw),'selected':number,
                        'subsets':len(means),'degree_feasible_subsets':len(feasible),
                        'exact_feasible_minimum':str(optimum)})
    return records,duals


def frontier_checks():
    count=0;ineligible=0
    for m in range(1,9):
        n=m+1;density=F(1,2)
        for theta in (F(9,10),F(19,20),F(99,100)):
            if theta**n<density:continue
            kappa=F(3,2);stretch=theta**n/density
            assert stretch>=1
            for fraction in (F(1,4),F(1,2),F(1)):
                time_cost=kappa*stretch*fraction
                effective=kappa*fraction
                multiplier=kappa-effective
                probability=theta**n/stretch
                time_mean=theta/(stretch*(m+2))
                spatial_mean=theta/F(m+2)
                primal=time_cost*time_mean+kappa*m*spatial_mean
                quantile_mean=effective*F(n,n+1)*theta
                dual=quantile_mean+multiplier*m*spatial_mean
                claimed=(kappa*m*theta+time_cost*density*theta**(-m))/F(m+2)
                assert probability==density and primal==dual==claimed
                assert primal-kappa*m*spatial_mean==time_cost*density*theta**(-m)/F(m+2)
                count+=1
            assert 2*theta**m>1
            ineligible+=1
    return count,ineligible


if __name__=='__main__':
    geometric=geometric_cdf_checks();moments=higher_moment_checks()
    quantiles=[]
    for knots in ((F(0),F(1),F(3)),(F(0),F(1),F(3),F(7)),
                  (F(0),F(-1),F(2),F(5)),(F(0),F(1),F(2),F(4),F(9))):
        enclosure=quantile_enclosure(knots)
        quantiles.append({'knots':[str(x) for x in knots],
                          'quantile_interval':[str(x) for x in enclosure[:2]],
                          'conditional_mean_interval':[str(x) for x in enclosure[2:]]})
    finite,duals=finite_selection_checks();frontier,ineligible=frontier_checks()
    result={'status':'passed','independent_polygon_cdf_moments':geometric,
            'higher_dimensional_integrated_moments':moments,
            'exact_rational_quantile_enclosures':quantiles,
            'finite_degree_feasible_selections':finite,'finite_dual_inequalities':duals,
            'exact_constrained_frontiers':frontier,'fast_branch_K2_ineligibility_checks':ineligible,
            'scope':'Joint costs and necessary coordinate-degree feasibility only. '
                    'No source/period determinant rank or approximation tuple is certified.'}
    Path(__file__).with_name('joint_quantile_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:len(v) if isinstance(v,list) else v for k,v in result.items()},indent=2))
