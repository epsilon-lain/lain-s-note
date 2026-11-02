"""Exact Taylor differential residual and weighted remainder certificates."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import sys

for ancestor in Path(__file__).resolve().parents:
    for candidate in [ancestor/'constructive-review',
                      ancestor/'pi-prescribed-cancellation-2026-11-02/constructive-review']:
        if (candidate/'audit_complete_truncation.py').exists():
            sys.path.insert(0,str(candidate))
            break
    else:
        continue
    break
else:
    raise RuntimeError('The preceding prescribed-cancellation checkpoint is required.')
from audit_complete_truncation import epsilon, beta_moment, truncated_polynomial
from audit_rodrigues import endpoint, mul, period


def polynomial_integral(p):
    return sum((F(v,2*j+1) for j,v in enumerate(p)),F(0))


def residual(c,K):
    return F(c*2**c*(-1)**K*comb(c+K,K),3**(c+K))


def shifted_even_moment(a,b,j,E):
    assert E%2==0
    return sum((F((-1)**(E-l)*comb(E,l)*2**l)*beta_moment(a,b,j+l)
                for l in range(E+1)),F(0))


def centered_moment(a,b,j,k):
    return sum((F((-1)**(k-l)*comb(k,l)*2**l)*beta_moment(a,b,j+l)
                for l in range(k+1)),F(0))


def check(m):
    a,b,c,K=m,2*m,4*m+1,13*m
    s=truncated_polynomial(c,K)
    g=residual(c,K)
    # Differential residual is a single highest-degree polynomial term.
    actual=[F(0)]*(K+1)
    for j,v in enumerate(s):
        actual[j]+=c*v
        if j:
            actual[j-1]+=j*v
            actual[j]+=j*v
    target=[g*comb(K,j)*2**j*(-1)**(K-j) for j in range(K+1)]
    assert actual==target
    # Independently check the proved pointwise inequality at rational points.
    for j in range(21):
        t=F(j,20)
        sk=sum((v*t**i for i,v in enumerate(s)),F(0))
        exact=(1+t)**(-c)
        pointbound=abs(g)*abs(2*t-1)**(K+1)/((2*K-c+3)*(1+t))
        assert abs(sk-exact)<=pointbound
    mu0=sum((v*beta_moment(a,b,j) for j,v in enumerate(s)),F(0))
    mu1=sum((v*beta_moment(a,b,j+1) for j,v in enumerate(s)),F(0))
    increments=[]
    for j,mu in [(0,mu0),(1,mu1)]:
        ww=[centered_moment(a,b,j,k) for k in range(K+1)]
        ff=[F(2,3)**c*F((-1)**k*comb(c+k-1,k),3**k)*ww[k]
            for k in range(K+1)]
        assert all(x>0 for x in ff)
        assert sum(ff,F(0))==mu
        increments.append(ff)
        for k in range(1,K):
            assert (2*k+6*m+2*j+3)*ww[k+1]==-(2*m-2*j+1)*ww[k]+2*k*ww[k-1]
            assert 9*(k+1)*(2*k+6*m+2*j+3)*ff[k+1]==3*(2*m-2*j+1)*(c+k)*ff[k]+2*(c+k)*(c+k-1)*ff[k-1]
    prefix0,prefix1=F(0),F(0)
    for k in range(K+1):
        assert increments[0][k]>increments[1][k]>0
        prefix0+=increments[0][k]
        prefix1+=increments[1][k]
        assert 0<prefix1<prefix0
        assert prefix1-F(6*m+1,2*m-1)*prefix0<0
    # Audit the weaker D>0 assumption and the K=2m profile as well.
    shortK=2*m
    short=truncated_polynomial(c,shortK)
    shortg=residual(c,shortK)
    for j in range(21):
        t=F(j,20)
        sk=sum((v*t**i for i,v in enumerate(short)),F(0))
        assert abs(sk-(1+t)**(-c))<=abs(shortg)*abs(2*t-1)**(shortK+1)/((2*shortK-c+3)*(1+t))
    for j in (0,1):
        assert shifted_even_moment(a,b,j,shortK)<=F(49,864)**m
    assert abs(shortg)<4*m+1
    h=[mu1,-mu0]
    beta,alpha=period(mul(endpoint(a,b),h),c)
    assert alpha
    candidate=-4*beta/alpha
    q=candidate.denominator
    # Telescoping of the denominator combination to one high Taylor residue.
    S,R=2*m-1,6*m+1
    target_n=mul(endpoint(a,b),[R,-S])
    target_n+=[F(0)]*(c+1-len(target_n))
    rr=[]
    for j in range(c):
        rr.append((target_n[j]-(2*j-2*c+1)*(rr[-1] if rr else 0))/(2*j+1))
    assert rr[-1]==target_n[c]
    assert sum(rr)==F(2**(c-1),comb(4*m,2*m))
    moment_poly=[F(0)]+rr
    power_poly=[F(comb(K,j)*2**j*(-1)**(K-j)) for j in range(K+1)]
    top=polynomial_integral(mul(moment_poly,power_poly))
    value=R*mu0-S*mu1
    assert value==2*sum(rr)*sum(s)-2*g*top
    E=2*((K+1)//2)
    z0=shifted_even_moment(a,b,0,E)
    z1=shifted_even_moment(a,b,1,E)
    M0,M1=beta_moment(a,b,0),beta_moment(a,b,1)
    assert 0<z0<M0 and 0<z1<M1
    old_error=epsilon(c,K)*(abs(mu1)*M0+abs(mu0)*M1)
    new_error=abs(g)* (abs(mu1)*z0+abs(mu0)*z1) / (2*K-c+3)
    assert new_error<old_error
    assert abs(g)/((2*K-c+3)*epsilon(c,K))==1+F(c-1,(K+2)*(2*K-c+3))
    assert z0<=F(1,64**m) and z1<=F(1,64**m)
    return dict(m=m,K=K,E=E,q_bits=q.bit_length(),
                denominator_combination_bits=value.denominator.bit_length(),
                positive_increment_checks=2*(K+1),
                all_prefix_nonzero_period_checks=K+1,
                pointwise_rational_checks=42,telescoping_identity=True,
                old_quality_upper=str(4*q*q*old_error/abs(alpha)),
                profiled_quality_upper=str(4*q*q*new_error/abs(alpha)),
                displayed_upper_improvement=float(old_error/new_error),
                displayed_profiled_quality_upper=float(4*q*q*new_error/abs(alpha)),
                profiled_upper_below_one=4*q*q*new_error<abs(alpha))


if __name__=='__main__':
    rows=[check(m) for m in range(1,7)]
    result=dict(status='passed',scales=list(range(1,7)),samples=rows,
                scope='Exact rational comparison of remainder certificates on the fixed generated fractions; no candidate search or record claim.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in r.items() if not k.endswith('quality_upper')} for r in rows],indent=2))
