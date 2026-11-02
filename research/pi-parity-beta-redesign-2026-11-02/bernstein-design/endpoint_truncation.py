"""Exact endpoint-centered kernel generator and beta-moment tail bills."""
from fractions import Fraction as F
from math import comb,factorial,gcd,lcm
from pathlib import Path
import json


def start(m):
    den=1
    for r in range(2*m+1): den*=2*m+2*r+1
    return F(factorial(2*m),2**(2*m+1)*den)


def ratio(m,k,j=0):
    return F((4*m+k+1)*(2*m+k+1),(k+1)*(6*m+2*j+2*k+3))


def moments(m,K):
    term=start(m)
    mu0=mu1=F(0)
    for k in range(K+1):
        mu0+=term
        mu1+=term*F(2*m+1,6*m+2*k+3)
        term*=ratio(m,k)
    return mu0,mu1,term


def machin(n=140):
    def atan(q):
        ss=F(0)
        for k in range(n): ss+=F((-1)**k,(2*k+1)*q**(2*k+1))
        nxt=F((-1)**n,(2*n+1)*q**(2*n+1))
        return min(ss,ss+nxt),max(ss,ss+nxt)
    aa=atan(5);bb=atan(239)
    return 16*aa[0]-4*bb[1],16*aa[1]-4*bb[0]


def independent_checks():
    checks=0
    for m in range(1,5):
        c=4*m+1
        for K in range(7):
            ss=[F(0)]*(K+1)
            for k in range(K+1):
                v=F(comb(c+k-1,k),2**(c+k))
                for r in range(k+1):ss[r]+=v*(-1)**r*comb(k,r)
            nn=[F(0)]*(m+2*m+K+1)
            for i in range(2*m+1):
                for r,v in enumerate(ss):nn[m+i+r]+=(-1)**i*comb(2*m,i)*v
            expected=[sum((F(v,2*(r+j)+1) for r,v in enumerate(nn)),F(0)) for j in (0,1)]
            got=moments(m,K)
            assert expected==list(got[:2]);checks+=2
    return checks


def case(m,K):
    mu0,mu1,next0=moments(m,K)
    aa=F(comb(2*m,m),64**m)
    rr=F(6*m+1,2*m-1)
    bb=F(1,(2*m-1)*comb(4*m,2*m))
    alpha=aa*(mu1-rr*mu0)
    beta=bb*mu0
    assert alpha<0 and mu0>mu1>0
    cand=-4*beta/alpha;q=cand.denominator
    LL=lcm(mu0.denominator,mu1.denominator)
    # No per-term dyadic denominator growth after beta integration.
    assert all((v.denominator & -v.denominator)<=2**(2*m+1) for v in (mu0,mu1))
    cap=1
    for z in range(1,6*m+2*K+4):cap=lcm(cap,z)
    while cap%2==0:cap//=2
    assert (2**(2*m+1)*cap)%LL==0
    out=dict(m=m,K=K,p=str(cand.numerator),q=str(q),q_bits=q.bit_length(),
             moment_denominator_bits=LL.bit_length(),mu0=str(mu0),mu1=str(mu1),
             nonzero_period_proved=True)
    r0=ratio(m,K+1);r1=ratio(m,K+1,1)
    if r0<1:
        tail0=next0/(1-r0)
        next1=next0*F(2*m+1,6*m+2*(K+1)+3)
        tail1=next1/(1-r1)
        # A monotone likelihood-ratio identity proves the signed integral is
        # positive, so dropping its negative second tail is a valid bound.
        integral_upper=mu1*tail0
        quality_upper=4*q*q*integral_upper/abs(alpha)
        out.update(remainder_quality_upper=str(quality_upper),displayed_remainder_upper=float(quality_upper))
    lo,hi=machin()
    assert cand>hi
    errlo=min(abs(lo-cand),abs(hi-cand)) if not lo<=cand<=hi else F(0)
    errhi=max(abs(lo-cand),abs(hi-cand))
    out.update(displayed_Q_lower=float(q*q*errlo),displayed_Q_upper=float(q*q*errhi))
    if 'remainder_quality_upper' in out: assert q*q*errhi<=F(out['remainder_quality_upper'])
    return out


if __name__=='__main__':
    count=independent_checks()
    rows=[case(m,K) for m in range(1,7) for K in (m,2*m,3*m,4*m)]
    result=dict(status='passed',construction='positive Taylor expansion at t=1; exact beta moments',
                independent_polynomial_moment_checks=count,
                scope='Bounded design control, not a broad candidate search. No target fraction or pi enters generation; Machin intervals only certify outputs.',cases=rows)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',prescribed_pairs=len(rows),independent_moment_checks=count,
                         smallest_displayed_actual_quality=min(r['displayed_Q_lower'] for r in rows)),indent=2))
