"""Independent bounded exact checks of new constructive mechanisms.

No author verifier or LLL module is imported. Asymptotic statements use the
written proofs, not extrapolation from these finite checks.
"""
from fractions import Fraction as F
from math import comb,factorial,gcd,lcm
from pathlib import Path
import json


def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def polyint(a):return sum((F(x,j+1) for j,x in enumerate(a)),F(0))


def evenint(a):return sum((F(x,2*j+1) for j,x in enumerate(a)),F(0))


def moment(m,j,k=0):
    # Direct polynomial expansion, independent of the positive recurrences.
    w=[F(0)]*(m+j)+[F((-1)**r*comb(2*m,r)) for r in range(2*m+1)]
    z=[F((-1)**(k-r)*comb(k,r)*2**r) for r in range(k+1)]
    return evenint(mul(w,z))


def centered(c,K):
    s=[F(0)]*(K+1)
    for k in range(K+1):
        a=F(2,3)**c*F((-1)**k*comb(c+k-1,k),3**k)
        for j in range(k+1):s[j]+=a*comb(k,j)*2**j*(-1)**(k-j)
    return s


def endpoint(c,K):
    s=[F(0)]*(K+1)
    for k in range(K+1):
        a=F(comb(c+k-1,k),2**(c+k))
        for j in range(k+1):s[j]+=a*comb(k,j)*(-1)**j
    return s


def integrated(m,s,j):return sum((v*moment(m,j+l) for l,v in enumerate(s)),F(0))


def oddpart(n):
    while n%2==0:n//=2
    return n


def lcm_to(n):
    a=1
    for j in range(1,n+1):a=lcm(a,j)
    return a


def numerator(R,alpha,c):
    n=[F(0)]*max(2*c+1,len(R)+1)
    for j,r in enumerate(R):
        if j:n[j-1]+=j*r
        n[j+1]+=(j-2*c+2)*r
    for j in range(c):n[2*j]+=alpha*comb(c-1,j)
    return n


def run():
    profile_points=profile_cases=endpoint_cases=telescoping_cases=0
    for m in range(1,5):
        c=4*m+1
        for K in sorted({2*m,4*m,13*m}):
            s=centered(c,K)
            gamma=F(c*2**c*(-1)**K*comb(c+K,K),3**(c+K))
            den=2*K-c+3
            assert den>0
            lhs=[F(0)]*(K+1)
            for j,v in enumerate(s):
                lhs[j]+=c*v
                if j:lhs[j-1]+=j*v;lhs[j]+=j*v
            assert lhs==[gamma*comb(K,j)*2**j*(-1)**(K-j) for j in range(K+1)]
            for r in range(21):
                t=F(r,20)
                value=sum((v*t**j for j,v in enumerate(s)),F(0))
                bound=abs(gamma)*abs(2*t-1)**(K+1)/(den*(1+t))
                assert abs(value-(1+t)**(-c))<=bound
                profile_points+=1
            E=2*((K+1)//2)
            Z=[moment(m,j,E) for j in (0,1)]
            assert 0<Z[1]<Z[0]
            mu=[integrated(m,s,j) for j in (0,1)]
            assert 0<mu[1]<mu[0]
            eps=F(2,3)**c*F(comb(c+K,K+1),3**(K+1))/(1-F(c+K+1,3*(K+2)))
            assert abs(gamma)/(den*eps)==1+F(c-1,(K+2)*den)
            if K==2*m:
                assert Z[0]<=F(7,128)**m
                assert abs(gamma)/den<=F(6*m+1,3)
            if K==13*m:assert Z[0]<=F(1,64)**m
            profile_cases+=1
        # Directly test every positive endpoint truncation in a bounded range.
        for K in range(4*m+1):
            s=endpoint(c,K)
            mu=[integrated(m,s,j) for j in (0,1)]
            assert 0<mu[1]<mu[0]
            expected=[]
            for j in (0,1):
                terms=[]
                for k in range(K+1):
                    prod=1
                    for r in range(2*m+k+1):prod*=2*m+2*j+2*r+1
                    term=F(comb(4*m+k,k)*factorial(2*m+k),2**(2*m+1)*prod)
                    terms.append(term)
                    assert (term.denominator&-term.denominator)<=2**(2*m+1)
                expected.append(sum(terms,F(0)))
            assert mu==expected
            LL=lcm(mu[0].denominator,mu[1].denominator)
            assert (2**(2*m+1)*oddpart(lcm_to(6*m+2*K+3)))%LL==0
            aa=F(comb(2*m,m),64**m)
            alpha=aa*(mu[1]-F(6*m+1,2*m-1)*mu[0])
            assert alpha<0
            if K>=3*m:
                ratio=lambda k,j:F((4*m+k+1)*(2*m+k+1),(k+1)*(6*m+2*j+2*k+3))
                for j in (0,1):assert 0<ratio(K+1,j)<ratio(K,j)<1
            if K==4*m:assert ratio(K,0)<=F(6,7)
            endpoint_cases+=1
        # Full-parity telescoping includes a nonzero R(0) boundary.
        R=[F((j%3)-1) for j in range(2*c)]
        alpha=F(3,2)
        K=2*m;s=centered(c,K)
        sx=[F(0)]*(2*K+1)
        for j,v in enumerate(s):sx[2*j]=v
        n=numerator(R,alpha,c)
        direct=polyint(mul(n,sx))
        g=F(c*2**c*(-1)**K*comb(c+K,K),3**(c+K))
        z=[F(0)]*(2*K+1)
        for j in range(K+1):z[2*j]=(-1)**(K-j)*comb(K,j)*2**j
        p=[F(0)]*(2*c-1)
        for j in range(c):p[2*j]=comb(c-1,j)
        rhs=2*sum(R,F(0))*sum(s,F(0))-R[0]*s[0]
        rhs-=2*g*polyint(mul([F(0)]+R,z))
        rhs+=alpha*polyint(mul(p,sx))
        assert direct==rhs
        assert n[2*c-1]==0
        telescoping_cases+=1
    # Saved parity-free c=5 benchmark, checked directly without search code.
    R=[0,-339,0,-1243,-2,-1630,-78,-745,-164,-59]
    n=numerator(R,339,5)
    hh=[0,-8,83,-286,376,-118,-59]
    assert mul([0,0,1,-2,1],hh)==n
    assert n[9]==0
    assert F(-4*sum(R),16*339)==F(355,113)
    assert sum(R)==-4260
    out=dict(status='passed',profile_cases=profile_cases,profile_rational_point_checks=profile_points,
             positive_endpoint_cases=endpoint_cases,full_parity_boundary_cases=telescoping_cases,
             parity_free_benchmark='355/113 at c=5',
             scope='Bounded exact identity audit; all-scale convergence and period-image statements use the written proofs, not these samples.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':run()
