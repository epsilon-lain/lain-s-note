"""Independent bounded audit of the explicit K=13m truncation generator."""
from fractions import Fraction as F
from math import comb,factorial,gcd,lcm
from pathlib import Path
import json
from audit_rodrigues import endpoint,mul,period


def epsilon(c,k):
    return F(2,3)**c*F(comb(c+k,k+1),3**(k+1))/(1-F(c+k+1,3*(k+2)))


def beta_moment(a,b,j):
    den=1
    for r in range(b+1):den*=2*a+2*j+2*r+1
    return F(2**b*factorial(b),den)


def truncated_polynomial(c,k):
    out=[F(0)]*(k+1)
    for r in range(k+1):
        factor=F(2,3)**c*F((-1)**r*comb(c+r-1,r),3**r)
        for j in range(r+1):out[j]+=factor*comb(r,j)*2**j*(-1)**(r-j)
    return out


def integrate_polynomial(n):
    return sum((F(x,2*j+1) for j,x in enumerate(n)),F(0))


def odd_lcm(n):
    x=1
    for j in range(1,n+1):x=lcm(x,j)
    while x%2==0:x//=2
    return x


def run():
    rows=[]
    for m in range(1,4):
        a,b,c,k=m,2*m,4*m+1,13*m
        n0=endpoint(a,b);n1=[F(0)]+n0
        beta0,alpha0=period(n0,c);beta1,alpha1=period(n1,c)
        assert beta0==0
        assert alpha0==F(comb(2*m,m),64**m)
        assert alpha1==alpha0*F(6*m+1,2*m-1)
        assert beta1==-F(1,(2*m-1)*comb(4*m,2*m))
        delta=alpha0*beta1-alpha1*beta0
        assert delta
        M0,M1=beta_moment(a,b,0),beta_moment(a,b,1)
        eps=epsilon(c,k)
        assert eps*(abs(alpha0)*M1+abs(alpha1)*M0)<=abs(delta)/2
        s=truncated_polynomial(c,k)
        mu0=sum((v*beta_moment(a,b,j) for j,v in enumerate(s)),F(0))
        mu1=sum((v*beta_moment(a,b,j+1) for j,v in enumerate(s)),F(0))
        h=[mu1,-mu0]
        assert integrate_polynomial(mul(mul(n0,h),s))==0
        beta,alpha=period(mul(n0,h),c)
        assert (beta,alpha)==(mu1*beta0-mu0*beta1,mu1*alpha0-mu0*alpha1)
        assert abs(alpha-delta)<=eps*(abs(alpha0)*M1+abs(alpha1)*M0)
        assert alpha
        assert abs(delta)/2<=abs(alpha)<=3*abs(delta)/2
        L=lcm(mu0.denominator,mu1.denominator)
        D0=lcm(*(x.denominator for x in [alpha0,alpha1,beta0,beta1]))
        assert (3**(c+k)*odd_lcm(2*(a+b+k+1)+1))%L==0
        assert (4**(c-1)*odd_lcm(c-1))%D0==0
        A,B=L*D0*alpha,L*D0*beta
        assert A.denominator==B.denominator==1
        candidate=-4*beta/alpha
        u0,u1=L*mu0,L*mu1
        explicit=F(2**(6*m+2)*u0,comb(4*m,2*m)*comb(2*m,m)*((6*m+1)*u0-(2*m-1)*u1))
        assert candidate==explicit
        q=candidate.denominator
        assert q==abs(A.numerator)//gcd(abs(A.numerator),4*abs(B.numerator))
        assert q<=L*D0*abs(alpha)
        assert q<=F(3,2)*comb(2*m,m)*L
        integral_upper=eps*(abs(mu1)*M0+abs(mu0)*M1)
        absolute_error_upper=4*integral_upper/abs(alpha)
        assert absolute_error_upper<=F(128,3)*(4*m*m-1)*F(2**33,3**23)**m
        quality_upper=q*q*absolute_error_upper
        coarse_upper=12*(L*D0)**2*abs(delta)*eps*(1+eps)*M0*M1
        assert quality_upper<=coarse_upper
        # Independently recover the rational primitive and clear its integer map.
        n=mul(n0,h);n+= [F(0)]*(c+1-len(n))
        rr=[]
        for j in range(c):
            prev=rr[-1] if rr else F(0)
            rr.append((n[j]-(2*j-2*c+1)*prev-alpha*comb(c-1,j))/ (2*j+1))
        assert rr[-1]==n[c]
        assert sum(rr,F(0))/2**(c-1)==beta
        scale=2**(c-1)*factorial(c-1)*L
        assert all((scale*x).denominator==1 for x in rr+[alpha])
        rows.append(dict(m=m,K=k,p=str(candidate.numerator),q=str(q),
                         q_bits=q.bit_length(),multiplier_denominator_bits=L.bit_length(),
                         alpha=str(alpha),delta=str(delta),
                         absolute_error_upper=str(absolute_error_upper),
                         absolute_error_upper_below_uniform_bound=True,
                         quality_upper=str(quality_upper),
                         quality_upper_less_than_one=quality_upper<1,
                         exact_integer_realization_scale=str(scale)))
    out=dict(status='passed',scales_checked=[1,2,3],rule='K=13m',samples=rows,
             scope='Finite exact generator and cost audit. The written proof establishes absolute-error convergence and divergent reduced denominators; neither it nor this finite check establishes quality tending to zero.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='samples'},indent=2))


if __name__=='__main__':run()
