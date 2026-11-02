"""Independent exact reconstruction and error-bound checks; standard library."""
from fractions import Fraction as F
from pathlib import Path
from math import gcd,lcm
import json

ROOT=Path(__file__).parent

def coefs(n):
    c=F(1);out=[]
    for k in range(n):
        out.append(c)
        c*=F(2*(k+1),2*k+3)
    return out

def approximant(x,n):
    rho=x*x/(1+x*x)
    h=x/(1+x*x)
    total=F()
    power=F(1)
    for c in coefs(n):
        total+=c*power
        power*=rho
    return h*total

def atan_bounds(k,n):
    p=F();power=F(1,k);r=F(1,k*k)
    for j in range(n):
        p+=(-1)**j*power/(2*j+1)
        power*=r
    b=p+(-1)**n*power/(2*n+1)
    return min(p,b),max(p,b)

def qbounds(fr,pi):
    lo,hi=pi
    e0=F() if lo<=fr<=hi else min(abs(fr-lo),abs(fr-hi))
    e1=max(abs(fr-lo),abs(fr-hi))
    return e0*fr.denominator**2,e1*fr.denominator**2

def main():
    data=json.loads((ROOT/'results.json').read_text())
    a5,a239=atan_bounds(5,91),atan_bounds(239,31)
    pi=(16*a5[0]-4*a239[1],16*a5[1]-4*a239[0])
    assert pi[0]<pi[1]
    saved_pi=tuple(map(F,data['pi_interval']))
    assert saved_pi[0]<=pi[0]<pi[1]<=saved_pi[1]
    improvements=0
    for row in data['rows']:
        a,b,n,L=[row[k] for k in ['a','b','n','delta_terms']]
        D=a*a+b*b;m=b*b+2*a*b-a*a;k=a*a+2*a*b-b*b
        assert abs(m)==1 and m*m+k*k==2*D*D
        assert gcd(a,b)==1 and (a-b)%2==1
        r,s,delta=F(b,a),F(a-b,a+b),F(m,k)
        zr,zs,zd=[x*x/(1+x*x) for x in [r,s,delta]]
        assert zr-zs==F(m,2*D) and zd==F(1,2*D*D)
        A=2**(n+1)*a*b**(2*n+1)
        B=(a*a-b*b)*(a-b)**(2*n)
        assert A!=B and gcd(A,B)==1
        if n==1: assert B-A==-m*D
        if n%4==1: assert (B-A)%D==0
        ar,ass,ad=[approximant(x,j) for x,j in [(r,n),(s,n),(delta,L)]]
        rr,rs=8*ar-4*ad,8*ass+4*ad
        C=(B*rr-A*rs)/(B-A)
        naive=4*(ar+ass)
        assert str(C)==row['fraction'] and str(naive)==row['naive_fraction']
        assert C.numerator==int(row['p']) and C.denominator==int(row['q'])
        d=lcm(rr.denominator,rs.denominator)
        P=B*int(rr*d)-A*int(rs*d);Q=d*(B-A)
        assert gcd(abs(P),abs(Q))==int(row['gcd'])
        H=lcm(*(c.denominator for c in coefs(max(n,L))))
        budget=H*(2*D)**max(n,2*L)
        assert budget%rr.denominator==0 and budget%rs.denominator==0
        # General synchronised tail bound, evaluated as a rational certificate.
        cn=coefs(n+1)[-1];cl=coefs(L+1)[-1]
        R=max(zr,zs)
        kr=r/(1+r*r)*zr**n
        hd=abs(delta/(1+delta*delta))
        bound=8*cn*B*kr/abs(B-A)/(2*D*(1-R)**2)
        bound+=4*(A+B)/abs(B-A)*hd*cl*zd**L/(1-zd)
        actual=max(abs(C-pi[0]),abs(C-pi[1]))
        assert actual<bound
        qq=qbounds(C,pi);qn=qbounds(naive,pi)
        improved=qq[1]<qn[0]
        assert improved==row['strictly_better_than_naive']
        improvements+=improved
        display=F(row['quality_display'])
        assert abs(display-sum(qq)/2)<max(F(1),qq[1])*F(1,10**15)
        assert qbounds(F(355,113),pi)[1]<qq[0]
    assert len(data['rows'])==336 and improvements==1
    report={'status':'PASS','independent_candidates':336,'shared_norm_checks':336,
            'actual_gcd_and_denominator_budget_checks':336,
            'synchronised_tail_error_bound_certificates':336,
            'quality_improvements_over_own_naive_baseline':improvements,
            'all_candidates_worse_than_355_over_113':True}
    (ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
