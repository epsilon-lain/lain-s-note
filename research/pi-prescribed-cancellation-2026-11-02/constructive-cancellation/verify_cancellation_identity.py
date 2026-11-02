"""Exact algebra audit; multiplier generation does not use pi."""
from fractions import Fraction as F
from math import comb, gcd, lcm, prod
from pathlib import Path
import json


def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j]+=x*y
    return out


def falling(x,k):
    return prod((x-j for j in range(k)),start=F(1))


def rodrigues(a,b,k,p):
    out=[F(0)]*(len(p)+k)
    for ell,coef in enumerate(p):
        for s in range(k+1):
            fac=coef*comb(k,s)*falling(F(2*a+2*k+2*ell-1,2),k-s)*(-1)**s*falling(F(b+k),s)
            for j in range(k-s+1):
                out[ell+s+j]+=fac*(-1)**j*comb(k-s,j)
    return out


def numerator(a,b,h):
    return [F(0)]*a+mul([F((-1)**j*comb(b,j)) for j in range(b+1)],h)


def period(n,c):
    # Independent coefficient-by-coefficient rational primitive reduction.
    n=[F(v) for v in n]+[F(0)]*(c+1-len(n))
    rr=[]
    for j in range(c):
        A=n[j];B=F(-comb(c-1,j))
        if j:
            A-=(2*j-2*c+1)*rr[-1][0]
            B-=(2*j-2*c+1)*rr[-1][1]
        rr.append((A/(2*j+1),B/(2*j+1)))
    alpha=(n[c]-rr[-1][0])/rr[-1][1]
    r=[x+alpha*y for x,y in rr]
    # Reconstruct every numerator coefficient, including the top one.
    for j in range(c+1):
        val=(2*j+1)*(r[j] if j<c else 0)
        val+=(2*j-2*c+1)*(r[j-1] if j else 0)
        val+=alpha*(comb(c-1,j) if j<c else 0)
        assert val==n[j]
    beta=sum(r)/2**(c-1)
    return alpha,beta


def weighted_moment(a,b,h,j):
    return sum(coef*(-1)**s*comb(b,s)/F(2*(a+j+ell+s)+1,2)
               for ell,coef in enumerate(h) for s in range(b+1))


def laurent_mul(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():out[i+j]=out.get(i+j,0)+x*y
    return out


def fourier_period(n,c):
    n=list(n)+[0]*(c+1-len(n))
    assert all(F(v).denominator==1 for v in n)
    n=[int(v) for v in n]
    top=n[c];nn=c-1;D=4**nn;L=lcm(*range(1,nn+1))
    row={}
    for j in range(c):
        coef=n[j]-top*comb(c,j)
        if not coef:continue
        term={0:1}
        for _ in range(j):term=laurent_mul(term,{-1:-1,0:2,1:-1})
        for _ in range(nn-j):term=laurent_mul(term,{-1:1,0:2,1:1})
        for k,v in term.items():row[k]=row.get(k,0)+coef*v
    assert all(row.get(k,0)==row.get(-k,0) for k in range(nn+1))
    A=row.get(0,0)
    B=2*D*L*top
    sine=(0,1,0,-1)
    for k in range(1,nn+1):B+=(L//k)*2*row.get(k,0)*sine[k%4]
    return F(A,D),F(B,2*D*L),A,B,L


def check():
    moment_checks=period_checks=nondegenerate=0
    for m in range(1,7):
        a,b,c,d=m,2*m,4*m+1,m+1
        for k in range(d+1):
            residuals=[[F(0)]*ell+[F(1)] for ell in range(d-k+1)]
            residuals.append([F((-1)**ell*(ell+1)) for ell in range(d-k+1)])
            for p in residuals:
                h=rodrigues(a,b,k,p)
                assert len(h)-1<=d
                assert all((2**k*x).denominator==1 for x in h)
                for j in range(k):
                    assert weighted_moment(a,b,h,j)==0
                    moment_checks+=1
                old=period(numerator(a,b,h),c)
                new=period(numerator(a+k,b+k,p),c+k)
                factor=prod(range(c,c+k))
                assert old==tuple(factor*x for x in new)
                period_checks+=1
                if old[0]:
                    assert -4*old[1]/old[0]==-4*new[1]/new[0]
                    nondegenerate+=1
    constant_checks=0
    for m in range(1,13):
        alpha,beta=period(numerator(m,2*m,[F(1)]),4*m+1)
        assert alpha==F(comb(2*m,m),64**m) and beta==0
        constant_checks+=1
    fourier_checks=0
    for c in range(3,18):
        for offset in range(3):
            n=[(-1)**(j+offset)*(j+offset+1) for j in range(c+1)]
            exact=period(n,c)
            alpha,beta,A,B,L=fourier_period(n,c)
            assert exact==(alpha,beta)
            if A:
                frac=F(-2*B,L*A)
                assert frac==-4*beta/alpha
                assert frac.denominator==abs(L*A)//gcd(2*B,L*A)
            fourier_checks+=1
    return dict(status='passed',maximum_family_m=6,
                exact_rational_moment_checks=moment_checks,
                shifted_period_identity_checks=period_checks,
                nonzero_period_reduced_fraction_equalities=nondegenerate,
                constant_kernel_mass_checks=constant_checks,
                independent_Fourier_period_checks=fourier_checks,
                scope='Exact bounded algebra audit. No quality-zero sequence or asymptotic coefficient bound is inferred.')


if __name__=='__main__':
    out=check()
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
