"""Exact period extraction and joint certificates for prescribed moment cancellation."""
from fractions import Fraction as F
from math import comb, gcd, lcm, factorial
from functools import reduce
from pathlib import Path
import json


def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return out


def falling(a,k):
    out=1
    for j in range(k): out*=a-j
    return out


def jacobi(m,L):
    # 2^L times the Rodrigues derivative; clear its full integer content.
    h=[0]*(L+1)
    for j in range(L+1):
        coeff=(-1)**j*comb(L,j)*2**j*falling(2*m+L,j)
        for k in range(L-j): coeff*=2*m+2*L-1-2*k
        for k in range(L-j+1):
            h[j+k]+=coeff*(-1)**k*comb(L-j,k)
    content=reduce(gcd,map(abs,h))
    return [x//content for x in h],F(2**L,content)


def extract(n,c):
    n=n+[0]*(c+1-len(n))
    nu=[sum(n[k]*comb(k,j)*(-1)**(k-j) for k in range(j,c+1)) for j in range(c+1)]
    rho=[]
    for j in range(c-1):
        prev=rho[-1] if rho else 0
        rho.append(F(nu[j]-(2*j-2*c+1)*prev,2*(c-j-1)))
    alpha=nu[c-1]+rho[-1]
    rho.append(F(nu[c]))
    rr=[sum(rho[j]*comb(j,k) for j in range(k,c)) for k in range(c)]
    r1=sum(rr)
    scale=lcm(*(x.denominator for x in rr),alpha.denominator)
    integers=[int(x*scale) for x in rr]+[int(alpha*scale)]
    content=reduce(gcd,map(abs,integers))
    integers=[x//content for x in integers]
    # Independently audit monomial differential identity.
    r=integers[:-1]; ai=integers[-1]
    ni=[(2*j+1)*(r[j] if j<c else 0)+(2*j-2*c+1)*(r[j-1] if j else 0)+ai*(comb(c-1,j) if j<c else 0) for j in range(c+1)]
    multiplier=F(scale,content)
    assert all(F(ni[j])==multiplier*n[j] for j in range(c+1))
    if alpha:
        rational=-r1/(2**(c-3)*alpha)
    else: rational=None
    return dict(alpha=alpha,r1=r1,rational=rational,integer_r=r,integer_alpha=ai,
                numerator_multiplier=multiplier,numerator_integer=ni)


def data(m):
    L=m+1;c=4*m+1
    h,K=jacobi(m,L)
    n=[0]*m+mul([(-1)**j*comb(2*m,j) for j in range(2*m+1)],h)
    d=extract(n,c)
    A=m+L;B=2*m+L;C=c+L
    ns=[0]*A+[(-1)**j*comb(B,j) for j in range(B+1)]
    ds=extract(ns,C)
    factor=K*falling(c+L-1,L)
    assert d['alpha']==factor*ds['alpha']
    # Equality of complete rational-plus-pi integrals after integration by parts.
    assert d['r1']/2**(c-1)==factor*ds['r1']/2**(C-1)
    assert d['rational']==ds['rational']
    assert d['alpha']!=0
    q=d['rational'].denominator
    shape=F(A**A*B**B,2**A*(A+B)**(A+B))
    upper=4*q*q*shape/abs(ds['alpha'])
    # Universal exponential-size denominator certificate by trig Fourier expansion.
    lc=lcm(*range(1,C))
    a0=ds['alpha']*2**(2*C-2)
    beta=ds['r1']/2**(C-1)
    bint=beta*2**(2*C-2)*lc
    assert a0.denominator==1 and bint.denominator==1
    assert q<=lc*abs(a0)
    nbar=[ns[j]-(ns[-1]*comb(C,j)) for j in range(C)]+[0]
    norm=sum(map(abs,nbar))
    assert abs(a0)<=2**(2*C-2)*norm
    assert abs(ds['alpha'])<=A+2*B
    coarse=2**(4*C-2)*lc*lc*(A+2*B)*shape
    assert upper<=coarse
    out=dict(m=m,L=L,c=c,H=h,R=d['integer_r'],alpha=d['integer_alpha'],
             reduced_fraction=str(d['rational']),q=q,rodrigues_factor=str(K),
             shifted_constant_kernel=dict(A=A,B=B,C=C,alpha=str(ds['alpha']),
                 boundary=str(beta),a0=int(a0),lcm=lc),
             exact_joint_quality_upper=str(upper),
             displayed_joint_quality_upper=float(upper),
             primitive_H_l1=sum(map(abs,h)),
             coarse_quality_upper=str(coarse),
             equivalence_verified=True)
    return out


if __name__=='__main__':
    result={'status':'passed','scope':'Explicit prescribed cancellation, exact arithmetic and finite joint certificates; no quality-zero theorem.',
            'samples':[data(m) for m in range(1,9)]}
    Path(__file__).with_name('explicit_periods.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([dict(m=s['m'],fraction=s['reduced_fraction'],q=s['q'],upper=s['displayed_joint_quality_upper']) for s in result['samples']],indent=2))
