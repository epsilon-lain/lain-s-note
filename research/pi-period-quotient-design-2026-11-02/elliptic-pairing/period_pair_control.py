"""Shared denominator approximants to h and theta*h at the Ramanujan point.

Construction uses integer source coefficients only; quality intervals come from
the positive Ramanujan series afterwards. Standard library arithmetic.
"""
from fractions import Fraction as F
from math import comb, gcd, lcm
from pathlib import Path
import json


def source(n):
    v=[1]
    for k in range(n):
        a=4*(4*k+1)**2*v[-1]
        assert a%((k+1)**2)==0
        v.append(a//((k+1)**2))
    s=[(42*k+5)*comb(2*k,k)**3 for k in range(n+1)]
    assert all(sum(v[j]*v[k-j] for j in range(k+1))==comb(2*k,k)**3 for k in range(n+1))
    return v,s


def solve(a,b):
    n=len(b)
    mat=[[F(x) for x in row]+[F(y)] for row,y in zip(a,b)]
    for j in range(n):
        k=next(k for k in range(j,n) if mat[k][j])
        mat[j],mat[k]=mat[k],mat[j]
        p=mat[j][j];mat[j]=[x/p for x in mat[j]]
        for k in range(n):
            if k!=j:
                p=mat[k][j]
                mat[k]=[x-p*y for x,y in zip(mat[k],mat[j])]
    return [mat[j][-1] for j in range(n)]


def convolution(q,a,n):
    return sum((q[j]*a[n-j] for j in range(min(n,len(q)-1)+1)),F(0))


def evaluate(p,x=F(1,4096)):
    z=F(0)
    for a in reversed(p):z=z*x+a
    return z


def mul(a,b):
    z=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):z[i+j]+=x*y
    return z


def compatibility_defect(polys,N):
    q,p0,p1=polys
    one=mul(q,[i*x for i,x in enumerate(p0)])
    two=mul(p0,[i*x for i,x in enumerate(q)])
    three=mul(p1,q)
    d=[a-b-c for a,b,c in zip(one,two,three)]
    assert not any(d[:N+1])
    assert any(d[N+1:])
    return d[N+1:]


def normalize(*polys):
    den=lcm(*(x.denominator for p in polys for x in p))
    ints=[[int(x*den) for x in p] for p in polys]
    content=0
    for p in ints:
        for x in p:content=gcd(content,x)
    return [[x//content for x in p] for p in ints]


def pair(v,s):
    # Q degree2s; P0,P1 degree2s; both defects begin at3s+1.
    r=2*s;target=3*s
    deriv=[k*x for k,x in enumerate(v)]
    rows=[];rhs=[]
    for a in (v,deriv):
        for k in range(r+1,target+1):
            rows.append([a[k-j] if k>=j else 0 for j in range(1,r+1)])
            rhs.append(-a[k])
    q=[F(1)]+solve(rows,rhs)
    pp=[[convolution(q,a,k) for k in range(r+1)] for a in (v,deriv)]
    assert all(convolution(q,a,k)==0 for a in (v,deriv) for k in range(r+1,target+1))
    qi,p0i,p1i=normalize(q,*pp)
    qv,p0v,p1v=[evaluate(p) for p in (qi,p0i,p1i)]
    return F(16)*qv*qv/(p0v*(5*p0v+84*p1v)),(qi,p0i,p1i)


def pade(a,N):
    r=N//2;d=N-r
    q=[F(1)]+solve([[a[k-j] if k>=j else 0 for j in range(1,d+1)] for k in range(r+1,N+1)],[-a[k] for k in range(r+1,N+1)])
    p=[convolution(q,a,k) for k in range(r+1)]
    assert all(convolution(q,a,k)==0 for k in range(r+1,N+1))
    qi,pi=normalize(q,p)
    return F(16)*evaluate(qi)/evaluate(pi),(qi,pi)


def pi_interval(N=200):
    lower=sum((F((42*k+5)*comb(2*k,k)**3,4096**k) for k in range(N+1)),F(0))
    nxt=F((42*(N+1)+5)*comb(2*(N+1),N+1)**3,4096**(N+1))
    upper=lower+F(32,31)*nxt
    return F(16)/upper,F(16)/lower


def quality(r,lo,hi):
    if r<lo:return (lo-r)*r.denominator**2,(hi-r)*r.denominator**2
    if r>hi:return (r-hi)*r.denominator**2,(r-lo)*r.denominator**2
    raise AssertionError('Need narrower pi enclosure')


def pair_certificate(polys,N,v):
    x=F(1,4096)
    q,p0,p1=polys
    qv=evaluate(q,x)
    r0,r1=evaluate(p0,x)/qv,evaluate(p1,x)/qv
    tails=[]
    for deriv in (False,True):
        radius=F(0)
        for j,a in enumerate(q):
            ell=N-j
            first=v[ell+1]*x**(ell+1)
            if deriv:first*=ell+1
            radius+=abs(a)*x**j*first*(F(32,31) if deriv else F(64,63))
        tails.append(radius/abs(qv))
    e0,e1=tails
    center=5*r0*r0+84*r0*r1
    radius=10*abs(r0)*e0+5*e0*e0+84*(abs(r0)*e1+abs(r1)*e0+e0*e1)
    assert center>radius>0
    error=16*radius/(center*(center-radius))
    candidate=16/center
    return candidate.denominator**2*error,dict(period_radius=str(e0),derivative_radius=str(e1),ramanujan_radius=str(radius))


def main():
    lo,hi=pi_interval()
    rows=[]
    for s in range(1,5):
        N=3*s;v,a=source(N)
        rp,polys=pair(v,s);rs,spolys=pade(a,N)
        lp,up=quality(rp,lo,hi);ls,us=quality(rs,lo,hi)
        cert,components=pair_certificate(polys,N,source(N+1)[0])
        defect=compatibility_defect(polys,N)
        assert len(defect)<=s
        assert up<cert
        rows.append(dict(s=s,source_coefficients=N+1,
                         shared_pair=dict(p=str(rp.numerator),q=str(rp.denominator),
                             coefficient_height=max(abs(x) for p in polys for x in p),
                             Q_lower=str(lp),Q_upper=str(up),displayed_Q=float(up),
                             rational_pair_Q_certificate=str(cert),certificate_components=components,
                             derivative_compatibility_defect_after_t_power=defect,
                             defect_t_power=N+1,
                             polynomials=polys),
                         scalar_pade=dict(p=str(rs.numerator),q=str(rs.denominator),
                             coefficient_height=max(abs(x) for p in spolys for x in p),
                             Q_lower=str(ls),Q_upper=str(us),displayed_Q=float(us),
                             polynomials=spolys)))
    out=dict(status='passed',generation='Integer recurrence and exact common-denominator linear equations; pi enclosure only after construction.',controls=rows,
             scope='Four fixed source budgets; no inference of asymptotic improvement.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([dict(s=r['s'],source_coefficients=r['source_coefficients'],pair_q=r['shared_pair']['q'],pair_Q=r['shared_pair']['displayed_Q'],scalar_q=r['scalar_pade']['q'],scalar_Q=r['scalar_pade']['displayed_Q']) for r in rows],indent=2))


if __name__=='__main__':main()
