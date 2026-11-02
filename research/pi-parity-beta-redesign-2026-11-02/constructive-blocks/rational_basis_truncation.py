"""Rational-basis truncation with exact rational/log(2) cancellation.

Candidate generation uses rational arithmetic only. Machin intervals are
applied afterwards to diagnose the fixed bounded controls.
"""
from fractions import Fraction as F
from math import comb,gcd,lcm,factorial
from functools import reduce,lru_cache
from pathlib import Path
import importlib.util,json

for ancestor in Path(__file__).resolve().parents:
    candidates = [ancestor/'constructive-review/audit_rodrigues.py',
                  ancestor/'pi-prescribed-cancellation-2026-11-02/constructive-review/audit_rodrigues.py']
    dependency = next((p for p in candidates if p.exists()), None)
    if dependency is not None:
        break
else:
    raise RuntimeError('The preceding prescribed-cancellation checkpoint is required.')
spec=importlib.util.spec_from_file_location('period_audit',dependency)
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)

def mul(a,b):
    z=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):z[i+j]+=x*y
    return z

def shifted(p):
    return [sum(p[k]*comb(k,j)*(-1)**(k-j) for k in range(j,len(p))) for j in range(len(p))]

@lru_cache(None)
def rational_moment(A,B,D):
    p=[0]*A+[(-1)**j*comb(B,j) for j in range(B+1)]
    coefficients=shifted(p)
    beta=F(0);log=F(0)
    for j,d in enumerate(coefficients):
        e=j-D+1
        if e==0:log+=d
        else:beta+=d*(F(2)**e-1)/e
    return beta,log

def pi_period(poly,c):
    # For degree<=2c-2 the odd part cannot produce log(2).
    assert len(poly)-1<=2*c-2
    even=poly[::2]
    beta,alpha=audit.period([F(x) for x in even],c)
    odd=poly[1::2]
    if odd:
        d=shifted(odd)
        for j,x in enumerate(d):
            e=j-c+1
            assert e!=0
            beta+=F(x,2)*(F(2)**e-1)/e
    return beta,alpha

def rational_pi(n=200):
    def at(q):
        x=sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(n)),F(0))
        y=x+F((-1)**n,(2*n+1)*q**(2*n+1))
        return min(x,y),max(x,y)
    a,b=at(5);c,d=at(239)
    return 16*a-4*d,16*b-4*c

def sample(m,K):
    c=4*m+1
    moments=[]
    for j in range(3):
        rat=log=F(0)
        for k in range(K+1):
            be,le=rational_moment(2*m+j+k,2*m+k,c+k)
            q=comb(c+k-1,k)
            rat+=q*be;log+=q*le
        moments.append((rat,log))
    a=[x[0] for x in moments];b=[x[1] for x in moments]
    h=[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
    scale=lcm(*(x.denominator for x in h))
    h=[int(x*scale) for x in h]
    content=reduce(gcd,map(abs,h))
    if not content:return dict(m=m,K=K,status='rank below two')
    h=[x//content for x in h]
    if next(x for x in h if x)!=abs(next(x for x in h if x)):h=[-x for x in h]
    assert sum(F(x)*y for x,y in zip(h,a))==0
    assert sum(F(x)*y for x,y in zip(h,b))==0
    base=[0]*(2*m)+[(-1)**j*comb(2*m,j) for j in range(2*m+1)]
    poly=mul(base,h)
    beta,alpha=pi_period(poly,c)
    if not alpha:return dict(m=m,K=K,H=[str(x) for x in h],status='zero pi period')
    candidate=-4*beta/alpha;q=candidate.denominator
    # z=x(1-x)/(1+x)<=7/40, an exact rational upper bound.
    delta=F(7,40);ratio=delta*F(c+K+1,K+2)
    if ratio>=1:return dict(m=m,K=K,status='tail ratio too large')
    eps=F(comb(c+K,K+1))*delta**(K+1)/(1-ratio)
    mass=F(factorial(2*m)**2,factorial(4*m+1))
    bound=4*q*q*sum(map(abs,h))*mass*eps/abs(alpha)
    lo,hi=rational_pi()
    actual_lo=max(F(0),lo-candidate,candidate-hi)*q*q
    actual_hi=max(abs(lo-candidate),abs(hi-candidate))*q*q
    assert actual_hi<bound
    return dict(m=m,K=K,c=c,status='passed',H=[str(x) for x in h],
                alpha=str(alpha),beta=str(beta),p=str(candidate.numerator),q=str(q),
                q_bits=q.bit_length(),quality_upper=str(bound),displayed_upper=float(bound),
                actual_quality_interval=[str(actual_lo),str(actual_hi)],displayed_actual_Q=float(actual_hi))

def run():
    rows=[]
    for m in range(1,4):
        c=4*m+1
        for K in [m,2*m,c,2*c]:
            rows.append(sample(m,K))
    out=dict(status='bounded controls complete',generation='Exact rational/log(2) crossproduct, no pi/target input',samples=rows,
             scope='Finite candidate generation and exact quality certificates; no all-scale nonzero-period or quality theorem inferred.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([{k:s.get(k) for k in ('m','K','status','q_bits','displayed_upper','displayed_actual_Q')} for s in rows],indent=2))

def rational_only_sample(m,K):
    c=4*m+K+3
    mu=[]
    for j in range(2):
        rat=log=F(0)
        for k in range(K+1):
            be,le=rational_moment(2*m+j+k,2*m+k,c+k)
            assert le==0
            rat+=comb(c+k-1,k)*be
        mu.append(rat)
    h=[mu[1],-mu[0]]
    scale=lcm(*(x.denominator for x in h))
    h=[int(x*scale) for x in h]
    content=reduce(gcd,map(abs,h));h=[x//content for x in h]
    assert sum(h[j]*mu[j] for j in range(2))==0
    base=[0]*(2*m)+[(-1)**j*comb(2*m,j) for j in range(2*m+1)]
    beta,alpha=pi_period(mul(base,h),c)
    if not alpha:return dict(m=m,K=K,c=c,status='zero pi period')
    candidate=-4*beta/alpha;q=candidate.denominator
    lo,hi=rational_pi()
    actual_hi=max(abs(lo-candidate),abs(hi-candidate))*q*q
    delta=F(7,40);ratio=delta*F(c+K+1,K+2)
    if ratio<1:
        eps=F(comb(c+K,K+1))*delta**(K+1)/(1-ratio)
        mass=F(factorial(2*m)**2,factorial(4*m+1))
        bound=4*q*q*sum(map(abs,h))*mass*eps/abs(alpha)
        assert actual_hi<bound
    else:bound=None
    return dict(m=m,K=K,c=c,status='passed',H=[str(x) for x in h],
                alpha=str(alpha),beta=str(beta),p=str(candidate.numerator),q=str(q),
                q_bits=q.bit_length(),quality_upper=str(bound),
                displayed_upper=float(bound) if bound is not None else None,
                actual_quality_upper=str(actual_hi),displayed_actual_Q=float(actual_hi),
                multiplier_scale_bits=scale.bit_length())

def rational_only_run():
    rows=[rational_only_sample(m,K) for m in range(1,4) for K in [2*m,4*m,8*m]]
    out=dict(status='bounded controls complete',generation='Higher pole makes all truncated moments rational',samples=rows,
             scope='Finite generation only. No convergence or quality-zero theorem inferred.')
    Path(__file__).with_name('rational_only_truncation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([{k:s.get(k) for k in ('m','K','c','status','q_bits','displayed_upper','displayed_actual_Q')} for s in rows],indent=2))

if __name__=='__main__':
    import sys
    rational_only_run() if '--rational-only' in sys.argv else run()
