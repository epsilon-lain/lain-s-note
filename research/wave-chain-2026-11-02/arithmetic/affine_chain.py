"""Exact affine-period arithmetic and finite critical-set optimization.

This verifier uses a rational pi enclosure only to rank already defined chain
candidates. It is not a pi-free new approximant generator or a novelty claim.
"""
from fractions import Fraction as F
from math import gcd,isqrt
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def egcd(a,b):
    oldr,r=a,b; olds,s=1,0; oldt,t=0,1
    while r:
        q=oldr//r
        oldr,r=r,oldr-q*r
        olds,s=s,olds-q*s
        oldt,t=t,oldt-q*t
    if oldr<0: return -oldr,-olds,-oldt
    return oldr,olds,oldt

def chain_data(P0,A0,P1,A1):
    h=gcd(P1,A1)
    assert h>0
    u,v=P1//h,A1//h
    g,s,t=egcd(u,v); assert g==1 and s*u+t*v==1
    r=s*P0+t*A0
    D=u*A0-v*P0
    assert D!=0,'collinear periods give a constant fraction'
    return {'h':h,'r':r,'D':D,'u':u,'v':v,'bezout':[s,t]}

def divisors(n):
    # Trial division is adequate for these fixtures. Factoring/height costs are
    # part of the general algorithm and are not claimed logarithmic here.
    n=abs(n)
    return sorted({x for d in range(1,isqrt(n)+1) if n%d==0 for x in (d,n//d)})

def congruence(h,r,d):
    g=gcd(h,d)
    if r%g: return None
    m=d//g
    a=0 if m==1 else (-r//g*pow(h//g,-1,m))%m
    return a,m

def interval_add(x,y): return x[0]+y[0],x[1]+y[1]
def interval_scale(x,a):
    z=(x[0]*a,x[1]*a)
    return min(z),max(z)
def interval_div(x,y):
    assert not y[0]<=0<=y[1], 'refine target enclosure to separate slope from zero'
    z=[a/b for a in x for b in y]
    return min(z),max(z)
def linear_interval(A,P,alpha):
    z=(A*alpha[0]-P,A*alpha[1]-P)
    return min(z),max(z)
def abs_interval(x):
    return (F(0) if x[0]<=0<=x[1] else min(abs(x[0]),abs(x[1])),
            max(abs(x[0]),abs(x[1])))
def floor(x): return x.numerator//x.denominator
def ceil(x): return -floor(-x)

def critical_intervals(P0,A0,P1,A1,alpha):
    r0,r1=linear_interval(A0,P0,alpha),linear_interval(A1,P1,alpha)
    assert not r1[0]<=0<=r1[1]
    roots=[interval_div(interval_scale(r0,-1),r1)]
    if A1:
        qroot=F(-A0,A1); roots.append((qroot,qroot))
        numerator=interval_scale(interval_add(interval_scale(r0,A1),interval_scale(r1,A0)),-1)
        roots.append(interval_div(numerator,interval_scale(r1,2*A1)))
    return roots

def candidate_set(P0,A0,P1,A1,alpha):
    data=chain_data(P0,A0,P1,A1)
    roots=critical_intervals(P0,A0,P1,A1,alpha)
    candidates=set();classes=[]
    for d in divisors(data['D']):
        cm=congruence(data['h'],data['r'],d)
        if cm is None:continue
        a,m=cm; local=set()
        for lo,hi in roots:
            jl=floor((lo-a)/m)-1
            jh=ceil((hi-a)/m)+1
            # Wide root brackets must be refined rather than silently guessed.
            assert jh-jl<1000,'refine the target interval'
            local.update(a+m*j for j in range(jl,jh+1))
        local={k for k in local if A0+k*A1!=0}
        candidates.update(local)
        classes.append({'d':d,'residue':a,'modulus':m,'candidate_k':sorted(local)})
    return sorted(candidates),data,classes,roots

def quality_interval(P0,A0,P1,A1,k,alpha):
    P,A=P0+k*P1,A0+k*A1
    if A==0:return None
    g=gcd(P,A)
    e=abs_interval(linear_interval(A,P,alpha))
    factor=F(abs(A),g*g)
    rat=F(P,A)
    return {'k':k,'p':rat.numerator,'q':rat.denominator,'g':g,
            'quality_lo':factor*e[0],'quality_hi':factor*e[1]}

def atan_interval(inv,n=60):
    s=sum((F((-1)**k,(2*k+1)*inv**(2*k+1)) for k in range(n)),F(0))
    nxt=F((-1)**n,(2*n+1)*inv**(2*n+1))
    return min(s,s+nxt),max(s,s+nxt)
def pi_interval():
    a,b=atan_interval(5);c,d=atan_interval(239)
    return 16*a-4*d,16*b-4*c

def encode(x):
    if isinstance(x,F): return {'exact':str(x),'display':float(x)}
    if isinstance(x,dict): return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    return x

def verify():
    gcd_checks=optimization_checks=congruence_checks=0
    # Deterministic broad sign, nonprimitive direction, and fixed-content cases.
    for P0,A0,P1,A1 in [(7,3,5,2),(-9,12,6,-4),(0,5,2,0),
                          (-924,-294,47,15),(24,36,-14,21),(1,2,0,3)]:
        data=chain_data(P0,A0,P1,A1)
        for k in range(-100,101):
            g=gcd(P0+k*P1,A0+k*A1)
            assert g==gcd(data['h']*k+data['r'],data['D'])
            gcd_checks+=1
            for d in divisors(data['D']):
                cm=congruence(data['h'],data['r'],d)
                condition=False if cm is None else (k-cm[0])%cm[1]==0
                assert condition==(g%d==0)
                congruence_checks+=1
        for alpha in [F(-11,7),F(7,5),F(22,7),F(355,113),F(13,3)]:
            if A1*alpha-P1==0:continue
            candidates,*_=candidate_set(P0,A0,P1,A1,(alpha,alpha))
            best=min(quality_interval(P0,A0,P1,A1,k,(alpha,alpha))['quality_hi'] for k in candidates)
            # The global candidate minimum must dominate each finite-window
            # brute check; if its selected k lies inside, equality must hold.
            brute=[quality_interval(P0,A0,P1,A1,k,(alpha,alpha)) for k in range(-200,201) if A0+A1*k]
            assert best<=min(r['quality_hi'] for r in brute)
            if any(abs(k)<=200 and quality_interval(P0,A0,P1,A1,k,(alpha,alpha))['quality_hi']==best for k in candidates):
                assert best==min(r['quality_hi'] for r in brute)
            optimization_checks+=1
    fixture=json.loads((HERE/'period_fixture.json').read_text())
    row0,row1=fixture['rows'];c=fixture['c'];M=fixture['base_multiplier']
    a0,a1=int(row0['alpha']),int(row1['alpha'])
    raw=(-M*sum(row0['R']),M*a0*2**(c-3),-sum(row1['R']),a1*2**(c-3))
    content=gcd(gcd(abs(raw[0]),abs(raw[1])),gcd(abs(raw[2]),abs(raw[3])))
    P0,A0,P1,A1=[x//content for x in raw]
    assert (P0,A0,P1,A1)==(313044,99645,-4948,-1575)
    alpha=pi_interval()
    candidates,data,classes,roots=candidate_set(P0,A0,P1,A1,alpha)
    rows=[quality_interval(P0,A0,P1,A1,k,alpha) for k in candidates]
    best=min(rows,key=lambda x:x['quality_hi'])
    assert (best['k'],best['p'],best['q'])==(3,355,113)
    assert all(best['quality_hi']<r['quality_lo'] for r in rows if r['k']!=best['k'])
    # Store an integer kernel representative to independently verify N and the
    # period identity; no division of kernel rows was used.
    R=[M*x+best['k']*y for x,y in zip(row0['R'],row1['R'])]
    coefficient=M*a0+best['k']*a1
    assert F(-sum(R),2**(c-3)*coefficient)==F(355,113)
    output={'status':'PASS','gcd_checks':gcd_checks,'congruence_checks':congruence_checks,
            'finite_bruteforce_optimization_controls':optimization_checks,
            'raw_periods':raw,'global_chain_content':content,
            'normalized_periods':[P0,A0,P1,A1],'chain_data':data,
            'admissible_divisor_classes':classes,'critical_intervals':roots,
            'number_distinct_candidates':len(candidates),'best':best,
            'best_integer_kernel':{'c':c,'alpha':coefficient,'R':R},
            'all_candidates':rows,
            'scope':'All integer k in this fixed affine chain are covered by the proof. This is a known 355/113 benchmark and target interval ranking, not a new approximation record or an infinite family quality-zero theorem.'}
    (HERE/'affine_chain_checks.json').write_text(json.dumps(encode(output),indent=2)+'\n')
    print(json.dumps(encode({k:v for k,v in output.items() if k in ['status','gcd_checks','congruence_checks','finite_bruteforce_optimization_controls','normalized_periods','number_distinct_candidates','best']}),indent=2))

if __name__=='__main__':verify()
