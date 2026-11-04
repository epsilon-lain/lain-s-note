"""Try small-norm conics from Apéry's pi-squared source and adjacent mixtures.

Inputs are the zeta(2) recurrence, n<=128, k=0..8, and two norm budgets.
No decimal pi or simple continued fraction of pi is used to generate candidates.
The small-norm condition, together with a source-error budget, would suffice
for Q->0 if successful conics existed with unbounded source denominator B.
This script makes a finite construction attempt, not that existence claim.
"""
from fractions import Fraction as F
from math import gcd,isqrt
from pathlib import Path
import hashlib,json,sys,time
from construct import display_interval
from verify import atan_taylor

sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parent
NMAX=128
MIXES=range(9)

def source(nmax):
    u=[F(1),F(3)];v=[F(0),F(5)]
    for n in range(1,nmax):
        c=11*n*n+11*n+3
        u.append((c*u[-1]+n*n*u[-2])/(n+1)**2)
        v.append((c*v[-1]+n*n*v[-2])/(n+1)**2)
    assert all(x.denominator==1 for x in u)
    return u,v

def root_convergents(A,B,qmax):
    """All ordinary convergents of sqrt(A/B), with positive q<=qmax."""
    D=A*B;s=isqrt(D)
    pm2,pm1,qm2,qm1=0,1,1,0
    if s*s==D:
        num,den=s,B
        while den:
            a,rem=divmod(num,den)
            p,q=a*pm1+pm2,a*qm1+qm2
            if q>qmax:return
            yield p,q
            pm2,pm1,qm2,qm1=pm1,p,qm1,q
            num,den=den,rem
        return
    P,Q=0,B
    while True:
        a=(s+P)//Q
        p,q=a*pm1+pm2,a*qm1+qm2
        if q>qmax:return
        yield p,q
        pm2,pm1,qm2,qm1=pm1,p,qm1,q
        P=a*Q-P
        assert (D-P*P)%Q==0
        Q=(D-P*P)//Q
        assert Q>0

def run():
    start=time.monotonic();u,v=source(NMAX+1)
    cases=[];hits=[];source_data=[];checked=0
    for n in range(1,NMAX+1):
        r0,r1=6*v[n]/u[n],6*v[n+1]/u[n+1]
        delta=abs(r0-r1)
        source_data.append({'n':n,'u':str(u[n]),'v':str(v[n]),
                            'pi_squared_endpoints':[str(min(r0,r1)),str(max(r0,r1))]})
        for k in MIXES:
            R=6*(v[n]+k*v[n+1])/(u[n]+k*u[n+1]);A,B=R.numerator,R.denominator
            budgets={'fixed25':min(25,B),'square_root':isqrt(B)}
            limits={name:isqrt((H*delta.denominator)//(B*delta.numerator)) for name,H in budgets.items()}
            total=0;best=None;case_hits=0
            for p,q in root_convergents(A,B,max(limits.values())):
                if not(3*q<=p<=4*q):continue
                total+=1;checked+=1
                norm=B*p*p-A*q*q
                if best is None or abs(norm)<best:best=abs(norm)
                for name,H in budgets.items():
                    if q>limits[name] or abs(norm)>H:continue
                    assert gcd(p,q)==1
                    tail=q*q*delta
                    assert B*tail<=H
                    # Exact identity: Q = |norm/B - q^2(pi^2-R)|/(pi+p/q).
                    bound=(F(abs(norm),B)+tail)/6
                    hits.append({'n':n,'mix':k,'budget':name,'A':str(A),'B':str(B),
                                 'p':str(p),'q':str(q),'q_digits':len(str(q)),
                                 'norm':str(norm),'H':str(H),'source_width':str(delta),
                                 'q_limit':str(limits[name]),'quality_upper_from_conic':str(bound)})
                    case_hits+=1
            cases.append({'n':n,'mix':k,'A_digits':len(str(A)),'B_digits':len(str(B)),
                          'q_limit_digits':{z:len(str(t)) for z,t in limits.items()},
                          'convergents_checked':total,'min_absolute_norm':str(best) if best is not None else None,
                          'hit_count_with_budgets':case_hits})
    # A separate pi interval only scores already generated candidates.
    maxq=max((int(h['q']).bit_length() for h in hits),default=1)
    l5,h5,_=atan_taylor(5,2*maxq+256)
    l239,h239,_=atan_taylor(239,2*maxq+256)
    lo,hi=16*l5-4*h239,16*h5-4*l239
    distinct={}
    for hit in hits:
        p,q=int(hit['p']),int(hit['q']);r=F(p,q)
        assert r<lo or r>hi
        qlo,qhi=(q*(q*lo-p),q*(q*hi-p)) if r<lo else (q*(p-q*hi),q*(p-q*lo))
        assert qhi<F(hit['quality_upper_from_conic'])
        hit['quality_interval']=[str(qlo),str(qhi)]
        hit['quality_decimal_enclosure']=display_interval(qlo,qhi)
        distinct[(p,q)]={'p':str(p),'q':str(q),'quality_decimal_enclosure':display_interval(qlo,qhi)}
    out={'construction':'Small-norm conics B*p^2-A*q^2=h from adjacent Apéry pi^2 sources.',
         'parameters':{'n_min':1,'n_max':NMAX,'mixes':list(MIXES),
                       'norm_budgets':{'fixed25':'min(25,B)','square_root':'floor(sqrt(B))'}},
         'sources':source_data,'cases':cases,'hits':hits,
         'summary':{'conics_tried':len(cases),'root_convergents_checked':checked,
                    'hits_with_budgets':len(hits),'distinct_fractions':list(distinct.values()),
                    'max_hit_n':max((h['n'] for h in hits),default=None),
                    'elapsed_seconds_display':round(time.monotonic()-start,6)},
         'scope':'New project experiment, classical Apéry/Pell ingredients. No new record or infinite successful family.'}
    (ROOT/'conic_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out['summary'],indent=2))
    return out

if __name__=='__main__':run()
