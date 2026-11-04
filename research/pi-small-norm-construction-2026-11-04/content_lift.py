"""Prescribe common content before testing the moving small-norm conics.

Adjacent primitive pi^2 sources are mixed with positive integer weights.  A
modular-kernel construction makes their common content explicit.  Rational-CF
convergents of the *known rational* kernel slope provide short positive weights.
The resulting primitive source height includes the entire weight cost.
"""
from fractions import Fraction as F
from math import gcd,isqrt
from pathlib import Path
import json,time
from conic_attempt import source,root_convergents,NMAX
from construct import display_interval
from verify import atan_taylor

ROOT=Path(__file__).resolve().parent

def bezout(a,b):
    r0,r1,s0,s1,t0,t1=a,b,1,0,0,1
    while r1:
        q=r0//r1
        r0,r1,s0,s1,t0,t1=r1,r0-q*r1,s1,s0-q*s1,t1,t0-q*t1
    return r0,s0,t0

def weights(A0,B0,A1,B1):
    G=abs(A0*B1-A1*B0)
    g,r,s=bezout(A1,B1);assert g==1
    k=(-(r*A0+s*B0))%G
    assert (A0+k*A1)%G==(B0+k*B1)%G==0
    choices={(1,0),(0,1)}
    if k:
        num,den=k,G;pm2,pm1,qm2,qm1=0,1,1,0
        while den:
            a,rem=divmod(num,den)
            for j in {1,a,max(1,a-1)}:
                if j>a:continue
                p,q=j*pm1+pm2,j*qm1+qm2
                x,y=q,k*q-G*p
                if x>0 and y>=0:
                    common=gcd(x,y);choices.add((x//common,y//common))
            p,q=a*pm1+pm2,a*qm1+qm2
            pm2,pm1,qm2,qm1=pm1,p,qm1,q
            num,den=den,rem
    candidates=[]
    for x,y in choices:
        aa,bb=x*A0+y*A1,x*B0+y*B1
        content=gcd(aa,bb)
        R=F(aa,bb)
        candidates.append((R.denominator,R.numerator,x,y,content))
    B,A,x,y,content=min(candidates)
    assert (x*A0+y*A1)//content==A and (x*B0+y*B1)//content==B
    return A,B,{'modulus':str(G),'kernel_slope_numerator':str(k),
                'weight0':str(x),'weight1':str(y),'content':str(content),
                'positive_weight_candidates':len(choices)}

def run():
    start=time.monotonic();u,v=source(NMAX+1);cases=[];hits=[];total=0
    for n in range(1,NMAX+1):
        r0,r1=6*v[n]/u[n],6*v[n+1]/u[n+1];delta=abs(r0-r1)
        A0,B0,A1,B1=r0.numerator,r0.denominator,r1.numerator,r1.denominator
        A,B,witness=weights(A0,B0,A1,B1)
        assert min(r0,r1)<=F(A,B)<=max(r0,r1)
        budgets={'fixed25':min(25,B),'square_root':isqrt(B)}
        limits={name:isqrt((H*delta.denominator)//(B*delta.numerator)) for name,H in budgets.items()}
        checked=0;best=None
        for p,q in root_convergents(A,B,max(limits.values())):
            if not 3*q<=p<=4*q:continue
            checked+=1;total+=1;norm=B*p*p-A*q*q
            best=abs(norm) if best is None else min(best,abs(norm))
            for name,H in budgets.items():
                if q<=limits[name] and abs(norm)<=H:
                    hits.append({'n':n,'budget':name,'p':str(p),'q':str(q),'norm':str(norm),
                                 'H':str(H),'q_limit':str(limits[name]),
                                 'quality_upper_from_conic':str((F(abs(norm),B)+q*q*delta)/6)})
        cases.append({'n':n,'A':str(A),'B':str(B),'source0_B':str(B0),'source1_B':str(B1),
                      'source_width':str(delta),'height_gain_over_source0':str(F(B0,B)),
                      'source0_B_squared_width':str(B0*B0*delta),
                      'new_B_squared_width':str(B*B*delta),
                      'content_witness':witness,'root_convergents_checked':checked,
                      'min_absolute_norm':str(best) if best is not None else None})
    l5,h5,_=atan_taylor(5,512);l239,h239,_=atan_taylor(239,512)
    lo,hi=16*l5-4*h239,16*h5-4*l239
    distinct=set()
    for hit in hits:
        p,q=int(hit['p']),int(hit['q']);r=F(p,q)
        if r<lo:ql,qh=q*(q*lo-p),q*(q*hi-p)
        else:
            assert r>hi
            ql,qh=q*(p-q*hi),q*(p-q*lo)
        assert qh<F(hit['quality_upper_from_conic'])
        hit['quality_interval']=[str(ql),str(qh)]
        hit['quality_decimal_enclosure']=display_interval(ql,qh)
        distinct.add((p,q))
    out={'construction':'Prescribed-content positive mixtures, followed by small-norm conics.',
         'n_max':NMAX,'cases':cases,'hits':hits,
         'summary':{'sources_preconditioned':len(cases),'root_convergents_checked':total,
                    'hits_with_budgets':len(hits),'distinct_fractions':[f'{p}/{q}' for p,q in sorted(distinct,key=lambda z:z[1])],
                    'last_source_original_B_digits':len(cases[-1]['source0_B']),
                    'last_source_new_B_digits':len(cases[-1]['B']),
                    'max_hit_n':max((h['n'] for h in hits),default=None),
                    'elapsed_seconds_display':round(time.monotonic()-start,6)},
         'scope':'All height gains include integer weight costs. No claim of globally minimal height, new pi record, or infinite small-norm solutions.'}
    (ROOT/'content_lift_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out['summary'],indent=2))

if __name__=='__main__':run()
