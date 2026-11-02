#!/usr/bin/env python3
"""Exact zero-endpoint phase quadrature and prescribed-gcd construction.

Candidate generation and selection use no pi value. A rational Machin interval
is accessed only afterward to certify the actual finite outcomes.
"""
from fractions import Fraction as F
from math import gcd,lcm
from pathlib import Path
from decimal import Decimal,localcontext
import json

def trap(n,fun):
    return sum(((1 if j in (0,n) else 2)*fun(F(j,n)) for j in range(n+1)),F())/(2*n)

def base(t):return 4/(1+t*t)

def basis_derivative(t,k,m=1):
    return (k*t**(k-1)*(1-t)-t**k)/(1+t*t)**m-2*m*t**(k+1)*(1-t)/(1+t*t)**(m+1)

def egcd(a,b):
    if b==0:return abs(a),(1 if a>=0 else -1),0
    g,x,y=egcd(b,a%b)
    return g,y,x-(a//b)*y

def bezout_mod(m,bs):
    g=m;cs=[1]
    for b in bs:
        h,x,y=egcd(g,b)
        cs=[x*c for c in cs]+[y]
        g=h
    assert cs[0]*m+sum(x*y for x,y in zip(cs[1:],bs))==g
    return g,cs

def atan_bounds(k,n):
    a=sum((F((-1)**j,(2*j+1)*k**(2*j+1)) for j in range(n)),F())
    b=a+F((-1)**n,(2*n+1)*k**(2*n+1))
    return min(a,b),max(a,b)

def pi_bounds():
    # Angle range is known: 0<4 atan(1/5)-atan(1/239)<pi/2;
    # rational tangent addition yields 1, fixing pi/4 without a digit input.
    z=F(1,5)
    z=2*z/(1-z*z);z=2*z/(1-z*z)
    assert (z-F(1,239))/(1+z*F(1,239))==1
    a,b=atan_bounds(5,70),atan_bounds(239,20)
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]

def interval_error(r,lo,hi):
    if r<lo:return lo-r,hi-r
    if r>hi:return r-hi,r-lo
    return F(),max(r-lo,hi-r)

def dec(x):
    with localcontext() as c:
        c.prec=18
        return str(Decimal(x.numerator)/Decimal(x.denominator))

def scored(r,lo,hi):
    a,b=interval_error(r,lo,hi);q=r.denominator
    return {'fraction':str(r),'p':str(r.numerator),'q':str(q),
            'error_interval':[str(a),str(b)],'quality_interval':[str(a*q*q),str(b*q*q)],
            'quality_display':dec((a+b)*q*q/2)}

def main():
    n=2
    t0=trap(n,base)
    bs=[trap(n,lambda t,k=k:basis_derivative(t,k)) for k in (1,2)]
    assert t0==F(31,10) and bs==[F(9,200),F(-13,200)]
    M=lcm(t0.denominator,*(x.denominator for x in bs))
    A=int(M*t0);B=[int(M*x) for x in bs]
    assert (M,A,B)==(200,620,[9,-13])
    m=5;predictor=t0+F(1,6*n*n)
    candidates=[];eligible=[]
    for c1 in range(-4,5):
        for c2 in range(-4,5):
            raw=A+B[0]*c1+B[1]*c2;r=F(raw,M)
            row={'c':[c1,c2],'raw_numerator':raw,'raw_denominator':M,'gcd':gcd(raw,M),
                 'fraction':str(r),'q':r.denominator,'prescribed_divisor':m,
                 'proxy':str(r.denominator**2*abs(r-predictor))}
            assert r==trap(n,lambda t:base(t)+c1*basis_derivative(t,1)+c2*basis_derivative(t,2))
            candidates.append(row)
            if raw%m==0 and r.denominator>t0.denominator:eligible.append(row)
    eligible.sort(key=lambda r:F(r['proxy']))
    selected=eligible[0]
    assert selected['c']==[4,2] and selected['fraction']=='63/20'
    assert len(eligible)>1 and F(eligible[0]['proxy'])<F(eligible[1]['proxy'])
    selected_before_pi=dict(selected)
    # Predictive selection is robust on this entire coarse target interval.
    eps=F(1,1000)
    selected_proxy=F(selected['proxy']);selected_q=selected['q']
    robust_margin=min(F(row['proxy'])-selected_proxy-(row['q']**2+selected_q**2)*eps
                      for row in eligible[1:])
    assert robust_margin>0
    baseline_proxy=t0.denominator**2*abs(t0-predictor)
    uniform_gain=(baseline_proxy-t0.denominator**2*eps)/(selected_proxy+selected_q**2*eps)
    assert uniform_gain==F(61,56)>1

    # Generic prescribed-divisor certificates, independent of pi.
    congruence_checks=0;certificates=[]
    for nn in (2,3,4,5):
        tt=trap(nn,base)
        vv=[trap(nn,lambda t,k=k:basis_derivative(t,k)) for k in (1,2,3)]
        MM=lcm(tt.denominator,*(x.denominator for x in vv));AA=int(MM*tt)
        BB=[int(MM*x) for x in vv]
        for mm in range(1,101):
            if MM%mm:continue
            g,bc=bezout_mod(mm,BB)
            if AA%g==0:
                cc=[-(AA//g)*x for x in bc[1:]]
                assert (AA+sum(b*c for b,c in zip(BB,cc)))%mm==0
            else:
                assert AA%g!=0
            congruence_checks+=1
        certificates.append({'N':nn,'M':MM,'A':AA,'B':BB})

    # Natural polynomial phase: H=-t²(1-t)/3, no pi selection.
    em=[]
    for nn in (1,2,4,8,16):
        tr=trap(nn,base)
        rr=trap(nn,lambda t:base(t)-(2*t-3*t*t)/3)
        assert rr==tr+F(1,6*nn*nn)
        em.append({'N':nn,'base':tr,'corrected':rr})

    # Only now certify actual quality, after candidate generation/selection.
    lo,hi=pi_bounds()
    ca,cb=atan_bounds(5,3),atan_bounds(239,1)
    coarse_lo,coarse_hi=16*ca[0]-4*cb[1],16*ca[1]-4*cb[0]
    assert predictor-eps<coarse_lo<coarse_hi<predictor+eps
    assert hi-lo<F(1,10**95)
    assert selected==selected_before_pi
    baseline=scored(t0,lo,hi);chosen=scored(F(selected['fraction']),lo,hi)
    bq=list(map(F,baseline['quality_interval']));cq=list(map(F,chosen['quality_interval']))
    assert cq[1]<bq[0] and int(chosen['q'])>int(baseline['q'])
    selected['actual_quality']=chosen
    for row in em:
        tr,rr=row.pop('base'),row.pop('corrected')
        row['base']=scored(tr,lo,hi);row['corrected']=scored(rr,lo,hi)
        if row['N']>=2:
            assert F(row['corrected']['quality_interval'][1])<F(row['base']['quality_interval'][0])
    report={'status':'PASS','candidate_generation_uses_pi':False,'selection_uses_pi':False,
            'selection_predictor':str(predictor),'candidate_count':len(candidates),
            'eligible_count':len(eligible),'candidates':candidates,'selected':selected,
            'baseline':baseline,'quality_improvement_factor_lower':str(bq[0]/cq[1]),
            'robust_target_interval':[str(predictor-eps),str(predictor+eps)],
            'robust_selection_margin':str(robust_margin),
            'coarse_pi_certificate':[str(coarse_lo),str(coarse_hi)],
            'uniform_quality_improvement_factor':str(uniform_gain),
            'quality_improvement_factor_display':dec(bq[0]/cq[1]),
            'generic_congruence_checks':congruence_checks,'integer_models':certificates,
            'natural_polynomial_phase':em,
            'scope':'Finite exact joint gains and prescribed-gcd interface; no new pi record or infinite quality-zero sequence.'}
    Path(__file__).with_name('results.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','candidate_count':len(candidates),'selected':selected['fraction'],
          'base_quality':baseline['quality_display'],'selected_quality':chosen['quality_display'],
          'improvement_factor':report['quality_improvement_factor_display'],
          'generic_congruence_checks':congruence_checks},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
