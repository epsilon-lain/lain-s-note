#!/usr/bin/env python3
"""Certify geometric integer-base functions by rectangle integration and layer waves.

Standard library only. No continued fractions or known good denominators are input.
All numerical certification uses exact rational intervals, not displayed floats.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
BITS=768
DEPTH=24
BASES=range(2,193)


def atan_bounds(q, n=180):
    s=sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(n)),F(0))
    t=F((-1)**n,(2*n+1)*q**(2*n+1))
    return min(s,s+t),max(s,s+t)


def dyadic_enclosure(lo,hi):
    den=2**BITS
    a=(lo*den).__floor__()
    b=(hi*den).__floor__()
    assert a==b
    return F(a,den),F(a+1,den)


def constants():
    a,b=atan_bounds(5)
    c,d=atan_bounds(239)
    plo,phi=16*a-4*d,16*b-4*c
    pi=dyadic_enclosure(plo-3,phi-3)
    den=2**BITS
    z=isqrt(2*den*den)
    sqrt2=(F(z-den,den),F(z+1-den,den))
    # Classical Liouville L=sum 10^(-k!), k>=1. First five terms.
    exps=[1,2,6,24,120]
    s=sum((F(1,10**e) for e in exps),F(0))
    liouville=dyadic_enclosure(s+F(1,10**720),s+F(2,10**720))
    return {'pi_minus_3':pi,'one_third':(F(1,3),F(1,3)),
            'sqrt2_minus_1':sqrt2,'liouville':liouville}


def digits(lo,hi,b,n):
    ds=[]; rs=[]
    for _ in range(n):
        lo,hi=b*lo,b*hi
        d=lo.__floor__()
        assert d==hi.__floor__(), ('uncertain digit',b,len(ds)+1)
        assert 0<=d<b
        lo,hi=lo-d,hi-d
        ds.append(d);rs.append((lo,hi))
    return ds,rs


def psi_interval(lo,hi):
    assert 0<=lo<=hi<=1
    p=lambda t:t*(1-t)
    return min(p(lo),p(hi)),F(1,4) if lo<=F(1,2)<=hi else max(p(lo),p(hi))


def geometric_interval(lo,hi,b,ds):
    # Integrate each actually filled rectangle. Keep the remaining active box.
    x=y=F(0);w=h=F(1);area=my=F(0)
    for i,d in enumerate(ds):
        if i%2==0:
            fw=w*d/b
            a=fw*h
            my+=a*(y+h/2)
            area+=a
            x+=fw;w/=b
        else:
            fh=h*d/b
            a=w*fh
            my+=a*(y+fh/2)
            area+=a
            y+=fh;h/=b
    tlo,thi=lo-area,hi-area
    assert 0<=tlo<=thi<=w*h
    prefix=b*(area-2*my)
    corners=[b*t*(1-2*v) for t in [tlo,thi] for v in [y,y+h]]
    return prefix+min(corners),prefix+max(corners)


def wave_interval(rs,b):
    lo=hi=F(0)
    for j in range(0,len(rs),2):
        a,c=psi_interval(*rs[j])
        d,e=psi_interval(*rs[j+1])
        scale=F(1,b**(3*(j//2)))
        lo+=scale*(a-e/b**2)
        hi+=scale*(c-d/b**2)
    return lo,hi+F(1,4*b**(3*(len(rs)//2)))


def main():
    const=constants()
    result={}
    for name,(alo,ahi) in const.items():
        rows=[]
        for b in BASES:
            ds,rs=digits(alo,ahi,b,DEPTH)
            gl,gh=geometric_interval(alo,ahi,b,ds)
            wl,wh=wave_interval(rs,b)
            # Independent integration and analytic layer-wave certificates agree.
            lo,hi=max(gl,wl,F(0)),min(gh,wh,F(1,4))
            assert lo<=hi
            rlo,rhi=rs[0]
            pl,ph=psi_interval(rlo,rhi)
            # Universal main term and its one-sided geometric correction.
            assert lo<=ph and hi>=pl-F(1,4*b*b)
            cl,ch=b*b*(pl-hi),b*b*(ph-lo)
            sl,sh=psi_interval(*rs[1])
            assert cl<=sh and ch>=sl-F(1,4*b)
            mid=(lo+hi)/2
            rows.append([b,float(mid),float((pl+ph)/2),
                         float(b*b*((pl+ph)/2-mid)),float((sl+sh)/2),
                         float((rlo+rhi)/2),float((rs[1][0]+rs[1][1])/2),float(hi-lo)])
        result[name]={'area_interval_exact':[[str(z.numerator),str(z.denominator)] for z in [alo,ahi]],
                      'rows':rows}
    data={'base_range_inclusive':[2,192],'digit_depth':DEPTH,'precision_bits':BITS,
          'definitions':{'J':'b*integral h_b(x)*(1-h_b(x)) dx',
                         'phase':'R_n=frac(b^n*alpha)','psi':'t*(1-t)',
                         'delta':'b^2*(psi(R_1)-J)'},
          'row_columns':['b','J_display','psi_R1_display','delta_display','psi_R2_display',
                         'R1_display','R2_display','certified_interval_width_display'],
          'checks':{'all_764_digit_prefixes_certified':True,
                    'all_764_independent_rectangle_and_wave_intervals_overlap':True,
                    'all_764_main_term_and_second_layer_bounds_consistent':True},
          'scope':'Integer bases only; displayed floats do not serve as certificates. No fitted coefficients or new approximation record.',
          'objects':result}
    for name,obj in result.items():
        best=min(obj['rows'],key=lambda r:r[1])
        print(name,'minimum displayed J in scanned bases:',best[:4])
    # A complete arithmetic progression of bases has a certified quadratic main curve.
    pl,ph=const['pi_minus_3']
    elo,ehi=F(1,7)-ph,F(1,7)-pl
    assert 0<elo<ehi and 192*ehi<1
    for row in result['pi_minus_3']['rows']:
        if row[0]%7==0:
            b=row[0]
            _,rs=digits(pl,ph,b,DEPTH)
            rlo,rhi=rs[0]
            assert rlo==1-b*ehi and rhi==1-b*elo
    data['checks']['all_27_multiples_of_7_have_certified_quadratic_phase']=True
    data['pi_mod_7_branch_definition']='epsilon=22/7-pi; R1=frac(r/7-b*epsilon) when b mod 7=r'
    # Explicit geometric conventions: same area can have different shapes.
    assert geometric_interval(F(1,4),F(1,4),2,[0,1])==(F(1,4),F(1,4))
    assert geometric_interval(F(1,4),F(1,4),4,[1])==(F(0),F(0))
    data['checks']['square_versus_vertical_strip_counterexample']=True
    (HERE/'results.json').write_text(json.dumps(data,ensure_ascii=False,separators=(',',':'))+'\n')
    print(json.dumps(data['checks'],indent=2))


if __name__=='__main__':
    main()
