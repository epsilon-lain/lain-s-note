"""Binomial filters of exact arctan Taylor sections, with shared denominator.

Generation uses no pi value. Strict pi intervals score candidates afterward.
The filter is the classical binomial/Euler transformation, not a novelty claim.
"""
from fractions import Fraction as F
from functools import cache
from math import comb, gcd, lcm
from pathlib import Path
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
import hashlib
import json
import sys

sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parent
LEGACY=ROOT.parent/'pi-small-norm-construction-2026-11-04'
sys.path.insert(0,str(LEGACY))
from verify import atan_taylor

FAMILIES={'gregory':[(4,1,1)],
          'two_angle':[(4,1,2),(4,1,3)],
          'machin':[(16,1,5),(-4,1,239)]}
NMAX=40

@cache
def section(s,t,n):
    c=F(s,t)
    return sum(((-1)**j*c**(2*j+1)/ (2*j+1) for j in range(n)),F(0))

@cache
def filtered(s,t,n,k):
    D=s*s+t*t
    return sum((F(comb(k,j)*s**(2*(k-j))*t**(2*j),D**k)*section(s,t,n+j)
                for j in range(k+1)),F(0))

@cache
def beta_moment(n,k):
    # Integral_0^1 x^(2n)(1-x^2)^k dx; exact positive beta product.
    answer=F(1,2*n+1)
    for j in range(1,k+1):
        answer*=F(2*j,2*n+2*j+1)
    return answer

def error_interval(s,t,n,k):
    c=F(s,t)
    base=c**(2*n+1)*(c*c/(1+c*c))**k*beta_moment(n,k)
    low,high=base/(1+c*c),base
    return (low,high) if n%2==0 else (-high,-low)

def raw_denominator(s,t,n,k):
    odd=1
    for j in range(n+k):odd=lcm(odd,2*j+1)
    return t**(2*n-1)*(s*s+t*t)**k*odd

def decimal_interval(lo,hi):
    with localcontext() as ctx:
        ctx.prec=18;ctx.rounding=ROUND_FLOOR
        a=str(Decimal(lo.numerator)/Decimal(lo.denominator))
        ctx.rounding=ROUND_CEILING
        b=str(Decimal(hi.numerator)/Decimal(hi.denominator))
    return [a,b]

def record_hash(record):
    return hashlib.sha256(json.dumps(record,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def run():
    generated=[]
    for family,angles in FAMILIES.items():
        for total in range(1,NMAX+1):
            for k in range(total):
                n=total-k
                r=sum((C*filtered(s,t,n,k) for C,s,t in angles),F(0))
                elo=ehi=F(0);D=1
                for C,s,t in angles:
                    al,au=error_interval(s,t,n,k)
                    if C<0:al,au=au,al
                    elo+=C*al;ehi+=C*au
                    d=raw_denominator(s,t,n,k)
                    assert d%filtered(s,t,n,k).denominator==0
                    D=lcm(D,d)
                assert D%r.denominator==0
                generated.append((family,total,n,k,r,elo,ehi,D))
    bits=2*max(r.denominator.bit_length() for _,_,_,_,r,_,_,_ in generated)+256
    l5,u5,_=atan_taylor(5,bits);l239,u239,_=atan_taylor(239,bits)
    pl,pu=16*l5-4*u239,16*u5-4*l239
    rows=[];full=[]
    for family,total,n,k,r,elo,ehi,D in generated:
        actual_lo,actual_hi=pl-r,pu-r
        assert elo<actual_lo<actual_hi<ehi
        q=r.denominator
        if actual_lo>0:ql,qh=q*q*actual_lo,q*q*actual_hi
        else:
            assert actual_hi<0
            ql,qh=-q*q*actual_hi,-q*q*actual_lo
        f={'family':family,'total_order':total,'n':n,'k':k,
           'p':str(r.numerator),'q':str(q),'raw_denominator':str(D),
           'raw_gcd':str(D//q),'error_interval':list(map(str,(elo,ehi))),
           'quality_interval':list(map(str,(ql,qh)))}
        full.append(f)
        rows.append({'family':family,'total_order':total,'n':n,'k':k,
                     'q_bits':q.bit_length(),'raw_denominator_bits':D.bit_length(),
                     'raw_gcd_bits':(D//q).bit_length(),
                     'actual_quality_decimal':decimal_interval(ql,qh),
                     'record_sha256':record_hash(f)})
    summaries={};selected={}
    for family in FAMILIES:
        own=[i for i,r in enumerate(rows) if r['family']==family]
        best=min(own,key=lambda i:F(full[i]['quality_interval'][1]))
        for i in own:
            if (full[i]['p'],full[i]['q'])!=(full[best]['p'],full[best]['q']):
                assert F(full[best]['quality_interval'][1])<F(full[i]['quality_interval'][0])
        improvements=[];wins=[]
        for total in range(1,NMAX+1):
            ids=[i for i in own if rows[i]['total_order']==total]
            base=next(i for i in ids if rows[i]['k']==0)
            win=min(ids,key=lambda i:F(full[i]['quality_interval'][1]))
            for i in ids:
                if (full[i]['p'],full[i]['q'])!=(full[win]['p'],full[win]['q']):
                    assert F(full[win]['quality_interval'][1])<F(full[i]['quality_interval'][0])
            if F(full[win]['quality_interval'][1])<F(full[base]['quality_interval'][0]):
                improvements.append(total)
            wins.append({'total_order':total,'winner_k':rows[win]['k'],
                         'winner_record_sha256':rows[win]['record_sha256'],
                         'baseline_record_sha256':rows[base]['record_sha256']})
            if total in (1,3,40):
                selected[str((family,total,'winner'))]=full[win]
                selected[str((family,total,'baseline'))]=full[base]
        selected[str((family,'best'))]=full[best]
        summaries[family]={'candidate_count':len(own),'best_parameters':{x:rows[best][x] for x in ('total_order','n','k')},
                           'best_fraction':full[best]['p']+'/'+full[best]['q'],
                           'best_actual_quality_decimal':rows[best]['actual_quality_decimal'],
                           'quality_improved_total_orders':improvements,'winners':wins}
        win=next(i for i in own if rows[i]['total_order']==NMAX and
                 rows[i]['record_sha256']==wins[-1]['winner_record_sha256'])
        base=next(i for i in own if rows[i]['total_order']==NMAX and rows[i]['k']==0)
        wl,wu=map(F,full[win]['quality_interval']);bl,bu=map(F,full[base]['quality_interval'])
        summaries[family]['matched_order_40']={'winner_n':rows[win]['n'],'winner_k':rows[win]['k'],
            'winner_q_bits':rows[win]['q_bits'],'baseline_q_bits':rows[base]['q_bits'],
            'quality_ratio_decimal':decimal_interval(wl/bu,wu/bl),
            'winner_quality_decimal':rows[win]['actual_quality_decimal']}
    family_names=list(FAMILIES)
    compact=[[family_names.index(r['family']),r['total_order'],r['n'],r['k'],r['record_sha256']] for r in rows]
    out={'schema':2,'construction':'Binomial/Euler arctan Taylor filters with shared rational denominator.',
         'parameters':{'nmax':NMAX,'n_min':1,'families':FAMILIES,'pi_score_bits':bits},
         'family_names':family_names,'candidate_columns':['family_id','total_order','n','k','record_sha256'],
         'candidates':compact,'selected_witnesses':selected,'summary':summaries,
         'scope':'Finite classical-filter diagnostic. Better than matched Taylor section is not a new pi record or infinite q-normalized success.'}
    (ROOT/'fresh_probe_results.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(json.dumps({f:{k:v for k,v in s.items() if k!='winners'} for f,s in summaries.items()},indent=2))
    return out

if __name__=='__main__':run()
