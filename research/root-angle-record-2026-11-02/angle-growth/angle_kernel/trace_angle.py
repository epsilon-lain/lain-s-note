"""Golden-parameter trace spinor: exact finite rational candidates.

Generation uses only the fixed integer kernel and node count. Pi enclosures
are used after generation for scoring. Python standard library only.
"""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
from decimal import Decimal, localcontext
import json

ROOT=Path(__file__).resolve().parent

def atan_bounds(b, n):
    s=sum((F((-1)**j,(2*j+1)*b**(2*j+1)) for j in range(n)),F(0))
    t=F((-1)**n,(2*n+1)*b**(2*n+1))
    return min(s,s+t),max(s,s+t)

def pi_bounds():
    a,b=atan_bounds(5,200);c,d=atan_bounds(239,70)
    return 16*a-4*d,16*b-4*c

def kernel(u):
    b=u*(1-u)
    return (1-6*b+12*b*b-4*b**3)/(1-6*b+10*b*b-8*b**3+4*b**4)

def flat6(u):
    b=u*(1-u)
    return 4*(1-4*b)*(1-b)**2/(1-6*b+10*b*b-8*b**3+4*b**4)

def base(u):
    return 4/(1+u*u)

def trapezoid(fun,n):
    return ((fun(F(0))+fun(F(1)))/2+sum((fun(F(j,n)) for j in range(1,n)),F(0)))/n

def score(v,lo,hi):
    q=v.denominator
    a,b=v-hi,v-lo
    assert a*b>0
    l,u=min(abs(a),abs(b)),max(abs(a),abs(b))
    return q*q*l,q*q*u

def dec(v,prec=30):
    with localcontext() as ctx:
        ctx.prec=prec
        return str(Decimal(v.numerator)/Decimal(v.denominator))

def pack(v):return [str(v.numerator),str(v.denominator)]

def main():
    lo,hi=pi_bounds()
    nodes=list(range(1,33))+[64,128]
    rows=[]
    for n in nodes:
        trace=trapezoid(kernel,n)
        six=trapezoid(flat6,n)
        old=trapezoid(base,n)
        vals={'trace_trapezoid':trace,'standard_trapezoid':old,
              'flat6_trapezoid':six,'standard_EM1':old+F(1,6*n*n)}
        if n%2==0:
            vals['trace_Simpson']=(4*trace-trapezoid(kernel,n//2))/3
            vals['standard_Simpson']=(4*old-trapezoid(base,n//2))/3
            vals['trace_Richardson4']=(16*trace-trapezoid(kernel,n//2))/15
            vals['flat6_Richardson6']=(64*six-trapezoid(flat6,n//2))/63
        for method,v in vals.items():
            l,u=score(v,lo,hi)
            rows.append({'n':n,'method':method,'fraction':pack(v),'reduced':gcd(v.numerator,v.denominator)==1,
                         'q_bits':v.denominator.bit_length(),'quality_lower':pack(l),'quality_upper':pack(u),
                         'quality_display':dec((l+u)/2),'error_display':dec(((l+u)/2)/v.denominator**2)})
    best=min(rows,key=lambda r:F(*map(int,r['quality_upper'])))
    trace_rows=[r for r in rows if r['method'].startswith('trace_')]
    besttrace=min(trace_rows,key=lambda r:F(*map(int,r['quality_upper'])))
    flat_rows=[r for r in rows if r['method'].startswith('flat6_')]
    bestflat=min(flat_rows,key=lambda r:F(*map(int,r['quality_upper'])))
    gain=[]
    flat_gain=[]
    flat_em1_gain=[]
    for n in nodes:
        a=next(r for r in rows if r['n']==n and r['method']=='trace_trapezoid')
        b=next(r for r in rows if r['n']==n and r['method']=='standard_EM1')
        if F(*map(int,a['quality_upper']))<F(*map(int,b['quality_lower'])):
            gain.append(n)
        c=next(r for r in rows if r['n']==n and r['method']=='flat6_trapezoid')
        if F(*map(int,c['quality_upper']))<F(*map(int,a['quality_lower'])):
            flat_gain.append(n)
        if F(*map(int,c['quality_upper']))<F(*map(int,b['quality_lower'])):
            flat_em1_gain.append(n)
    output={'description':'Fixed golden-parameter conjugate trace direction; no pi digits in candidate generation.',
            'node_counts':nodes,'pi_width':pack(hi-lo),'rows':rows,
            'best_all':best,'best_trace':besttrace,'trace_trapezoid_beats_standard_EM1_at_same_n':gain,
            'best_flat6':bestflat,
            'flat6_trapezoid_beats_standard_EM1_at_same_n':flat_em1_gain,
            'flat6_trapezoid_beats_trace_at_same_n':flat_gain,
            'scope':'Finite comparison only; no infinite quality-zero result.'}
    (ROOT/'trace_angle_results.json').write_text(json.dumps(output,indent=2))
    print(json.dumps({'candidate_count':len(rows),'best_trace':{k:besttrace[k] for k in ['n','method','fraction','quality_display']},
                      'best_all':{k:best[k] for k in ['n','method','fraction','quality_display']},
                      'trace_gain_vs_standard_EM1_nodes':gain,'flat6_gain_vs_trace_nodes':flat_gain,
                      'flat6_gain_vs_standard_EM1_nodes':flat_em1_gain},indent=2))

if __name__=='__main__':main()
