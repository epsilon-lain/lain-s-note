#!/usr/bin/env python3
"""Combine rational curve templates, geometric stability, and reduced costs.

Uses exact Fraction arithmetic. Requires the previous base-function compute.py,
but no external Python packages. The constructor-template compilation stage
does not read pi; pi bounds are used only for subsequent validation.
"""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import importlib.util
import hashlib
import json

HERE=Path(__file__).resolve().parent
OLD=next(p for p in [HERE.parent/'pi-base-function-2026-11-03',
                    HERE.parent/'base_function_research'] if (p/'compute.py').exists())
spec=importlib.util.spec_from_file_location('certified_geometry',OLD/'compute.py')
geom=importlib.util.module_from_spec(spec)
spec.loader.exec_module(geom)


def psi(t): return t*(1-t)


def rational_J(a,q,b):
    """Sum every geometric layer exactly using the finite residue orbit."""
    u=F(a,q);a,q=u.numerator,u.denominator
    assert 0<=u<1 and b>=2
    states=[];seen={};s=a
    while s not in seen:
        seen[s]=len(states);states.append(s)
        s=(b*b*s)%q
    pre=seen[s];period=len(states)-pre
    terms=[psi(F((b*s)%q,q))-psi(F((b*b*s)%q,q))/b**2 for s in states]
    z=F(1,b**3)
    head=sum((z**j*terms[j] for j in range(pre)),F(0))
    cycle=sum((z**j*terms[pre+j] for j in range(period)),F(0))
    out=head+z**pre*cycle/(1-z**period)
    assert 0<=out<=F(1,4)
    for j,s in enumerate(states):
        assert F(s,q).denominator==q//gcd(q,b**(2*j))
    return out,{'preperiod':pre,'period':period,'states':len(states)}


def J_interval(lo,hi,b):
    ds,rs=geom.digits(lo,hi,b,24)
    gl,gh=geom.geometric_interval(lo,hi,b,ds)
    wl,wh=geom.wave_interval(rs,b)
    low,high=max(gl,wl,F(0)),min(gh,wh,F(1,4))
    assert low<=high
    return low,high


def compile_constructor_template():
    # Uses the repository's saved rational construction result only.
    data=json.loads((HERE/'constructor_summary.json').read_text())
    row=next(r for r in data if r['m']==2 and not r['forced'])
    p,q=int(row['p']),int(row['q'])
    upper=F(row['strict_quality_upper'])
    assert F(3)<F(p,q)<4 and gcd(p,q)==1
    a=p-3*q
    bases=list(range(2,513))+[1024,2048,4096,8192,q*q]
    curve=[]
    for b in bases:
        v,orbit=rational_J(a,q,b)
        curve.append([b,float(v),float(b*upper/q**2),orbit['states']])
    return p,q,upper,bases,curve


def quality_row(alo,ahi,b,m):
    n=b**m
    # Integer part 3 stays exact even at bases 2 and 3.
    xl,xh=n*(3+alo),n*(3+ahi)
    P=(xl+F(1,2)).__floor__()
    assert P==(xh+F(1,2)).__floor__()
    if xl>P: tl,th=xl-P,xh-P
    else:
        assert xh<P
        tl,th=P-xh,P-xl
    assert 0<tl<th<F(1,2)
    k=b**(m-1)
    d=(k*alo).__floor__()
    assert d==(k*ahi).__floor__()
    ul,uh=k*alo-d,k*ahi-d
    jl,jh=J_interval(ul,uh,b)
    g=gcd(P,n)
    scale=F(n,g*g)
    zl,zh=scale*jl,scale*jh
    cl,ch=scale*tl,scale*th
    # Strong interval checks of the multiplicative theorem, not float checks.
    assert zh<cl
    assert ch<F(2*b,b-1)*zl
    return {'b':b,'m':m,'P':P,'label_denominator':n,'gcd':g,
            'p_reduced':P//g,'q_reduced':n//g,
            'Z_display':float((zl+zh)/2),'C_display':float((cl+ch)/2),
            'ratio_Z_over_C_display':float((zl+zh)/(cl+ch)),
            'multiplicative_bound_certified':True}


def main():
    # Stage 1: compile an entire curve from old rational constructor output.
    p,q,U,bases,curve=compile_constructor_template()
    # Stage 2: independently validate against certified pi.
    const=geom.constants();pl,ph=const['pi_minus_3']
    u=F(p-3*q,q)
    assert max(abs(u-pl),abs(u-ph))<U/q**2
    templates=[(22,7),(333,106),(p,q)]
    scan=[];exact_count=0
    for b in bases:
        jl,jh=J_interval(pl,ph,b)
        vals=[]
        for pp,qq in templates:
            v,orb=rational_J(pp-3*qq,qq,b)
            gl,gh=J_interval(F(pp-3*qq,qq),F(pp-3*qq,qq),b)
            assert gl<=v<=gh
            exact_count+=1
            delta=max(abs(F(pp,qq)-3-pl),abs(F(pp,qq)-3-ph))
            assert max(abs(jl-v),abs(jh-v))<=b*delta
            vals.append(float(v))
        scan.append([b,float((jl+jh)/2),*vals,float(jh-jl)])
    # Closed branch formulas cover infinite congruence classes; sample check here.
    for b in bases:
        v,_=rational_J(1,7,b)
        if b%7==0: assert v==0
        if b%7 in [1,6]: assert v==F(6*b*(b+1),49*(b*b+b+1))
    # A global shape coincidence cannot by itself recover the number.
    for b in range(2,65):
        assert rational_J(1,3,b)[0]==rational_J(2,3,b)[0]
    quality=[quality_row(pl,ph,b,m) for b in range(2,65) for m in range(1,5)]
    examples=[quality_row(pl,ph,b,m) for b,m in [(7,1),(14,1),(113,1),(226,1),(10,2),(12769,1)]]
    output={
        'scope':'Finite exact integration experiments plus proved general bridges; no new approximation record.',
        'source_commit':'ce9ee54ef02e430d29b84e1dbb1b95c5be23041b',
        'source_constructor_path':'research/period-content-design-2026-11-03/global_summary.json',
        'source_constructor_sha256':hashlib.sha256((HERE/'constructor_summary.json').read_bytes()).hexdigest(),
        'geometry_source_sha256':hashlib.sha256((OLD/'compute.py').read_bytes()).hexdigest(),
        'constructor':{'p':p,'q':q,'strict_quality_upper_exact':str(U),
                       'curve_uniform_error_upper_for_b_at_most_q_squared_exact':str(U),
                       'compilation_reads_pi':False,'old_integral_certificate_rebuilt_this_round':False,
                       'curve_columns':['b','J_template_display','strict_error_radius_display','state_count'],
                       'curve':curve},
        'scan_columns':['b','J_pi_display','J_1_over_7_display','J_15_over_106_display',
                        'J_16_over_113_display','J_pi_interval_width_display'],
        'scan':scan,'quality_examples':examples,'quality_scan':quality,
        'checks':{'template_closed_forms_in_independent_integral_intervals':exact_count,
                  'cross_carry_stability_checks':len(bases)*len(templates),
                  'multiplicative_quality_checks':len(quality)+len(examples),
                  'rational_denominator_state_identities':True,
                  'mod_7_closed_branches':True,'reflection_symmetry_63_bases':True}}
    (HERE/'results.json').write_text(json.dumps(output,ensure_ascii=False,separators=(',',':'))+'\n')
    print(json.dumps({'checks':output['checks'],'examples':examples},indent=2))


if __name__=='__main__': main()
