#!/usr/bin/env python3
"""Primitive-denominator design and positive-error certificates. Standard library."""
from fractions import Fraction as F
from math import gcd, lcm, factorial
from decimal import Decimal, localcontext
from pathlib import Path
import hashlib
import json
import sys

sys.set_int_max_str_digits(0)

LARGE_M = [3,5,7,9,11,13,25,27,49,81,125,243,625,729,2187,6561]

def prime_power(m):
    for p in range(3,m+1,2):
        if m%p == 0:
            x,r=m,0
            while x%p == 0:
                x//=p;r+=1
            return (p,r) if x==1 and all(p%d for d in range(2,int(p**.5)+1)) else None
    return None

def beta(a,b):
    assert a>0 and b>0
    return F(factorial(a-1)*factorial(b-1),factorial(a+b-1))

def atan_interval(k,n):
    a=sum((F((-1)**j,(2*j+1)*k**(2*j+1)) for j in range(n)),F())
    b=a+F((-1)**n,(2*n+1)*k**(2*n+1))
    return min(a,b),max(a,b)

def pi_interval():
    a,b=atan_interval(5,70),atan_interval(239,20)
    u=2*F(1,5)/(1-F(1,5)**2)
    v=2*u/(1-u*u)
    assert (v-F(1,239))/(1+v*F(1,239))==1
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]

def absolute_interval(r,bounds):
    lo,hi=bounds
    assert r<lo or r>hi
    return min(abs(r-lo),abs(r-hi)),max(abs(r-lo),abs(r-hi))

def display(r,digits=15):
    with localcontext() as c:
        c.prec=digits
        return str(Decimal(r.numerator)/Decimal(r.denominator))

def log10_display(r):
    assert r>0
    with localcontext() as c:
        c.prec=18
        return str((Decimal(r.numerator)/Decimal(r.denominator)).log10())

def make_sums(n):
    values=[F()]
    for j in range(n): values.append(values[-1]+F((-1)**j,2*j+1))
    return values

def main():
    sums=make_sums(max(LARGE_M)+1)
    towers=[]
    for m in LARGE_M:
        p,r=prime_power(m)
        old,new=4*sums[m],4*sums[m]-F(1,m)
        a,b,c=4*sums[m-1],old,4*sums[m+1]
        D=lcm(a.denominator,b.denominator,c.denominator)
        A,B,C=[x.numerator*(D//x.denominator) for x in (a,b,c)]
        delta=A-2*B+C
        N=A*C-B*B
        actual_g=gcd(abs(N),abs(D*delta))
        x=B-A
        ratio=(c-b)/(b-a)
        s,rr=ratio.denominator,ratio.numerator
        assert x%s==0
        h=x//s
        assert delta==h*(rr-s) and N==h*(A*(rr-s)-s*s*h)
        assert actual_g%abs(h)==0
        assert F(N,D*delta)==new
        expected=old.denominator//(3 if p==3 and r%2==0 else 1)
        assert new.denominator==expected and old.denominator%m==0
        lower=F(1,2*m*(m+1)*(2*m+1));upper=F(1,4*m**3)
        relative_upper=F(new.denominator**2,old.denominator**2)*upper/(F(1,m)-upper)
        towers.append({
            'm':m,'prime':p,'exponent':r,
            'old_numerator':str(old.numerator),'old_denominator':str(old.denominator),
            'new_numerator':str(new.numerator),'new_denominator':str(new.denominator),
            'old_denominator_bits':old.denominator.bit_length(),
            'new_denominator_bits':new.denominator.bit_length(),
            'denominator_reduction_factor':old.denominator//new.denominator,
            'positive_error_lower':str(lower),'positive_error_upper':str(upper),
            'quality_ratio_strict_upper':str(relative_upper),
            'quality_gain_strict_lower':str(1/relative_upper),
            'old_Q_log10_lower':log10_display(old.denominator**2*(F(1,m)-upper)),
            'new_Q_log10_lower':log10_display(new.denominator**2*lower),
            'predicted_common_factor_bits':abs(h).bit_length(),
            'actual_common_factor_bits':actual_g.bit_length(),
        })
    # Higher-order candidates are all fixed before any pi interval is computed.
    higher=[]
    for m in range(3,101,2):
        factor=prime_power(m)
        if not factor:continue
        base=4*sums[m]-F(1,m)
        value=base
        for k in range(1,min(8,2*m-1)+1):
            value+=F((-1)**(k-1),2**(k-1))*beta(2*m-k,2*k+1)
            lower=F(1,2**k)*beta(2*m-k,2*k+3)
            upper=2*lower
            higher.append({'m':m,'k':k,'fraction':str(value),'p':str(value.numerator),
                           'q':str(value.denominator),'base_q':str(base.denominator),
                           'error_lower':str(lower),'error_upper':str(upper),
                           'error_sign':(-1)**(k+1)})
    universal=[]
    for m in range(1,1001):
        old=4*sums[m]
        new=old-F((-1)**(m-1),m)
        U=F(1,6) if m==1 else beta(2*m-1,3)-beta(2*m-1,5)/2
        assert new.denominator<=m*old.denominator
        bound=F(new.denominator**2,old.denominator**2)*U/(F(1,m)-U)
        assert bound<F(1,4)
        if m&(m-1)==0:
            assert new.denominator==m*old.denominator
            universal.append({'m':m,'old_q_bits':old.denominator.bit_length(),
                              'new_q_bits':new.denominator.bit_length(),
                              'denominator_growth':m,'certified_ratio_upper':str(bound),
                              'ratio_upper_display':display(bound)})
    pi_bounds=pi_interval()
    assert pi_bounds[1]-pi_bounds[0]<F(1,10**95)
    for row in higher:
        m=row['m'];value=F(row['fraction'])
        base=4*sums[m]-F(1,m)
        elo,ehi=absolute_interval(value,pi_bounds)
        blo,bhi=absolute_interval(base,pi_bounds)
        assert F(row['error_lower'])<elo<=ehi<F(row['error_upper'])
        assert (1 if value>pi_bounds[1] else -1)==row['error_sign']
        q=value.denominator;qb=base.denominator
        lo,hi=F(q*q,qb*qb)*elo/bhi,F(q*q,qb*qb)*ehi/blo
        row['quality_ratio_interval']=[str(lo),str(hi)]
        row['quality_ratio_display']=display((lo+hi)/2)
        row['Q_display']=display(q*q*(elo+ehi)/2)
        row['classification']='improved' if hi<1 else ('worse' if lo>1 else 'unresolved')
    main_checks=0
    for row in towers:
        m=row['m'];old=F(int(row['old_numerator']),int(row['old_denominator']))
        new=F(int(row['new_numerator']),int(row['new_denominator']))
        elo,ehi=absolute_interval(new,pi_bounds)
        blo,bhi=absolute_interval(old,pi_bounds)
        assert F(row['positive_error_lower'])<elo<=ehi<F(row['positive_error_upper'])
        actual_upper=F(new.denominator**2,old.denominator**2)*ehi/blo
        assert actual_upper<F(row['quality_ratio_strict_upper'])
        main_checks+=1
    improving=[r for r in higher if r['classification']=='improved']
    best=min(higher,key=lambda r:F(r['quality_ratio_interval'][1]))
    out={'schema':1,'scope':'All-index fourfold net gain and infinite prime-power unbounded gains; no absolute pi record or quality-zero claim.',
         'candidate_generation_uses_pi':False,'prime_power_rows':towers,'high_order_rows':higher,
         'dyadic_sharpness_rows':universal,
         'pi_interval':[str(x) for x in pi_bounds],
         'summary':{'tower_rows':len(towers),'tower_exact_checks':main_checks,
                    'universal_quarter_bound_checks':1000,'dyadic_sharpness_cases':len(universal),
                    'higher_order_rows':len(higher),'higher_order_improved':len(improving),
                    'higher_order_worse':sum(r['classification']=='worse' for r in higher),
                    'higher_order_unresolved':sum(r['classification']=='unresolved' for r in higher),
                    'best_relative_row':{'m':best['m'],'k':best['k'],'ratio':best['quality_ratio_display'],
                                         'fraction':best['fraction'],'Q':best['Q_display']}}}
    path=Path(__file__).with_name('results.json')
    path.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(out['summary'],ensure_ascii=False,indent=2))

if __name__=='__main__':main()
