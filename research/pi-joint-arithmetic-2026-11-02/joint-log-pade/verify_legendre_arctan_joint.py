"""Exact finite checks for LEGENDRE_ARCTAN_JOINT_COST.md.

Python standard library only. No floating point pi, numerical quadrature,
root finder, or guessed unreduced-denominator cost is used.
"""
from fractions import Fraction as F
from math import comb, factorial, gcd
from pathlib import Path
import json


def trim(p):
    while len(p)>1 and p[-1]==0:
        p.pop()
    return p


def add(p,q):
    r=[F(0)]*max(len(p),len(q))
    for i,x in enumerate(p): r[i]+=x
    for i,x in enumerate(q): r[i]+=x
    return trim(r)


def scale(p,s):
    return trim([s*x for x in p])


def mul(p,q):
    r=[F(0)]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        for j,y in enumerate(q): r[i+j]+=x*y
    return trim(r)


def integrate_symmetric(p):
    return sum(F(2,i+1)*x for i,x in enumerate(p) if i%2==0)


def polynomials(last):
    result=[[F(1)],[F(0),F(1)]]
    for n in range(1,last):
        result.append(scale(add(scale([F(0)]+result[-1],2*n+1),
                                scale(result[-2],-n)),F(1,n+1)))
    return result[:last+1]


def second_kind(p):
    result=[F(0)]*max(1,len(p)-1)
    for degree,x in enumerate(p):
        for r in range(degree):
            if r%2==0: result[degree-1-r]+=x/F(r+1)
    return trim(result)


def imag_phase_value(p,a,phase):
    total=F(0)
    for degree,x in enumerate(p):
        if not x: continue
        assert (degree-phase)%2==0
        sign=1 if ((degree-phase)//2)%2==0 else -1
        total+=x*a**degree*sign
    return total


def sequence(a,last):
    ell=[F(1),a]
    w=[F(0),F(1)]
    for n in range(1,last):
        ell.append(((2*n+1)*a*ell[n]+n*ell[n-1])/(n+1))
        w.append(((2*n+1)*a*w[n]+n*w[n-1])/(n+1))
    return ell[:last+1],w[:last+1]


def atan_interval(a,terms):
    """Strict alternating-series interval for atan(1/a), a>1."""
    assert a>1
    x=1/a
    power=x
    total=F(0)
    for k in range(terms):
        total+=((-1)**k)*power/(2*k+1)
        power*=x*x
    neighbor=total+((-1)**terms)*power/(2*terms+1)
    return min(total,neighbor),max(total,neighbor)


def scaled_interval(interval,c):
    assert c>0
    return c*interval[0],c*interval[1]


def abs_error_interval(target,r):
    low,high=target
    if r<low: return low-r,high-r
    if r>high: return r-high,r-low
    raise AssertionError('Target interval needs greater precision')


def power10(exponent):
    return F(10**exponent) if exponent>=0 else F(1,10**(-exponent))


def compact_interval(low,high,digits=10):
    """Returned integers times 10^exponent strictly enclose the input."""
    assert 0<low<=high
    exponent=len(str(low.numerator))-len(str(low.denominator))
    while low<power10(exponent): exponent-=1
    while low>=power10(exponent+1): exponent+=1
    unit=power10(exponent-digits)
    a=low/unit
    b=high/unit
    integer_low=a.numerator//a.denominator
    integer_high=b.numerator//b.denominator+1
    assert integer_low*unit<=low
    assert integer_high*unit>high
    return {'lower_integer':str(integer_low),
            'upper_integer':str(integer_high),
            'power_of_ten':exponent-digits,
            'interpretation':'lower_integer*10^power_of_ten <= value < upper_integer*10^power_of_ten'}


def polynomial_checks():
    ps=polynomials(12)
    orthogonal=norms=rodrigues=phases=0
    for n,p in enumerate(ps):
        # Independent factorial coefficient expansion of Rodrigues formula.
        expected=[F(0)]*(n+1)
        for k in range(n//2+1):
            expected[n-2*k]=F((-1)**k*factorial(2*n-2*k),
                                   2**n*factorial(k)*factorial(n-k)*factorial(n-2*k))
        assert p==trim(expected)
        rodrigues+=1
        for k in range(n):
            assert integrate_symmetric([F(0)]*k+p)==0
            orthogonal+=1
        assert integrate_symmetric(mul(p,p))==F(2,2*n+1)
        norms+=1
        wn=second_kind(p)
        for a in (F(1),F(5),F(239),F(3,2)):
            ell,w=sequence(a,n)
            assert imag_phase_value(p,a,n)==ell[n]
            if n: assert imag_phase_value(wn,a,n-1)==w[n]
            else: assert wn==[F(0)] and w[n]==0
            phases+=1
    return {'rodrigues_coefficients':rodrigues,'orthogonality':orthogonal,
            'squared_norms':norms,'endpoint_real_phase_and_W':phases}


def family(a,c,target,last):
    ell,w=sequence(a,last)
    target=scaled_interval(target,c)
    rows=[]
    identity_checks=error_checks=two_step_checks=0
    selected={0,1,2,3,5,10,20,50,100}
    for n in range(last+1):
        r=c*w[n]/ell[n]
        q=r.denominator
        assert gcd(r.numerator,q)==1 and q>0
        error=abs_error_interval(target,r)
        signed=(target[0]-r,target[1]-r)
        if n%2==0: assert signed[0]>0
        else: assert signed[1]<0
        lower=c*a/((a*a+1)*(2*n+1)*ell[n]**2)
        upper=c/(a*(2*n+1)*ell[n]**2)
        assert lower<error[0] and error[1]<upper
        quality=(q*q*error[0],q*q*error[1])
        row={'n':n,'fraction':r,'q':q,'error':error,'quality':quality}
        rows.append(row)
        error_checks+=1
        if n:
            assert ell[n]*w[n-1]-ell[n-1]*w[n]==F((-1)**n,n)
            previous=rows[n-1]
            difference=abs(r-previous['fraction'])
            assert difference==c/(n*ell[n]*ell[n-1])
            assert q*previous['q']>=n*ell[n]*ell[n-1]/c
            cross=r.numerator*previous['q']-previous['fraction'].numerator*q
            assert abs(cross)>=1
            theoretical=F(a*a*n*n,(a*a+1)**2*(4*n*n-1))
            assert quality[0]*previous['quality'][0]>theoretical
            identity_checks+=1
        if n>=2:
            previous2=rows[n-2]
            difference2=abs(r-previous2['fraction'])
            assert difference2==c*(2*n-1)*a/(n*(n-1)*ell[n]*ell[n-2])
            assert q*previous2['q']>=1/difference2
            theoretical2=F(n*n*(n-1)**2,
                           (2*n-1)**2*(a*a+1)**2*(2*n+1)*(2*n-3))
            assert quality[0]*previous2['quality'][0]>theoretical2
            two_step_checks+=1
    best=min(rows,key=lambda row:row['quality'][1])
    assert all(best['quality'][1]<row['quality'][0] for row in rows if row is not best)
    def export(row):
        return {'n':row['n'],'p':str(row['fraction'].numerator),'q':str(row['q']),
                'q_digits':len(str(row['q'])),
                'absolute_error':compact_interval(*row['error']),
                'normalized_quality':compact_interval(*row['quality'])}
    return {'a':str(a),'c':str(c),'last_n':last,
            'wronskian_reduced_denominator_and_adjacent_quality_checks':identity_checks,
            'same_parity_reduced_denominator_and_quality_checks':two_step_checks,
            'strict_analytic_error_and_alternating_sign_checks':error_checks,
            'unique_minimum_quality_in_checked_range':export(best),
            'selected_rows':[export(row) for row in rows if row['n'] in selected]},(ell,w)


def machin_combinations(seq5,seq239,pi_interval,last):
    rows=[]
    for n in range(1,last+1):
        r=16*seq5[1][n]/seq5[0][n]-4*seq239[1][n]/seq239[0][n]
        error=abs_error_interval(pi_interval,r)
        q=r.denominator
        rows.append({'n':n,'r':r,'q':q,'error':error,
                     'quality':(q*q*error[0],q*q*error[1])})
    best=min(rows,key=lambda row:row['quality'][1])
    assert all(best['quality'][1]<row['quality'][0] for row in rows if row is not best)
    def export(row):
        return {'n':row['n'],'p':str(row['r'].numerator),'q':str(row['q']),
                'q_digits':len(str(row['q'])),
                'absolute_error':compact_interval(*row['error']),
                'normalized_quality':compact_interval(*row['quality'])}
    return {'last_n':last,'strict_error_certificates':len(rows),
            'unique_minimum_quality_in_checked_range':export(best),
            'selected_rows':[export(row) for row in rows if row['n'] in (1,2,3,5,10,20,50,100)],
            'scope':'Same-index Machin combination; no adjacent-product theorem asserted for this combined sequence.'}


def main():
    atan5=atan_interval(F(5),500)
    atan239=atan_interval(F(239),200)
    # Exact tangent identities behind Machin's formula.
    tangent2=2*F(1,5)/(1-F(1,25))
    assert tangent2==F(5,12)
    tangent4=2*tangent2/(1-tangent2*tangent2)
    assert tangent4==F(120,119)
    assert (tangent4-F(1,239))/(1+tangent4*F(1,239))==1
    assert 4*F(5,26)-F(1,239)>0
    assert 4*F(1,5)<1
    pi_interval=(16*atan5[0]-4*atan239[1],16*atan5[1]-4*atan239[0])
    records=[]
    one,_=family(F(1),F(4),scaled_interval(pi_interval,F(1,4)),100)
    five,seq5=family(F(5),F(1),atan5,100)
    two39,seq239=family(F(239),F(1),atan239,100)
    rational,_=family(F(3,2),F(2,3),atan_interval(F(3,2),300),30)
    records.extend((one,five,two39,rational))
    result={'status':'passed','arithmetic':'exact fractions; Python standard library only',
            'scope':'All-n joint-quality result is analytic. Finite ranges do not exclude good subsequences or decide pi bad approximability.',
            'polynomial_checks':polynomial_checks(),
            'pi_interval':compact_interval(*pi_interval),
            'interval_provenance':{'atan_1_over_5_terms':500,'atan_1_over_239_terms':200,
                                   'machin_identity_tangent_checks':'passed'},
            'families':records,
            'same_index_machin':machin_combinations(seq5,seq239,pi_interval,100)}
    path=Path(__file__).with_name('legendre_arctan_joint_verification.json')
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'polynomial_checks':result['polynomial_checks'],
                      'families':[{'a':r['a'],'c':r['c'],'last_n':r['last_n'],
                                   'best':r['unique_minimum_quality_in_checked_range']} for r in records],
                      'machin_best':result['same_index_machin']['unique_minimum_quality_in_checked_range']},indent=2))


if __name__=='__main__':
    main()
