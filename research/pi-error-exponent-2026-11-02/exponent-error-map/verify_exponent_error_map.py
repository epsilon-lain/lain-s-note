"""Exact finite C/q^c certificates for the fixed reduced Legendre sequence.

Only the Python standard library is used. Rational c=k/l is checked by
raising positive quantities to l. No floating point logarithms or powers.
The all-n asymptotic/rate formulas are proved in the companion note.
"""
from fractions import Fraction as F
from math import comb, factorial, gcd
from pathlib import Path
import json


def valuation2(x):
    assert x>0
    return (x&-x).bit_length()-1


def f_value(n):
    m,r=divmod(n,4)
    return (2*m,2*m,2*m+2,2*m+3+valuation2(m+1))[r]


def arctan_interval(a,terms):
    s=sum((F((-1)**k,(2*k+1)*a**(2*k+1)) for k in range(terms)),F(0))
    t=F((-1)**terms,(2*terms+1)*a**(2*terms+1))
    return min(s,s+t),max(s,s+t)


def pi_interval():
    # The exact tangent identity (5+i)^4(239-i)=114244(1+i),
    # with angles between 0 and 1<pi/2, proves the selected Machin branch.
    z=(1,0)
    for _ in range(4):z=(5*z[0]-z[1],z[0]+5*z[1])
    z=(239*z[0]+z[1],239*z[1]-z[0])
    assert z==(114244,114244)
    assert 4*(F(1,5)-F(1,375))-F(1,239)>0 and F(4,5)<1
    a,b=arctan_interval(5,380),arctan_interval(239,120)
    lo,hi=16*a[0]-4*b[1],16*a[1]-4*b[0]
    assert 3<lo<hi<F(22,7) and hi-lo<F(1,10**520)
    return lo,hi


def e_interval():
    lo=sum((F(1,factorial(k)) for k in range(161)),F(0))
    hi=lo+F(1,factorial(161))/(1-F(1,162))
    assert 2<lo<hi<F(11,4) and hi-lo<F(1,10**280)
    return lo,hi


def classify(error_lo,error_hi,q,c,budget):
    k,l=c.numerator,c.denominator
    lo=error_lo**l*q**k
    hi=error_hi**l*q**k
    if hi<budget[0]**l:return 'passes'
    if lo>budget[1]**l:return 'fails'
    raise AssertionError('Interval too broad for this finite comparison')


def exponent_bracket(error_lo,error_hi,q):
    assert q>1
    scale=32
    def holds(k):return classify(error_lo,error_hi,q,F(k,scale),(F(1),F(1)))=='passes'
    lo,hi=0,128
    assert holds(lo) and not holds(hi)
    while hi-lo>1:
        mid=(lo+hi)//2
        if holds(mid):lo=mid
        else:hi=mid
    return {'strict_lower':str(F(lo,scale)),'strict_upper':str(F(hi,scale))}


def main():
    max_n=100
    a,b=[1,1],[0,1]
    for n in range(1,max_n):
        a.append((2*n+1)*a[n]+n*n*a[n-1])
        b.append((2*n+1)*b[n]+n*n*b[n-1])
    ell_values=[F(a[n],factorial(n)) for n in range(max_n+1)]
    lcms=[1]
    for n in range(1,max_n+1):
        lcms.append(lcms[-1]*n//gcd(lcms[-1],n))
    lcm_divisibility_checks=0
    for m in range(1,max_n//2+1):
        assert comb(2*m,m)%(lcms[2*m]//lcms[m])==0
        lcm_divisibility_checks+=1
        if 2*m+1<=max_n:
            assert ((2*m+1)*comb(2*m,m))%(lcms[2*m+1]//lcms[m])==0
            lcm_divisibility_checks+=1
    plo,phi=pi_interval()
    budgets={'1/10':(F(1,10),F(1,10)),
             '1':(F(1),F(1)),'10':(F(10),F(10)),
             'Euler e':e_interval()}
    cs=tuple(map(F,('1/2','3/4','1','4/3','3/2','2','5/2','3')))
    # Rigorous comparison 3/4<c_PNT: lambda>12/5, sqrt(2)<10/7,
    # and e<11/4 imply lambda^5>e^3*2^(3/2).
    assert F(7,5)**2<2<F(10,7)**2
    assert F(12,5)**5>F(11,4)**3*F(20,7)
    counts={label:{str(c):{'passes':0,'fails':0} for c in cs} for label in budgets}
    mapping_checks=0
    saturation_checks=0
    samples=[]
    for n in range(1,max_n+1):
        nf=factorial(n);g=gcd(a[n],4*b[n])
        p,q=4*b[n]//g,a[n]//g
        ell,R=F(a[n],nf),F(g,nf)
        assert gcd(p,q)==1 and q>0 and F(q)==ell/R
        numerator,denominator=R.numerator,R.denominator
        d=valuation2(numerator)-valuation2(denominator)
        E=numerator//(numerator&-numerator)
        D=denominator//(denominator&-denominator)
        assert E%2 and D%2 and gcd(E,D)==1
        assert R==F(2)**d*F(E,D)
        assert d==f_value(n)-(n-n.bit_count())-(n==3)
        assert 2*d+n>=0
        lcm_odd=lcms[n]//(lcms[n]&-lcms[n])
        assert lcm_odd%D==0
        assert R*R*lcms[n]**2*2**n>=1
        assert F(b[n],nf)==sum((ell_values[j-1]*ell_values[n-j]/j
                                for j in range(1,n+1)),F(0))
        saturation_checks+=1
        r=F(p,q)
        elo,ehi=(r-phi,r-plo) if n%2 else (plo-r,phi-r)
        assert 0<elo<ehi
        for label,budget in budgets.items():
            for c in cs:
                decision=classify(elo,ehi,q,c,budget)
                counts[label][str(c)][decision]+=1
                if budget[0]==budget[1]:
                    C=budget[0];k,l=c.numerator,c.denominator
                    # Original score and exact reduced-cost threshold agree.
                    assert elo**l*q**k/C**l==ell**k*(elo/C)**l/R**k
                    assert ehi**l*q**k/C**l==ell**k*(ehi/C)**l/R**k
                    # Exact net-odd version, retaining numerator and deficit.
                    assert (F(E,D)*F(2)**d)**k==R**k
                    mapping_checks+=1
        if n in (3,4,8,16,32,64,100):
            samples.append({'n':n,'p':str(p),'q':str(q),
                            'R_numerator':str(R.numerator),
                            'R_denominator':str(R.denominator),
                            'two_adic_cost_d':d,
                            'twice_beta':2*d+n,
                            'positive_odd_excess_E':str(E),
                            'odd_deficit_D':str(D),
                            'observed_exponent':exponent_bracket(elo,ehi,q),
                            'C_1_results':{str(c):classify(elo,ehi,q,c,budgets['1']) for c in cs},
                            'Euler_e_results':{str(c):classify(elo,ehi,q,c,budgets['Euler e']) for c in cs}})
    result={'status':'passed','arithmetic':'Python standard-library exact rational intervals',
            'maximum_n':max_n,
            'rational_exponent_error_comparisons':max_n*len(cs)*len(budgets),
            'exact_reduced_cost_threshold_comparisons':mapping_checks,
            'exact_two_part_gcd_saturation_checks':saturation_checks,
            'exact_convolution_and_odd_deficit_checks':max_n,
            'elementary_lcm_divisibility_checks':lcm_divisibility_checks,
            'three_quarters_strictly_below_PNT_sufficient_exponent':True,
            'pi_interval_width_less_than':'1e-520',
            'Euler_e_interval_width_less_than':'1e-280',
            'comparison_counts':counts,'samples':samples,
            'scope':'Finite exact classifications for the actual reduced sequence. The all-n exponent threshold, constants, and limsup formulas require the analytic proofs, not finite extrapolation. No global irrationality-measure or bad-approximability conclusion.'}
    Path(__file__).with_name('exponent_error_map_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('comparison_counts','samples')},indent=2))


if __name__=='__main__':
    main()
