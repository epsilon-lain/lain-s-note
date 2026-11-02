"""Exact finite counterexamples to oddpart(g_n) dividing oddpart(n!)."""
from fractions import Fraction as F
from math import factorial, gcd
from pathlib import Path
import json


def valuation(x,p):
    assert x>0
    value=0
    while x%p==0:
        x//=p
        value+=1
    return value


def odd_part(x):
    return x//(x&-x)


def main():
    a,b=[1,1],[0,1]
    for n in range(1,176):
        a.append((2*n+1)*a[n]+n*n*a[n-1])
        b.append((2*n+1)*b[n]+n*n*b[n-1])
    first=None
    for n in range(1,73):
        g=gcd(a[n],4*b[n])
        if odd_part(F(g,factorial(n)).numerator)>1:
            first=n
            break
    assert first==72
    cases=[]
    for n,p,expected_gv,expected_fv in ((72,17,5,4),(107,73,2,1),(176,23,8,7)):
        g=gcd(a[n],4*b[n]); nf=factorial(n); ratio=F(g,nf)
        gv,fv=valuation(g,p),valuation(nf,p)
        assert gv==expected_gv and fv==expected_fv
        assert odd_part(nf)%odd_part(g)!=0
        assert odd_part(ratio.numerator)==p
        cases.append({'n':n,'prime':p,'a_n':str(a[n]),'b_n':str(b[n]),
                      'gcd':str(g),'factorial':str(nf),
                      'gcd_valuation':gv,'factorial_valuation':fv,
                      'reduced_ratio_numerator':str(ratio.numerator),
                      'reduced_ratio_denominator':str(ratio.denominator),
                      'reduced_ratio_odd_numerator':p})
    result={'status':'passed','first_counterexample_in_n_1_through_72':first,
            'arithmetic':'Python standard-library exact integers/rationals',
            'cases':cases,
            'scope':'Finite exact counterexamples refute the universal odd-divisibility shortcut. They do not refute a weaker polynomial odd-excess bound or establish an infinite good approximation subsequence.'}
    Path(__file__).with_name('odd_excess_counterexamples.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'first_counterexample':first,
                      'valuations':[(c['n'],c['prime'],c['gcd_valuation'],c['factorial_valuation']) for c in cases]},indent=2))


if __name__=='__main__':
    main()
