from fractions import Fraction as F
from decimal import Decimal,localcontext
from pathlib import Path
from math import gcd
import json

HERE=Path(__file__).resolve().parent

def atan(b,n):
    s=sum((F((-1)**j,(2*j+1)*b**(2*j+1)) for j in range(n)),F(0))
    t=F((-1)**n,(2*n+1)*b**(2*n+1))
    return min(s,s+t),max(s,s+t)

def quality(p,q,L,U):
    l=q*q*L-p*q;u=q*q*U-p*q
    assert l*u>0
    return min(abs(l),abs(u)),max(abs(l),abs(u))

def pack(v):return [str(v.numerator),str(v.denominator)]
def dec(v):
    with localcontext() as c:
        c.prec=80
        return str(Decimal(v.numerator)/Decimal(v.denominator))

def main():
    a,b=atan(5,500);c,d=atan(239,180)
    L,U=16*a-4*d,16*b-4*c
    l,u=L,U;coeff=[]
    while True:
        x=l.numerator//l.denominator;y=u.numerator//u.denominator
        if x!=y:break
        coeff.append(x)
        l,u=1/(u-x),1/(l-x)
    bestlo,besthi=quality(355,113,L,U)
    p0,p1,q0,q1=0,1,1,0
    checked=[]
    target=None
    for i,a in enumerate(coeff[:-1]):
        p0,p1=p1,a*p1+p0;q0,q1=q1,a*q1+q0
        assert gcd(p1,q1)==1
        lo,hi=quality(p1,q1,L,U)
        row={'index_zero_based':i,'partial_quotient':a,'next_partial_quotient':coeff[i+1],
             'p':str(p1),'q':str(q1),'q_digits':len(str(q1)),
             'quality_lower':pack(lo),'quality_upper':pack(hi),'quality_display':dec((lo+hi)/2)}
        if i>3:
            if hi<bestlo:
                target=row
                break
            assert lo>besthi
        checked.append(row)
    assert target
    result={'pi_lower':pack(L),'pi_upper':pack(U),'pi_width':dec(U-L),
            'machin_terms':[500,180],'common_prefix_length':len(coeff),
            'baseline_quality':dec((bestlo+besthi)/2),
            'first_better_convergent':target,'previous_convergent':checked[-1],
            'all_intervening_convergents_checked':len(checked)-4,
            'continued_fraction_prefix_through_next':coeff[:target['index_zero_based']+2],
            'coverage':'Any primitive fraction with quality below baseline <1/2 must be a convergent by Legendre.'}
    assert F(22887386169930294,10**19) < F(*map(int,target['quality_lower']))
    assert F(*map(int,target['quality_upper'])) < F(22887386169930295,10**19)
    # Check the selected explicit integers again with a narrower independent enclosure.
    aa,bb=atan(5,520);cc,dd=atan(239,190)
    LL,UU=16*aa-4*dd,16*bb-4*cc
    pp,qq=int(target['p']),int(target['q'])
    ql,qu=quality(pp,qq,LL,UU)
    baseline_lower,_=quality(355,113,LL,UU)
    assert gcd(pp,qq)==1 and qu<baseline_lower
    assert F(22887386169930294,10**19)<ql<qu<F(22887386169930295,10**19)
    result['direct_recheck']={'machin_terms':[520,190], 'quality_lower':pack(ql),
                             'quality_upper':pack(qu), 'strictly_beats_355_113':True}
    (HERE/'first_improvement.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:result[k] for k in ['pi_width','common_prefix_length','baseline_quality','first_better_convergent','all_intervening_convergents_checked']},indent=2))

if __name__=='__main__':main()
