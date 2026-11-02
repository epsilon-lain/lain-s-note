"""Small exact checks supporting the record-source review.

This does not rerun Zuuv, does not certify its complete output, and does not
claim the abstract guard counterexample is an execution trace of Zuuv.
"""
from fractions import Fraction as F
from pathlib import Path
from decimal import Decimal, localcontext
import random,json

ROOT=Path(__file__).resolve().parent

def convergent(terms):
    p0,p1,q0,q1=0,1,1,0
    for a in terms:p0,p1,q0,q1=p1,a*p1+p0,q1,a*q1+q0
    return p1,q1,p0,q0

def cf(v):
    a=[]
    while v.denominator!=1:
        t=v.numerator//v.denominator;a.append(t);v=1/(v-t)
    a.append(v.numerator)
    return a

def main():
    rng=random.Random(20261102)
    cases=0
    for length in range(2,22):
        for _ in range(50):
            terms=[rng.randrange(0,5)]+[rng.randrange(1,30) for _ in range(length-1)]
            A=rng.randrange(1,50);tail=F(A)+F(rng.randrange(1,100),101)
            p,q,pp,qq=convergent(terms)
            x=F(p*tail+pp,q*tail+qq)
            Q=q*q*abs(x-F(p,q))
            assert Q==1/(tail+F(qq,q))
            assert F(1,A+2)<Q<F(1,A)
            assert cf(x)[:len(terms)]==terms
            cases+=1
    # A positive correction tail alone is not a sufficient CF-prefix guard.
    x=F(7,5);terms=[1,1];p,q,pp,qq=convergent(terms)
    tail=(pp-x*qq)/(x*q-p)
    assert tail==F(2,3)>0 and cf(x)[:2]!=terms
    # Latest-record arithmetic is conditional on the reported coefficient.
    # Prior log checks are retained separately; no author files are needed here.
    A=845188380288
    with localcontext() as ctx:
        ctx.prec=50
        lower=str(Decimal(1)/Decimal(A+2));upper=str(Decimal(1)/Decimal(A))
    out={'exact_cf_identity_checks':cases,
         'same_sign_guard_counterexample':{'target':'7/5','hypothetical_prefix':[1,1],
           'actual_prefix':[1,2],'positive_correction_tail':'2/3',
           'scope':'Necessary-guard illustration, not a reproduced Zuuv failure.'},
         'record_quality_conditional_on_data':{'lower_exact':[1,A+2],'upper_exact':[1,A],
           'lower_display':lower,'upper_display':upper},
         'limits':'No full 653520000000-term rerun or independent input-pi certificate.'}
    (ROOT/'record_identity_checks.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
