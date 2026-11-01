"""Finite exact checks for RESULT.md; the universal results use its proofs."""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json


def union_length(intervals):
    intervals = sorted(intervals)
    if not intervals:
        return F(0)
    lo, hi = intervals[0]
    length = F(0)
    for left, right in intervals[1:]:
        if left <= hi:
            hi = max(hi, right)
        else:
            length += hi-lo
            lo, hi = left, right
    return length+hi-lo


def success_intervals(A, r, B):
    intervals = []
    for b in range(1, B+1):
        for a in range(A*b, (A+1)*b+1):
            if gcd(a, b) != 1:
                continue
            center = F(a, b)
            radius = 1/(b*(a+b*r))
            lo, hi = max(F(A), center-radius), min(F(A+1), center+radius)
            if lo < hi:
                intervals.append((lo, hi))
    return intervals


def main():
    # A unimodular pair from pi; algebraic identity checks use exact rational
    # placeholder alpha, so they require no floating approximation of pi.
    P,Q,R,S = 355,113,333,106
    assert abs(P*S-Q*R) == 1
    r = F(S,Q)
    identity_count = 0
    for alpha in (F(7,3),F(29263,100),F(355,113)):
        x = (P*alpha+R)/(Q*alpha+S)
        for a in range(1,18):
            for b in range(1,12):
                p,q = a*P+b*R,a*Q+b*S
                direct = q*abs(q*x-p)
                predicted = (a+b*r)*abs(a-b*alpha)/(alpha+r)
                assert direct == predicted
                assert gcd(p,q) == gcd(a,b)
                if gcd(a,b)==1:
                    assert (direct < Q*abs(Q*x-P)) == ((a+b*r)*abs(a-b*alpha)<1)
                identity_count += 1

    cases=[]
    for A in (2,3,7,292,20776):
        for r in (F(1,2),F(106,113)):
            for B in (1,2,10,50):
                intervals=success_intervals(A,r,B)
                length=union_length(intervals)
                harmonic=sum((F(1,b) for b in range(1,B+1)),F(0))
                bound=2*harmonic/A
                assert 0<=length<=1 and length<=bound
                if bound<1:
                    assert length<1
                cases.append({'A':A,'r':str(r),'B':B,
                              'intervals':len(intervals),
                              'length_display_only':float(length),
                              'bound_exact':str(bound),
                              'inequality_checked_exactly':True})
    result={'status':'passed','identity_checks':identity_count,
            'coverage_checks':cases,
            'scope':'Finite exact sanity checks. Universal coverage and no-amplification statements are proved in RESULT.md. No probabilistic claim about pi.'}
    out=Path(__file__).with_name('verification.json')
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'passed','identity_checks':identity_count,
                      'coverage_checks':len(cases)}))


if __name__=='__main__':
    main()
