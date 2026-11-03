#!/usr/bin/env python3
"""Independent bounded exhaustive checks of the determinant-divisor one-step rule."""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json
from experiment import atan_bounds, ceil, decimal

HERE = Path(__file__).resolve().parent


def divisors_square(n):
    """Factor n, then double exponents; do not trial-divide n*n."""
    out, p = [1], 2
    while p*p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            out = [v*p**k for v in out for k in range(2*e+1)]
        p += 1
    if n > 1:
        out = [v*n**k for v in out for k in range(3)]
    return sorted(out)


def solve(x, lo, hi):
    p, q = x.numerator, x.denominator
    sign = 1 if x < lo else -1
    rl, rh = (lo-x,hi-x) if sign == 1 else (x-hi,x-lo)
    inverse = pow(p,-1,q) if q > 1 else 0
    rows = []
    assert x < lo or x > hi, 'Need a narrower target interval to certify side.'
    for det in divisors_square(q):
        residue = (-sign*det*inverse) % q
        def pick(r):
            t = det/(q*r)
            return residue + q*max(0, ceil((t-residue)/q))
        Q = pick(rl)
        assert Q == pick(rh) and Q > 0
        assert (p*Q+sign*det) % q == 0 and q*Q % det == 0
        d = q*Q//det
        new = F((p*Q+sign*det)//q,Q)
        assert new == x + F(sign,d)
        cl = new.denominator**2*(rl-F(1,d))
        ch = new.denominator**2*(rh-F(1,d))
        assert cl > 0
        rows.append((cl,ch,d,det,new))
    best = min(rows,key=lambda row:row[1])
    assert all(best[1] < row[0] or best[2] == row[2] for row in rows)
    return sign,rl,rh,best,rows


def exhaustive_check(x,lo,hi):
    sign,rl,rh,best,rows = solve(x,lo,hi)
    q = x.denominator
    dmin = ceil(1/rl)
    assert dmin == ceil(1/rh)
    cl,B,d_best,det,new = best
    # Since gcd(p*d+sign*q,q*d) divides q^2, reduced Q>=d/q.
    # Hence quality >= (rl*d^2-d)/q^2. Beyond M it already exceeds B.
    def excluded(d):
        return rl*d*d-d > B*q*q
    M = dmin
    while not excluded(M):
        M *= 2
    low,high = dmin,M
    while low+1 < high:
        mid = (low+high)//2
        if excluded(mid):
            high = mid
        else:
            low = mid
    M = high
    assert d_best < M
    rn,rd = rl.numerator,rl.denominator
    bn,bd = B.numerator,B.denominator
    count = 0
    for d in range(dmin,M):
        if d == d_best:
            continue
        rawp = x.numerator*d+sign*q
        rawq = q*d
        Q = rawq//gcd(rawp,rawq)
        # Exact comparison B < Q^2 * (rl-1/d), without Fraction reductions.
        assert bn*rd*d < Q*Q*(rn*d-rd)*bd, (x,d,d_best)
        count += 1
    return count,det,dmin,M


def main():
    a,b = atan_bounds(5,100)
    c,d = atan_bounds(239,100)
    lo,hi = 16*a-4*d,16*b-4*c
    samples = []
    for x in [F(3),F(25,8),F(22,7),F(333,106),F(355,113),F(103993,33102)]:
        sign,rl,rh,best,rows = solve(x,lo,hi)
        cl,ch,d,det,new = best
        samples.append({'start':str(x),'sign':sign,'branches':len(rows),'best_unit_denominator':str(d),
            'branch_determinant':str(det),'best_fraction':str(new),'quality_lower_approx':decimal(cl),
            'quality_upper_approx':decimal(ch)})
    seen = set()
    total,nonunitdet = 0,[]
    for q in range(1,41):
        center = int(lo*q)
        for p in range(max(1,center-2),center+4):
            x = F(p,q)
            if x in seen:
                continue
            seen.add(x)
            count,det,dmin,M = exhaustive_check(x,lo,hi)
            total += count
            if det != 1:
                nonunitdet.append({'start':str(x),'best_branch_determinant':det})
    out = {'samples':samples,
        'exhaustive_check':{'distinct_rationals':len(seen),'other_unit_denominators_tested_exactly':total,
                            'beyond_finite_range_excluded_by_bound':True,
                            'non_unit_best_branch_examples':nonunitdet[:20]},
        'proof_note':'Global one-candidate-per-determinant proof must handle later gcd cancellation by descending to the reduced determinant branch; monotonicity alone is insufficient.'}
    (HERE/'finite_update_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__ == '__main__':
    main()
