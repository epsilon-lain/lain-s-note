#!/usr/bin/env python3
"""Exact pi intervals and bounded unit-fraction experiments; standard library only."""
from fractions import Fraction as F
from decimal import Decimal, localcontext
from math import gcd
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def atan_bounds(inv, n=800):
    s = sum((F((-1)**k, (2*k+1)*inv**(2*k+1)) for k in range(n)), F(0))
    t = F((-1)**n, (2*n+1)*inv**(2*n+1))
    return min(s, s+t), max(s, s+t)


def decimal(x):
    with localcontext() as ctx:
        ctx.prec = 35
        return str(Decimal(x.numerator)/Decimal(x.denominator))


def ceil(x):
    return (x.numerator+x.denominator-1)//x.denominator


def continued_fraction_prefix(lo, hi, count):
    """Certify each coefficient directly from exact bounds; no saved CF data."""
    coefficients = []
    for _ in range(count):
        a = lo.numerator//lo.denominator
        assert a == hi.numerator//hi.denominator
        coefficients.append(a)
        assert lo > a
        lo, hi = 1/(hi-a), 1/(lo-a)
    return coefficients


def certificate(s, lo, hi):
    q = s.denominator
    return {'p': str(s.numerator), 'q': str(q), 'q_decimal_digits': len(str(q)),
            'error_lower_approx': decimal(lo-s),
            'error_upper_approx': decimal(hi-s),
            'quality_lower_approx': decimal(q*q*(lo-s)),
            'quality_upper_approx': decimal(q*q*(hi-s)),
            'quality_interval_width_approx': decimal(q*q*(hi-lo))}


def main():
    a, b = atan_bounds(5)
    c, d = atan_bounds(239)
    lo, hi = 16*a-4*d, 16*b-4*c
    s = F(3)
    greedy = []
    for step in range(1, 10):
        dl, dh = ceil(1/(hi-s)), ceil(1/(lo-s))
        assert dl == dh
        old = s
        s += F(1, dl)
        assert lo > s
        row = certificate(s, lo, hi)
        row.update({'step': step, 'unit_denominator': str(dl),
                    'raw_denominator_cancellation_gcd': str(gcd(old.numerator*dl+old.denominator, old.denominator*dl))})
        greedy.append(row)

    cf = continued_fraction_prefix(lo, hi, 800)
    p0,p1,q0,q1 = 0,1,1,0
    lower_convergents = {}
    for n,a in enumerate(cf[:800]):
        p0,p1 = p1,a*p1+p0
        q0,q1 = q1,a*q1+q0
        if n%2 == 0:
            lower_convergents[(p1,q1)] = n
    s = F(3)
    selected = []
    best_lo, best_hi = lo-3, hi-3
    record_rows = []
    for step in range(1, 5001):
        p,q = s.numerator,s.denominator
        residue = (-pow(p,-1,q))%q if q > 1 else 0
        def pick(t):
            return residue+q*max(0,ceil((t-residue)/q))
        kl,kh = pick(1/(q*(hi-s))),pick(1/(q*(lo-s)))
        assert kl == kh
        k = kl
        new = s+F(1,q*k)
        assert new.denominator == k
        assert new.numerator*q-p*k == 1
        assert new > s and new < lo
        s = new
        row = certificate(s,lo,hi)
        row.update({'step': step, 'convergent_index_if_matched': lower_convergents.get((s.numerator,s.denominator))})
        if step <= 4 or row['convergent_index_if_matched'] in {2, 4, 306, 430}:
            selected.append(row)
        quality_lo, quality_hi = (lo-s)*k*k, (hi-s)*k*k
        if quality_hi < best_lo:
            best_lo, best_hi = quality_lo, quality_hi
            record_rows.append(row)
        else:
            assert quality_lo > best_hi, 'Record comparison unresolved at current precision.'
        if len(str(k)) > 350:
            break
    out = {'pi_method': 'Exact alternating-series bounds, 800 terms each, Machin pi=16atan(1/5)-4atan(1/239)',
           'pi_interval_width_approx': decimal(hi-lo),
           'greedy': greedy,
           'max_cancellation_positive_walk': {
               'formula': 'p/q -> p/q + 1/(q*k), k = min positive integer >=1/[q*(pi-p/q)] and k == -p^(-1) mod q; new reduced denominator exactly k.',
               'steps_completed': step, 'last_q_decimal_digits': len(str(k)),
               'selected_rows': selected, 'quality_record_rows': record_rows},
           'cautions': ['Decimal fields are approximate displays; all checks and intervals are exact Fraction arithmetic.',
                        'The positive max-cancellation walk rediscovers one-sided continued fraction intermediates; it is not a new shortcut or record.']}
    (HERE/'results.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'greedy_first4':greedy[:4], 'walk_steps':step,
                      'walk_last_q_digits':len(str(k)), 'walk_records':record_rows},indent=2))


if __name__ == '__main__':
    main()
