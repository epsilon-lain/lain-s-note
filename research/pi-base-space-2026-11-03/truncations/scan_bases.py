#!/usr/bin/env python3
"""Certified radix truncation experiment: pi, sqrt(2), decimal Liouville constant.

Requires gmpy2. All digit/quality comparisons use exact rational endpoints.
Logarithmic display and plot data use MPFR at 100 bits.
"""
import argparse
import json
import math
from pathlib import Path
import gmpy2 as g

HERE = Path(__file__).resolve().parent
CURVE_BASES = {2, 3, 10, 16, 20, 32, 40, 50, 64}


def target_bounds(bits):
    with g.context(precision=bits, round=g.RoundDown):
        a5l = g.atan(g.mpfr(g.mpq(1,5)))
        a239l = g.atan(g.mpfr(g.mpq(1,239)))
    with g.context(precision=bits, round=g.RoundUp):
        a5h = g.atan(g.mpfr(g.mpq(1,5)))
        a239h = g.atan(g.mpfr(g.mpq(1,239)))
    with g.context(precision=bits, round=g.RoundDown):
        pl = g.mpq(16*a5l-4*a239h)
    with g.context(precision=bits, round=g.RoundUp):
        ph = g.mpq(16*a5h-4*a239l)
    scale = g.mpz(1)<<bits
    sr = g.isqrt(2*scale*scale)
    sl, sh = g.mpq(sr,scale), g.mpq(sr+1,scale)
    partial = sum((g.mpq(1,g.mpz(10)**math.factorial(k)) for k in range(1,7)),g.mpq(0))
    next_term = g.mpq(1,g.mpz(10)**5040)
    # Later exponents are 40320,362880,...; their sum is <10^-5040.
    ll, lh = partial+next_term, partial+2*next_term
    return {'pi':(pl,ph), 'sqrt2':(sl,sh), 'liouville_decimal':(ll,lh)}


def rational_hex(x):
    return [x.numerator.digits(16),x.denominator.digits(16)]


def log10(x):
    assert x>0
    return float(g.log10(g.mpfr(x)))


def floor_exact(x):
    # g.floor(mpq) can round via MPFR: use exact integer division instead.
    return x.numerator//x.denominator


def row(m, p, raw_q, divisor, cl, ch, z, censored):
    q=raw_q//divisor
    return {'m':m,'p_hex':(p//divisor).digits(16),'q_hex':q.digits(16),
            'reduced_q_decimal_digits':len(q.digits(10)),
            'log10_reduced_q':log10(q),'log10_raw_denominator':log10(raw_q),
            'log10_gcd_cancellation':log10(divisor),
            'log10_C_lower':log10(cl),'log10_C_upper':log10(ch),
            'following_zeros_observed':z,'zero_run_right_censored':censored}


def scan(target,b,digits,lookahead):
    lo,hi=target
    length=digits+lookahead
    raw=g.mpz(b)**length
    a=floor_exact(lo*raw)
    assert a==floor_exact(hi*raw),(b,'insufficient digit precision')
    integer=floor_exact(lo)
    assert integer==floor_exact(hi)
    frac_digits=[]
    for _ in range(length):
        a,d=divmod(a,b)
        frac_digits.append(int(d))
    frac_digits.reverse()
    assert a==integer
    zeros=[0]*(length+1)
    for j in range(length-1,-1,-1):
        zeros[j]=(zeros[j+1]+1) if frac_digits[j]==0 else 0
    p,rawq=g.mpz(integer),g.mpz(1)
    best,best_noninteger,best_deep=None,None,None
    states={}
    curves=[]
    max_excess=None
    max_ratio=None
    max_cancel=None
    overlaps=0
    for m,d in enumerate(frac_digits[:digits],start=1):
        p=p*b+d
        rawq*=b
        divisor=g.gcd(p,rawq)
        q=rawq//divisor
        error_lo=lo-g.mpq(p,rawq)
        error_hi=hi-g.mpq(p,rawq)
        assert error_lo>0
        cl,ch=q*q*error_lo,q*q*error_hi
        z=zeros[m]
        censored=z==length-m
        value=None
        key=(p//divisor,q)
        def candidate():
            nonlocal value
            if value is None:
                value=row(m,p,rawq,divisor,cl,ch,z,censored)
            return value
        for label,allowed in [('all',True),('q_gt_1',q>1),('m_ge_32',m>=32)]:
            if not allowed:
                continue
            old=states.get(label)
            if old is None or ch<old[0]:
                states[label]=(cl,ch,key,candidate())
            elif cl<=old[1] and key!=old[2]:
                overlaps+=1
                raise AssertionError((b,m,'unresolved quality comparison'))
        excess=z-m
        if max_excess is None or excess>max_excess['zero_excess_z_minus_m']:
            max_excess={'m':m,'following_zeros_observed':z,
                        'zero_excess_z_minus_m':excess,'right_censored':censored}
        ratio=z/m
        if max_ratio is None or ratio>max_ratio['zero_ratio_z_over_m']:
            max_ratio={'m':m,'following_zeros_observed':z,
                       'zero_ratio_z_over_m':ratio,'right_censored':censored}
        if max_cancel is None or divisor>max_cancel[0]:
            max_cancel=(divisor,{'m':m,'log10_gcd_cancellation':log10(divisor),
                'log10_raw_denominator':log10(rawq),'log10_reduced_q':log10(q)})
        if b in CURVE_BASES:
            curves.append({'m':m,'log10_error_lower':log10(error_lo),
                'log10_error_upper':log10(error_hi),'log10_reduced_q':log10(q),
                'log10_C_lower':log10(cl),'log10_C_upper':log10(ch),
                'following_zeros_observed':z,'zero_run_right_censored':censored,
                'log10_gcd_cancellation':log10(divisor)})
    result={'base':b,'fractional_digits_scored':digits,'additional_digits_checked':lookahead,
        'best_reduced_C':states['all'][3],
        'best_reduced_C_q_gt_1':states['q_gt_1'][3],
        'best_reduced_C_m_ge_32':states['m_ge_32'][3],
        'max_zero_excess':max_excess,'max_zero_ratio':max_ratio,
        'maximum_gcd_cancellation':max_cancel[1],
        'unresolved_rank_comparisons':overlaps}
    return result,curves


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--digits',type=int,default=2000)
    parser.add_argument('--max-base',type=int,default=64)
    parser.add_argument('--lookahead',type=int,default=64)
    args=parser.parse_args()
    assert 32<=args.digits<=2000 and 10<=args.max_base<=64 and 1<=args.lookahead<=128
    bits=math.ceil((args.digits+args.lookahead)*math.log2(args.max_base))+256
    targets=target_bounds(bits)
    out={'parameters':vars(args),'binary_precision_pi_sqrt2':bits,
         'target_exact_bounds_hex':{name:{'lower':rational_hex(lo),'upper':rational_hex(hi)}
                                    for name,(lo,hi) in targets.items()},
         'targets':{},
         'scope':['Finite truncation scan; does not cover all rational approximations.',
                  'Every integer base was scanned with identical positional budget, not identical denominator size.',
                  'Selected minimum over63bases is descriptive, not out-of-sample predictive performance.',
                  'Zero runs reaching the inspected endpoint are lower bounds and are explicitly censored.',
                  'Liouville constant is sum10^(-k!) for k>=1; the special decimal structure is known in advance.',
                  'C comparison uses reduced denominators and exact rational intervals; logarithms are display values.']}
    curveout={'parameters':vars(args),'targets':{}}
    with g.context(precision=100):
        for name,target in targets.items():
            rows=[]
            curveout['targets'][name]={}
            for b in range(2,args.max_base+1):
                summary,curves=scan(target,b,args.digits,args.lookahead)
                rows.append(summary)
                if curves:
                    curveout['targets'][name][str(b)]=curves
            out['targets'][name]=rows
            print(name,'done',flush=True)
    (HERE/'base_scan_results.json').write_text(json.dumps(out,indent=2)+'\n')
    (HERE/'selected_curves.json').write_text(json.dumps(curveout,separators=(',',':'))+'\n')
    print('RESULTS_READY',flush=True)


if __name__=='__main__':
    main()
