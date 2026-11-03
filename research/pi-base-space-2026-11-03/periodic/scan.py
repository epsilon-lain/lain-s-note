#!/usr/bin/env python3
"""Generate pi approximants from bounded-base, short purely periodic templates.

No continued fraction input. For each fixed base b and length ell, read exactly
the first ell fractional digits, then repeat that word forever. Exact Machin
intervals certify digits, reduced fractions, matching lengths, and quality.
Standard library only. This is a bounded search, not a universal predictor.
"""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
from math import gcd, log2, log
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def atan_bounds(inv, n=100):
    s = sum((F((-1)**k, (2*k+1)*inv**(2*k+1)) for k in range(n)), F(0))
    term = F((-1)**n, (2*n+1)*inv**(2*n+1))
    return min(s, s+term), max(s, s+term)


def decimal(x):
    with localcontext() as ctx:
        ctx.prec = 35
        return str(Decimal(x.numerator)/Decimal(x.denominator))


def exact_pair(x):
    return [hex(x.numerator), hex(x.denominator)]


def digits(lo, hi, base, count):
    integer = lo.numerator//lo.denominator
    assert integer == hi.numerator//hi.denominator
    lo, hi = lo-integer, hi-integer
    out = []
    for _ in range(count):
        lo, hi = base*lo, base*hi
        digit = lo.numerator//lo.denominator
        assert digit == hi.numerator//hi.denominator, 'Increase --atan-terms.'
        assert 0 <= digit < base
        out.append(digit)
        lo, hi = lo-digit, hi-digit
    return integer, out


def quality(candidate, lo, hi):
    q = candidate.denominator
    if candidate < lo:
        return q*q*(lo-candidate), q*q*(hi-candidate), 'below'
    assert candidate > hi, 'Increase --atan-terms to decide candidate side.'
    return q*q*(candidate-hi), q*q*(candidate-lo), 'above'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-base', type=int, default=512)
    parser.add_argument('--max-period', type=int, default=4)
    parser.add_argument('--digits', type=int, default=32)
    parser.add_argument('--atan-terms', type=int, default=100)
    parser.add_argument('--outdir', type=Path, default=HERE)
    args = parser.parse_args()
    assert 2 <= args.max_base <= 4096 and 1 <= args.max_period <= 16
    assert args.digits >= max(args.max_period, 12)
    args.outdir.mkdir(exist_ok=True, parents=True)
    a,b = atan_bounds(5,args.atan_terms)
    c,d = atan_bounds(239,args.atan_terms)
    lo,hi = 16*a-4*d,16*b-4*c
    internal, special, best_by_base = [], {}, []
    for base in range(2,args.max_base+1):
        integer, stream = digits(lo,hi,base,args.digits)
        assert integer == 3
        word_value = 0
        base_rows = []
        for ell in range(1,args.max_period+1):
            word_value = base*word_value + stream[ell-1]
            raw_q = base**ell-1
            cancellation = gcd(word_value,raw_q)
            candidate = F(3*raw_q+word_value,raw_q)
            assert candidate.denominator == raw_q//cancellation
            cl,ch,side = quality(candidate,lo,hi)
            word = stream[:ell]
            matched = 0
            for n,digit in enumerate(stream):
                if digit != word[n%ell]:
                    break
                matched += 1
            row = {
                'base':base,'period':ell,'word_digits':word,
                'p':str(candidate.numerator),'q':str(candidate.denominator),
                'raw_q':str(raw_q),'gcd_cancelled':str(cancellation),
                'quality_lower_approx':decimal(cl),'quality_upper_approx':decimal(ch),
                'side':side,'literal_prefix_match_digits':matched,
                'match_is_only_lower_bound':matched==len(stream),
                'predicted_digits_after_training_word':matched-ell,
                'base_cost_bits':log2(base),'word_cost_bits':ell*log2(base),
                'template_cost_bits_base_plus_word':(ell+1)*log2(base),
                'reduced_denominator_bits_log2':log2(candidate.denominator),
                'cancellation_gain_base_digits':2*log(cancellation)/log(base),
                'digit_match_quality_upper_bound_approx':decimal(F(candidate.denominator**2,base**matched)),
            }
            assert ch < F(candidate.denominator**2,base**matched)
            item = (candidate,cl,ch,row)
            internal.append(item)
            base_rows.append(item)
        best = min(base_rows,key=lambda x:x[2])
        best_by_base.append(best[3])
        if base in [8,10,16,114]:
            special[str(base)] = {
                'first_12_fractional_digits':stream[:12],
                'one_digit_periodic_template':base_rows[0][3],
                'best_of_periods_1_to_limit':best[3]}
    # Collapse representations of the same rational before ranking.
    grouped = {}
    for cand,cl,ch,row in internal:
        if cand not in grouped:
            grouped[cand]=[cl,ch,row,[]]
        grouped[cand][3].append({'base':row['base'],'period':row['period']})
        if row['template_cost_bits_base_plus_word'] < grouped[cand][2]['template_cost_bits_base_plus_word']:
            grouped[cand][2]=row
    ranked = sorted(grouped.values(),key=lambda x:x[1])
    for left,right in zip(ranked,ranked[1:]):
        assert left[1] < right[0], 'Ranking unresolved; increase --atan-terms.'
    top=[]
    for cl,ch,row,reps in ranked[:20]:
        out=dict(row)
        out['representations_in_scan']=reps
        out['quality_exact_bounds_hex']={'lower':exact_pair(cl),'upper':exact_pair(ch)}
        top.append(out)
    predictions=sorted((x[3] for x in internal),key=lambda x:(-x['predicted_digits_after_training_word'],x['base'],x['period']))[:20]
    result={
        'method':'Pure periodic templates from only first ell digits; all digit and score decisions certified by exact Machin intervals; no CF input.',
        'budget':{'bases_inclusive':[2,args.max_base],'periods_inclusive':[1,args.max_period],
                  'digits_certified_per_base':args.digits,'atan_terms_per_input':args.atan_terms},
        'pi_bounds_hex':{'lower':exact_pair(lo),'upper':exact_pair(hi)},
        'candidate_templates':len(internal),'unique_rational_candidates':len(grouped),
        'top_quality_unique':top,'top_literal_predictions':predictions,
        'special_bases':special,
        'cautions':[
            'Bases and periods have a fixed finite budget; the best observed base is selected on these same data, not a blind prediction.',
            'Matched digit continuation is descriptive. Neither a long match nor this finite scan proves an infinite property of pi.',
            'Template costs are an explicit simple accounting convention, not a claim of optimal Kolmogorov complexity.',
            'The true score uses the reduced rational denominator; raw repeated-word denominator and digit length alone can overcharge severely.',
            'All decimals are display approximations. Exact pi bounds and p/q determine exact score bounds for every row.',
            'Any p/q=a+k/q becomes one-digit periodic in base q+1; unbounded base choice can simply encode the unknown rational denominator.'
        ]}
    (args.outdir/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    (args.outdir/'best_by_base.json').write_text(json.dumps(best_by_base,indent=2)+'\n')
    with (args.outdir/'all_candidates.jsonl').open('w') as f:
        for _,_,_,row in internal:
            f.write(json.dumps(row)+'\n')
    best=top[0]
    print(json.dumps({'templates':len(internal),'unique':len(grouped),
                     'best':{k:best[k] for k in ['base','period','word_digits','p','q','quality_upper_approx','literal_prefix_match_digits']},
                     'special_digits':{k:v['first_12_fractional_digits'] for k,v in special.items()}},indent=2))


if __name__=='__main__':
    main()
