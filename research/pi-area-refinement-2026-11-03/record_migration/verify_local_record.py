#!/usr/bin/env python3
"""Certify a local pi record's first three denominator-radix area digits.

Needs gmpy2. No search or global-record numerator/denominator is assumed.
Huge integers are summarized by lengths, hex prefixes/suffixes and SHA256.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import gmpy2 as g

HERE = Path(__file__).resolve().parent


def floor_exact(x):
    return x.numerator // x.denominator


def rat_pair(pair):
    return g.mpq(g.mpz(pair[0],16),g.mpz(pair[1],16))


def summary(z):
    h = z.digits(16)
    return {'bits':z.bit_length(), 'decimal_digits':len(z.digits(10)),
            'hex_prefix':h[:24], 'hex_suffix':h[-24:],
            'sha256_canonical_hex_utf8':hashlib.sha256(h.encode()).hexdigest()}


def display(x):
    with g.context(precision=180):
        return str(g.mpfr(x))


def machin_bounds(bits):
    with g.context(precision=bits,round=g.RoundDown):
        a5lo = g.atan(g.mpfr(g.mpq(1,5)))
        a239lo = g.atan(g.mpfr(g.mpq(1,239)))
    with g.context(precision=bits,round=g.RoundUp):
        a5hi = g.atan(g.mpfr(g.mpq(1,5)))
        a239hi = g.atan(g.mpfr(g.mpq(1,239)))
    with g.context(precision=bits,round=g.RoundDown):
        lo = g.mpq(16*a5lo-4*a239hi)
    with g.context(precision=bits,round=g.RoundUp):
        hi = g.mpq(16*a5hi-4*a239lo)
    return lo,hi


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    candidates=[
        HERE.parents[1]/'pi-search-2026-11-03/snapshots/terms-000010184980-0dfca99fac67/records/convergent_000000453292.json',
        HERE.parents[1]/'analysis_59e4/snapshot/records/convergent_000000453292.json',
    ]
    default_record=next((p for p in candidates if p.exists()),None)
    parser.add_argument('--record',type=Path,default=default_record,required=default_record is None)
    parser.add_argument('--output',type=Path,default=HERE/'local_results.json')
    args=parser.parse_args()
    started=time.monotonic()
    raw=args.record.read_bytes()
    rec=json.loads(raw)
    p,q=g.mpz(rec['p_hex'],16),g.mpz(rec['q_hex'],16)
    assert q>1 and g.gcd(p,q)==1
    assert q.bit_length()==775089 and rec['convergent_index_zero_based']==453292
    bits=4*q.bit_length()+512
    lo,hi=machin_bounds(bits)
    assert 3<lo<hi<4
    x=g.mpq(p,q)
    assert x<lo
    Clo=q*q*(lo-x)
    Chi=q*q*(hi-x)
    assert 0<Clo<Chi<1
    assert rat_pair(rec['quality_lower_hex'])<=Clo<Chi<=rat_pair(rec['quality_upper_hex'])
    rlo,rhi=lo-3,hi-3
    digits=[]
    for _ in range(3):
        rlo,rhi=q*rlo,q*rhi
        d=floor_exact(rlo)
        assert d==floor_exact(rhi)
        assert 0<=d<q
        digits.append(d)
        rlo,rhi=rlo-d,rhi-d
    A=p-3*q
    assert digits[0]==A and digits[1]==0
    assert digits[2]==floor_exact(q*Clo)==floor_exact(q*Chi)
    d3=digits[2]
    assert d3>0
    # Literal normalized occupancy at the third local refinement.
    occupancy=g.mpq(d3,q)
    assert occupancy<Clo<Chi<occupancy+g.mpq(1,q)
    # First two fractional digits already identify the record p/q.
    assert 3+g.mpq(digits[0],q)+g.mpq(digits[1],q*q)==x
    low3=x+g.mpq(d3,q**3)
    assert low3<lo<hi<low3+g.mpq(1,q**3)
    # Parent active area after second refinement is q^-2; pi fills C of it.
    assert q*q*(lo-x)==Clo and q*q*(hi-x)==Chi
    q10=q.digits(10)
    out={
        'input_record_sha256':hashlib.sha256(raw).hexdigest(),
        'record_index_zero_based':rec['convergent_index_zero_based'],
        'following_partial_quotient':rec['following_partial_quotient'],
        'scope':'Verified local record from uploaded snapshot, not the 845188380288 public record.',
        'pi_method':'Machin formula with directed MPFR bounds, including rounded atan inputs',
        'precision_bits':bits,'elapsed_seconds':time.monotonic()-started,
        'side':'below_pi','p':summary(p),'q':summary(q),
        'q_decimal_prefix':q10[:24],'q_decimal_suffix':q10[-24:],
        'radix':'q (width contraction 1/q)',
        'fractional_digits':[
            {'position':1,'value_summary':summary(A),'identity':'p-3*q'},
            {'position':2,'value':'0'},
            {'position':3,'value_summary':summary(d3),'identity':'floor(q*C), C=q^2*(pi-p/q)'}],
        'C_display_bounds':[display(Clo),display(Chi)],
        'third_level_filled_fraction_display':display(occupancy),
        'third_level_filled_fraction_error_upper_bound':'1/q',
        'literal_zero_tail_length_after_first_fractional_digit':1,
        'third_level_digit_fraction_is_small_but_digit_is_not_zero':True,
        'checks':{'strict_pi_side':True,'coprime_p_q':True,'uploaded_quality_interval_contains_new_interval':True,
                  'first_three_digits_certified_by_exact_integer_division':True,
                  'digits_equal_closed_form':True,'third_prefix_and_successor_strictly_bracket_pi':True},
        'interpretation':'The first fractional cell has area 1/q. Its first subcell is unresolved and no q^-2 subcells are fully filled. At the next refinement, roughly C=7.694e-8 of that q^-2 cell is filled. Huge radix packs the excess information into occupancy, not many literal zero levels.',
    }
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['record_index_zero_based','precision_bits','elapsed_seconds','side','q','C_display_bounds','literal_zero_tail_length_after_first_fractional_digit','checks']},indent=2))


if __name__=='__main__':
    main()
