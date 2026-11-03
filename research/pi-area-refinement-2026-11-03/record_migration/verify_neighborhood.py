#!/usr/bin/env python3
"""Certify integer radices q+h (-32<=h<=32) and m*q (1<=m<=8).

Uses fresh directed Machin bounds and exact rational/integer arithmetic.
Output contains summaries only, never the 233,325-digit denominator.
Requires gmpy2 and adjacent verify_local_record.py.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import gmpy2 as g

from verify_local_record import HERE, display, floor_exact, machin_bounds, summary


def log10_display(x):
    assert x > 0
    with g.context(precision=180):
        return str(g.log10(g.mpfr(x)))


def small_or_summary(z):
    return int(z) if z.bit_length() < 40 else summary(z)


def certify(B, lo, hi):
    blo, bhi = B*lo, B*hi
    P = floor_exact(blo+g.mpq(1,2))
    assert P == floor_exact(bhi+g.mpq(1,2))
    elo, ehi = blo-P, bhi-P
    assert -g.mpq(1,2) < elo < ehi < g.mpq(1,2)
    if elo > 0:
        side, dlo, dhi = 'below_pi', elo, ehi
    else:
        assert ehi < 0
        side, dlo, dhi = 'above_pi', -ehi, -elo
    common = g.gcd(P,B)
    rawlo, rawhi = B*dlo, B*dhi
    redlo, redhi = rawlo/common**2, rawhi/common**2
    return P, common, rawlo, rawhi, redlo, redhi, {
        'side': side,
        'gcd_P_B': small_or_summary(common),
        'distance_to_integer_display_bounds': [display(dlo), display(dhi)],
        'C_raw_display_bounds': [display(rawlo), display(rawhi)],
        'C_reduced_display_bounds': [display(redlo), display(redhi)],
        'log10_C_raw_display_bounds': [log10_display(rawlo), log10_display(rawhi)],
        'log10_C_reduced_display_bounds': [log10_display(redlo), log10_display(redhi)],
        'C_raw_gt_1_certified': bool(rawlo > 1),
        'C_reduced_gt_1_certified': bool(redlo > 1),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    candidates=[
        HERE.parents[1]/'pi-search-2026-11-03/snapshots/terms-000010184980-0dfca99fac67/records/convergent_000000453292.json',
        HERE.parents[1]/'analysis_59e4/snapshot/records/convergent_000000453292.json',
    ]
    default_record=next((p for p in candidates if p.exists()),None)
    parser.add_argument('--record',type=Path,default=default_record,required=default_record is None)
    parser.add_argument('--output', type=Path, default=HERE/'neighborhood_results.json')
    args = parser.parse_args()
    started = time.monotonic()
    raw = args.record.read_bytes()
    rec = json.loads(raw)
    p, q = g.mpz(rec['p_hex'],16), g.mpz(rec['q_hex'],16)
    assert rec['convergent_index_zero_based'] == 453292
    assert q.bit_length() == 775089 and g.gcd(p,q) == 1
    bits = 4*q.bit_length()+512
    lo, hi = machin_bounds(bits)
    assert 3 < lo < hi < 4 and g.mpq(p,q) < lo
    Clo, Chi = q*(q*lo-p), q*(q*hi-p)
    assert 0 < Clo < Chi < 1
    neighbors = []
    for h in range(-32,33):
        B = q+h
        P, common, rawlo, rawhi, redlo, redhi, row = certify(B,lo,hi)
        row = {'h':h, 'nearest_numerator_offset_from_p':int(P-p), **row}
        # q*pi=p+epsilon, so the local numerator shift equals nearest(h*pi).
        if h:
            assert rawlo > 1
            assert redlo > 1
            row['different_reduced_fraction_from_p_over_q'] = bool(P*q != p*B)
            assert row['different_reduced_fraction_from_p_over_q']
            hlo,hhi = (h*lo,h*hi) if h>0 else (h*hi,h*lo)
            assert floor_exact(hlo+g.mpq(1,2)) == floor_exact(hhi+g.mpq(1,2)) == P-p
        else:
            assert P == p and common == 1
        neighbors.append(row)
    multiples = []
    for m in range(1,9):
        B = m*q
        P, common, rawlo, rawhi, redlo, redhi, row = certify(B,lo,hi)
        assert P == m*p and common == m
        assert P//common == p and B//common == q
        assert rawlo == m*m*Clo and rawhi == m*m*Chi
        assert redlo == Clo and redhi == Chi
        multiples.append({'m':m, **row, 'reduces_to_same_p_over_q':True})
    nonzero = [r for r in neighbors if r['h']]
    out = {
        'input_record_sha256':hashlib.sha256(raw).hexdigest(),
        'record_index_zero_based':rec['convergent_index_zero_based'],
        'scope':'Uploaded local record only; no global-record denominator is assumed.',
        'q':summary(q), 'precision_bits':bits,
        'method':'Outward Machin/MPFR pi bounds; exact nearest-integer certification, exact gcd and rational quality intervals. Displayed decimals/logarithms are approximations of certified endpoints.',
        'C_definition':'C_raw=B*abs(B*pi-P), C_reduced=C_raw/gcd(P,B)^2, P=nearest(B*pi)',
        'neighbors_q_plus_h':neighbors,
        'multiples_m_times_q':multiples,
        'checks':{
            'all_nearest_integers_and_residual_signs_certified':True,
            'all_64_nonzero_h_have_C_raw_gt_1':all(r['C_raw_gt_1_certified'] for r in nonzero),
            'all_64_nonzero_h_have_C_reduced_gt_1':all(r['C_reduced_gt_1_certified'] for r in nonzero),
            'all_8_multiples_reduce_to_same_p_over_q':True,
            'all_neighbor_numerator_shifts_equal_nearest_h_times_pi':True,
        },
        'elapsed_seconds':time.monotonic()-started,
    }
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'checks':out['checks'],'elapsed_seconds':out['elapsed_seconds'],
        'neighbor_gcds':{str(r['h']):r['gcd_P_B'] for r in neighbors if r['gcd_P_B']!=1},
        'sample_neighbors':[r for r in neighbors if r['h'] in [-32,-7,-1,0,1,7,32]],
        'multiples':[{'m':r['m'],'gcd':r['gcd_P_B'],'log10_C_raw':r['log10_C_raw_display_bounds'][0], 'same_fraction':r['reduces_to_same_p_over_q']} for r in multiples]},indent=2))


if __name__ == '__main__':
    main()
