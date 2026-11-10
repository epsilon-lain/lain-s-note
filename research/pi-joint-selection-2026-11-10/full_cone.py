"""Globally optimize the uniform-tail triangle certificate over all positive mixes.

This is a separate, enlarged-domain experiment: no fixed25 or square_root H
cap is imposed. Every rational source in the adjacent Apéry interval is
allowed, including rational-square sources. No decimal pi or pi CF table is
used to generate or select candidates. See THEORY.md for finite completeness.
"""
from fractions import Fraction as F
from math import gcd, isqrt
from pathlib import Path
import argparse
import hashlib
import json
import sys

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parent
LEGACY = ROOT.parent / 'pi-small-norm-construction-2026-11-04'
sys.path.insert(0, str(LEGACY))
from conic_attempt import source, root_convergents
from construct import display_interval
from verify import atan_taylor


def envelope(p, q, lower, upper):
    square = F(p*p, q*q)
    nearest = max(lower, min(square, upper))
    return q*q*(upper-lower+abs(square-nearest))/6, nearest


def witness(r0, r1, R, p, q):
    A0, B0, A1, B1 = r0.numerator, r0.denominator, r1.numerator, r1.denominator
    A, B = R.numerator, R.denominator
    det = A0*B1-A1*B0
    sign = 1 if det > 0 else -1
    xx = sign*(A*B1-B*A1)
    yy = sign*(B*A0-A*B0)
    assert xx >= 0 and yy >= 0 and xx+yy > 0
    common = gcd(xx, yy)
    x, y = xx//common, yy//common
    aa, bb = x*A0+y*A1, x*B0+y*B1
    g = gcd(aa, bb)
    assert abs(det) % common == 0 and g == abs(det)//common
    assert (aa//g, bb//g) == (A, B)
    assert gcd(x, y) == gcd(A, B) == gcd(p, q) == 1
    h = B*p*p-A*q*q
    return dict(zip(('x','y','g','A','B','h'), map(str, (x,y,g,A,B,h))))


def mathematical_record(case):
    omitted = {'certificate_decimal_enclosure','actual_quality_interval',
               'actual_quality_decimal_enclosure'}
    return {k:v for k,v in case.items() if k not in omitted}


def record_digest(record):
    return hashlib.sha256(json.dumps(record,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def run(nmax=256):
    u, v = source(nmax+1)
    cases = []
    for n in range(1, nmax+1):
        r0, r1 = 6*v[n]/u[n], 6*v[n+1]/u[n+1]
        lower, upper = sorted((r0,r1))
        delta = upper-lower
        assert 9 < lower < upper <= 10
        C0 = (upper-9)/6
        cap = 6*C0/delta
        qlimit = isqrt(cap.numerator//cap.denominator)
        candidates = set()
        checked = 0
        for endpoint in (lower, upper):
            for p, q in root_convergents(endpoint.numerator, endpoint.denominator, qlimit):
                if 3*q <= p <= 4*q:
                    candidates.add((p,q))
                    checked += 1
        assert (3,1) in candidates
        value, q, p, R = min((envelope(p,q,lower,upper)[0], q, p,
                             envelope(p,q,lower,upper)[1]) for p,q in candidates)
        assert value <= C0
        w = witness(r0,r1,R,p,q)
        assert value == (F(abs(int(w['h'])),int(w['B']))+q*q*delta)/6
        cases.append({'n':n, 'source0':str(r0), 'source1':str(r1),
                      'L':str(lower), 'U':str(upper), 'delta':str(delta),
                      'C0':str(C0), 'q_limit':str(qlimit),
                      'endpoint_convergents_checked':checked,
                      'unique_candidates':len(candidates),
                      'p':str(p), 'q':str(q), 'certificate':str(value),
                      'certificate_decimal_enclosure':display_interval(value,value),
                      'R':str(R), 'witness':w,
                      'position':'below' if F(p*p,q*q)<lower else
                                 ('above' if F(p*p,q*q)>upper else 'inside')})
    bits = max(2*int(c['q']).bit_length()+256 for c in cases)
    l5,h5,_ = atan_taylor(5,bits)
    l239,h239,_ = atan_taylor(239,bits)
    pl,pu = 16*l5-4*h239,16*h5-4*l239
    for c in cases:
        p,q = int(c['p']),int(c['q'])
        if F(p,q) < pl:
            ql,qh = q*(q*pl-p),q*(q*pu-p)
        else:
            assert F(p,q) > pu
            ql,qh = q*(p-q*pu),q*(p-q*pl)
        assert 0 < ql < qh < F(c['certificate'])
        c['actual_quality_interval'] = [str(ql),str(qh)]
        c['actual_quality_decimal_enclosure'] = display_interval(ql,qh)
    best = min(cases,key=lambda c:(F(c['certificate']),int(c['q']),int(c['p'])))
    actual_best = min(cases,key=lambda c:F(c['actual_quality_interval'][1]))
    # The interval widths are much smaller than all cross-candidate gaps.
    for c in cases:
        if (c['p'],c['q']) != (actual_best['p'],actual_best['q']):
            assert F(actual_best['actual_quality_interval'][1]) < F(c['actual_quality_interval'][0])
    fractions = []
    fraction_ids = {}
    compact_cases = []
    witness_records = []
    previous_fraction = None
    for c in cases:
        key = (c['p'],c['q'])
        if key not in fraction_ids:
            fraction_ids[key] = len(fractions)
            fractions.append({'id':fraction_ids[key], 'p':c['p'], 'q':c['q'],
                              'actual_quality_interval':c['actual_quality_interval'],
                              'actual_quality_decimal_enclosure':c['actual_quality_decimal_enclosure']})
        exact = mathematical_record(c)
        compact_cases.append({'n':c['n'], 'fraction_id':fraction_ids[key],
                              'q_limit_bits':int(c['q_limit']).bit_length(),
                              'endpoint_convergents_checked':c['endpoint_convergents_checked'],
                              'unique_candidates':c['unique_candidates'], 'position':c['position'],
                              'certificate_decimal_enclosure':c['certificate_decimal_enclosure'],
                              'exact_record_sha256':record_digest(exact)})
        if key != previous_fraction or c['n'] == nmax:
            witness_records.append(exact)
        previous_fraction = key
    out = {'schema':2, 'experiment':'Unrestricted positive-source cone; global uniform-tail certificate optimum.',
           'parameters':{'n_min':1,'n_max':nmax,'fraction_domain':'coprime p,q; q>0; 3<=p/q<=4',
                         'fixed_norm_budgets':False},
           'source_commit':'88c331c410290c5eb915776ecd2415bffff20dc0',
           'digest_encoding':'sha256(json.dumps(mathematical_record,sort_keys=True,separators=(comma,colon)).encode(utf-8)); record fields are illustrated in witness_records',
           'cases':compact_cases, 'distinct_fractions':fractions, 'witness_records':witness_records,
           'summary':{'source_intervals':nmax,
                      'endpoint_convergents_checked':sum(c['endpoint_convergents_checked'] for c in cases),
                      'unique_candidates_per_interval_sum':sum(c['unique_candidates'] for c in cases),
                      'distinct_selected_fractions':len({(c['p'],c['q']) for c in cases}),
                      'inside_selected':sum(c['position']=='inside' for c in cases),
                      'best_certificate_n':best['n'], 'best_certificate':best['certificate'],
                      'best_certificate_decimal_enclosure':best['certificate_decimal_enclosure'],
                      'best_certificate_fraction':best['p']+'/'+best['q'],
                      'best_actual_fraction':actual_best['p']+'/'+actual_best['q'],
                      'best_actual_n':actual_best['n'],
                      'best_actual_quality_decimal_enclosure':actual_best['actual_quality_decimal_enclosure'],
                      'scoring_bits':bits},
           'scope':'Global only for this uniform-tail triangle certificate at each listed finite source interval. No global pi-record or infinite-family claim.'}
    (ROOT/'full_cone_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out['summary'],indent=2))
    return out


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n-max',type=int,default=256)
    args=parser.parse_args()
    assert args.n_max >= 1
    run(args.n_max)
