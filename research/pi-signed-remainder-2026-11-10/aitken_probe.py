"""Signed local-tail cancellation with full reduced-denominator cost.

This diagnostic constructs Aitken transforms of the classical Apéry pi^2
source. It generates with rational recurrence data only. Independent rational
Machin/Taylor bounds are used afterwards for exact error and cost comparisons.
No infinite-success claim and no search through the pi continued fraction.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / 'pi-small-norm-construction-2026-11-04'
sys.path.insert(0, str(OLD))
from conic_attempt import source
from construct import display_interval
from verify import atan_taylor

POINTS = [1, 2, 3, 4, 8, 16, 32, 64, 128]


def absolute_error_interval(value, lo, hi):
    if value < lo:
        return lo-value, hi-value
    assert value > hi
    return value-hi, value-lo


def aitken(r0, r1, r2):
    d0, d1 = r1-r0, r2-r1
    assert d0*d1 < 0
    rho = -d1/d0
    assert 0 < rho < 1
    T = r0+d0/(1+rho)
    assert T == (r0*r2-r1*r1)/(r0-2*r1+r2)
    return T, d0, rho


def run():
    u, v = source(max(POINTS)+3)
    raw = []
    for n in POINTS:
        rs = [6*v[j]/u[j] for j in range(n, n+4)]
        T, d, rho = aitken(*rs[:3])
        nextT, _, next_rho = aitken(*rs[1:])
        # Exact ratio-variation identity, not an asymptotic estimate.
        variation = d*rho*(1/(1+rho)-1/(1+next_rho))
        assert nextT-T == variation
        raw.append((n, rs[0], T, rho, next_rho, variation, rs[2]))

    # Construct all candidates before accessing the independent scoring bound.
    bits = 3000
    l5, h5, _ = atan_taylor(5, bits)
    l239, h239, _ = atan_taylor(239, bits)
    pl, ph = 16*l5-4*h239, 16*h5-4*l239
    lo, hi = pl*pl, ph*ph
    assert 9 < lo < hi < 10
    cases = []
    for n, R, T, rho, next_rho, variation, matched in raw:
        el, eh = absolute_error_interval(R, lo, hi)
        tl, th = absolute_error_interval(T, lo, hi)
        gain_l, gain_h = tl/eh, th/el
        factor = F(T.denominator*T.denominator, R.denominator*R.denominator)
        cost_l, cost_h = factor*gain_l, factor*gain_h
        ml, mh = absolute_error_interval(matched, lo, hi)
        matched_factor = F(T.denominator*T.denominator,
                           matched.denominator*matched.denominator)
        matched_l, matched_h = matched_factor*tl/mh, matched_factor*th/ml
        assert matched_h < 1 or matched_l > 1
        assert gain_h < 1 and cost_l > 1
        exact = {'n': n, 'R': str(R), 'T': str(T), 'rho': str(rho),
                 'rho_next': str(next_rho), 'T_next_minus_T': str(variation),
                 'source_error': list(map(str, (el, eh))),
                 'transformed_error': list(map(str, (tl, th))),
                 'error_ratio': list(map(str, (gain_l, gain_h))),
                 'arithmetic_cost_ratio': list(map(str, (cost_l, cost_h))),
                 'same_highest_input_source': str(matched),
                 'same_highest_input_cost_ratio': list(map(str, (matched_l, matched_h)))}
        digest = hashlib.sha256(json.dumps(exact, sort_keys=True,
                                          separators=(',', ':')).encode()).hexdigest()
        cases.append({'n': n, 'source': str(R), 'transformed_source': str(T),
                      'source_denominator_digits': len(str(R.denominator)),
                      'transformed_denominator_digits': len(str(T.denominator)),
                      'error_ratio_decimal_enclosure': display_interval(gain_l, gain_h),
                      'B_squared_error_ratio_decimal_enclosure': display_interval(cost_l, cost_h),
                      'analytic_error_strictly_improved': gain_h < 1,
                      'arithmetic_cost_strictly_worse': cost_l > 1,
                      'same_highest_input_source': str(matched),
                      'same_highest_input_denominator_digits': len(str(matched.denominator)),
                      'same_highest_input_cost_ratio_decimal': display_interval(matched_l, matched_h),
                      'same_highest_input_cost_improved': matched_h < 1,
                      'exact_record_sha256': digest})
    out = {'schema': 1, 'source_commit': 'c6b02f113e8527b57622d7fb8f6f60eaece91d81',
           'method': 'Aitken transform of pi-squared rational sources; pi used only for strict scoring.',
           'parameters': {'n': POINTS, 'scoring_bits': bits}, 'cases': cases,
           'summary': {'candidates': len(cases), 'analytic_improvements': len(cases),
                       'arithmetic_cost_improvements': 0,
                       'same_highest_input_cost_improvements': sum(c['same_highest_input_cost_improved'] for c in cases),
                       'last_source_denominator_digits': cases[-1]['source_denominator_digits'],
                       'last_transformed_denominator_digits': cases[-1]['transformed_denominator_digits']},
           'scope': 'Finite diagnostic of analytic acceleration versus actual reduced pi-squared source cost; not an infinite exclusion or a pi approximation theorem.'}
    (ROOT/'aitken_results.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out['summary'], indent=2))


if __name__ == '__main__':
    run()
