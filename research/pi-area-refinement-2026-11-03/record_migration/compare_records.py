#!/usr/bin/env python3
"""Compare a certified small example, an actual local record, and a reported record.

The reported record's conclusion is conditional on its published CF coefficient.
This script does not reconstruct its numerator or denominator.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def dec(x):
    with localcontext() as ctx:
        ctx.prec = 45
        return str(Decimal(x.numerator) / Decimal(x.denominator))


def main():
    old = json.loads((HERE.parent / 'results.json').read_text())
    def rat(pair):
        return F(int(pair[0], 16), int(pair[1], 16))
    lo = rat(old['pi_exact_bounds_hex']['lower'])
    hi = rat(old['pi_exact_bounds_hex']['upper'])
    c355lo = 113**2 * (F(355, 113) - hi)
    c355hi = 113**2 * (F(355, 113) - lo)
    assert 0 < c355lo < c355hi < F(1, 113)
    A = 845188380288
    N = 383059989912
    assert N % 2 == 0
    # The following pair is a theorem IF the published next coefficient is A.
    glo, ghi = F(1, A + 2), F(1, A)
    assert 2**40 > A + 2
    local = json.loads((HERE / 'local_results.json').read_text())
    neighbors = json.loads((HERE / 'neighborhood_results.json').read_text())
    nonzero = [r for r in neighbors['neighbors_q_plus_h'] if r['h']]
    assert all(r['C_raw_gt_1_certified'] and r['C_reduced_gt_1_certified'] for r in nonzero)
    assert neighbors['checks']['all_8_multiples_reduce_to_same_p_over_q']
    minimum = min(nonzero, key=lambda r: Decimal(r['log10_C_reduced_display_bounds'][0]))
    data = {
        'definition': 'C=q^2*abs(pi-p/q), relative difference inside one q^-2 cell',
        '355_over_113': {
            'side': 'above_pi', 'C_display_bounds': [dec(c355lo), dec(c355hi)],
            'fractional_digits': [15,112,112,69], 'evidence': 'certified rational pi interval'},
        'local_record': {
            'side': 'below_pi', 'q_decimal_digits': local['q']['decimal_digits'],
            'C_display_bounds': local['C_display_bounds'],
            'fractional_digits_symbolic': ['p-3*q', '0', 'floor(q*C)'],
            'evidence': 'actual p/q and first three radix-q digits independently certified'},
        'reported_record': {
            'following_partial_quotient_A': str(A), 'convergent_index_zero_based': N,
            'side_if_report_correct': 'below_pi',
            'C_strict_exact_bounds': [f'1/{A+2}', f'1/{A}'],
            'C_display_bounds': [dec(glo), dec(ghi)],
            'ratio_C355_to_reported_C_display_bounds': [dec(c355lo/ghi), dec(c355hi/glo)],
            'q_decimal_digits_conservative_lower_bound': 57458998486,
            'actual_p_q_reconstructed_in_this_work': False,
            'fractional_digits_symbolic_if_report_correct': ['p-3*q', '0', 'floor(q*C)'],
            'evidence': 'conditional theorem using published A and index, not a giant-integer calculation'},
        'neighborhood': {
            'scope': 'local record only; integer B=q+h, -32<=h<=32',
            'all_64_nonzero_offsets_fail_C_less_than_1_after_reduction': True,
            'best_nonzero_offset_after_reduction': minimum['h'],
            'best_nonzero_log10_C_reduced': minimum['log10_C_reduced_display_bounds'][0],
            'multiples_m_1_through_8_all_reduce_to_same_fraction': True},
        'interpretation': 'Transfer of the representation succeeds. No broad favorable neighborhood or new predictor was found. Selecting q from a known good fraction already supplies the favorable scale.'}
    (HERE / 'comparison.json').write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
