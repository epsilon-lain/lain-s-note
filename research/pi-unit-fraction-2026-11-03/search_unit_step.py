#!/usr/bin/env python3
"""Certify the best single unit-fraction correction that stays on the same side.

Standard library only. Full finite enumeration needs factoring the starting
denominator; this small research CLI bounds it to 40 bits to avoid long runs.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path

from experiment import atan_bounds, decimal
from finite_update_check import solve


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--start', default='22/7', help='Reduced automatically; e.g. 22/7 or 355/113')
    parser.add_argument('--terms', type=int, default=128, help='Alternating atan terms per input; increase if a comparison is unresolved')
    parser.add_argument('--output', type=Path, help='Optional exact certificate JSON path')
    args = parser.parse_args()
    x = Fraction(args.start)
    if x.denominator.bit_length() > 40:
        parser.error('This trial-division version accepts starting denominators of at most 40 bits.')
    if not 16 <= args.terms <= 4096:
        parser.error('--terms must be between 16 and 4096.')
    a, b = atan_bounds(5, args.terms)
    c, d = atan_bounds(239, args.terms)
    lo, hi = 16*a-4*d, 16*b-4*c
    if not (x < lo or x > hi):
        parser.error('Target interval does not certify the side; increase --terms.')
    try:
        sign, rl, rh, best, rows = solve(x, lo, hi)
    except AssertionError as exc:
        raise SystemExit('A comparison was unresolved; increase --terms. '+str(exc))
    cl, ch, unit_d, det, new = best
    assert new == x+Fraction(sign, unit_d)
    assert cl > 0 and all(ch < row[0] or unit_d == row[2] for row in rows)
    def rational(v):
        # Hex avoids Python's decimal conversion cap for certified large integers.
        return {'numerator_hex': hex(v.numerator), 'denominator_hex': hex(v.denominator)}
    out = {
        'start': str(x), 'sign': sign, 'result': str(new),
        'unit_denominator': str(unit_d),
        'branch_determinant': str(det), 'divisor_candidates': len(rows),
        'old_quality_approx': decimal(x.denominator**2*(rl+rh)/2),
        'new_quality_approx': decimal((cl+ch)/2),
        'strict_improvement_certified': ch < x.denominator**2*rl,
        'pi_exact_bounds': {'lower': rational(lo), 'upper': rational(hi)},
        'new_quality_exact_bounds': {'lower': rational(cl), 'upper': rational(ch)},
        'scope': 'Global minimum over one unit correction toward pi that remains strictly on the original side; not a multi-step or global rational approximation optimum.',
        'pi_method': f'16 atan(1/5)-4 atan(1/239), {args.terms} alternating-series terms, exact rational remainder bounds'
    }
    if args.output:
        args.output.write_text(json.dumps(out, indent=2)+'\n')
    display = {k:v for k,v in out.items() if not k.endswith('exact_bounds')}
    print(json.dumps(display, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
