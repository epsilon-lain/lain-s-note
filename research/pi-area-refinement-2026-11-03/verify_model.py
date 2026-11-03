#!/usr/bin/env python3
"""Exact checks for the ordered area-refinement model; standard library only."""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def atan_bounds(inv, count=100):
    total = sum((F((-1)**k, (2*k+1)*inv**(2*k+1)) for k in range(count)), F(0))
    next_term = F((-1)**count, (2*count+1)*inv**(2*count+1))
    return min(total, total+next_term), max(total, total+next_term)


def decimal(x):
    with localcontext() as ctx:
        ctx.prec = 40
        return str(Decimal(x.numerator)/Decimal(x.denominator))


def exact(x):
    return [hex(x.numerator), hex(x.denominator)]


def digits(lo, hi, base, count):
    integer = lo.numerator//lo.denominator
    assert integer == hi.numerator//hi.denominator
    lo, hi = lo-integer, hi-integer
    out = []
    for _ in range(count):
        lo, hi = base*lo, base*hi
        digit = lo.numerator//lo.denominator
        assert digit == hi.numerator//hi.denominator
        out.append(digit)
        lo, hi = lo-digit, hi-digit
    return integer, out


def rectangle_area(rect):
    x0, y0, x1, y1 = rect
    return (x1-x0)*(y1-y0)


def refine(parent, base, digit, level):
    """Fill ordered children, then recurse inside child digit+1.

    Odd levels cut vertically, even levels horizontally. The returned filled
    rectangle is a union of digit full equal-area children, not one new digit bar.
    """
    x0, y0, x1, y1 = parent
    assert 0 <= digit < base
    if level % 2:
        step = (x1-x0)/base
        filled = (x0, y0, x0+digit*step, y1)
        active = (x0+digit*step, y0, x0+(digit+1)*step, y1)
    else:
        step = (y1-y0)/base
        filled = (x0, y0, x1, y0+digit*step)
        active = (x0, y0+digit*step, x1, y0+(digit+1)*step)
    assert rectangle_area(filled) == digit*rectangle_area(parent)/base
    assert rectangle_area(active) == rectangle_area(parent)/base
    return filled, active


def overlap_area(a, b):
    x0, y0 = max(a[0],b[0]), max(a[1],b[1])
    x1, y1 = min(a[2],b[2]), min(a[3],b[3])
    return max(F(0),x1-x0)*max(F(0),y1-y0)


def main():
    a,b = atan_bounds(5)
    c,d = atan_bounds(239)
    lo,hi = 16*a-4*d, 16*b-4*c
    integer, stream = digits(lo,hi,113,12)
    assert integer == 3 and stream[:6] == [15,112,112,69,57,4]
    parent = (F(0), F(0), F(1), F(1))
    cells, rows = [], []
    prefix, K = F(0), None
    upper = F(16,113)
    template = (F(0),F(0),upper,F(1))
    for level, digit in enumerate(stream, 1):
        filled, parent = refine(parent,113,digit,level)
        assert all(overlap_area(filled,old)==0 for old in cells)
        cells.append(filled)
        prefix += F(digit,113**level)
        assert sum(map(rectangle_area,cells),F(0)) == prefix
        assert rectangle_area(parent) == F(1,113**level)
        assert overlap_area(filled,template) == rectangle_area(filled)
        candidate_K = (upper-prefix)*113**level
        assert candidate_K.denominator == 1
        if K is not None:
            assert candidate_K == 113*K-digit
        K = candidate_K
        gap_lo = F(K-1,113**level)
        gap_hi = F(K,113**level)
        assert gap_lo < F(355,113)-hi < F(355,113)-lo < gap_hi
        rows.append({'level':level,'digit':digit,'unfilled_units_K':str(K),
                     'active_cell_area':str(rectangle_area(parent)),
                     'gap_lower_from_prefix':str(gap_lo),'gap_upper_from_prefix':str(gap_hi),
                     'filled_rectangle_exact':list(map(exact,filled)),
                     'next_active_rectangle_exact':list(map(exact,parent))})
    assert [int(r['unfilled_units_K']) for r in rows[:4]] == [1,1,1,44]
    assert all(F(r['gap_upper_from_prefix']) < F(1) for r in rows)
    # Exact discontinuity counterexample to the old canonical digit-strip metric.
    contrasts = []
    for count in [2,3,4,6,8]:
        x, y = F(1,2)-F(1,10**count), F(1,2)
        xd = [4]+[9]*(count-1)
        yd = [5]+[0]*(count-1)
        assert sum((F(v,10**k) for k,v in enumerate(xd,1)),F(0)) == x
        assert sum((F(v,10**k) for k,v in enumerate(yd,1)),F(0)) == y
        old = sum((F(abs(a-b),10**k) for k,(a,b) in enumerate(zip(xd,yd),1)),F(0))
        assert old == F(1,5)-F(1,10**count)
        contrasts.append({'decimal_positions':count,'x':decimal(x),'y':'0.5',
                          'numeric_difference':decimal(y-x),
                          'old_strip_symmetric_difference':decimal(old),
                          'new_ordered_fill_symmetric_difference':decimal(y-x)})
    out = {'method':'Exact Fraction arithmetic and alternating Machin interval',
           'pi_exact_bounds_hex':{'lower':exact(lo),'upper':exact(hi)},
           'base':113,'width_contraction':'1/113', 'pi_integer_part':integer,
           'pi_fractional_digits':stream,'refinement':rows,
           'gap_355_over_113_minus_pi':{'lower_approx':decimal(F(355,113)-hi),
                                       'upper_approx':decimal(F(355,113)-lo)},
           'old_metric_counterexamples':contrasts,
           'checks':{'digit_prefix_certified':True,'all_cell_areas_exact':True,
                     'filled_cells_pairwise_disjoint':True,'prefix_area_matches_digits':True,
                     'all_filled_cells_inside_355_template':True,
                     'all_prefix_gap_bounds_certified':True},
           'scope':'A definition and representation audit, not a new approximation record or predictor.'}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'digits':stream,'first_four_gap_units':[r['unfilled_units_K'] for r in rows[:4]],
                      'gap':out['gap_355_over_113_minus_pi'],'checks':out['checks'],
                      'old_metric_counterexamples':contrasts},indent=2))


if __name__ == '__main__':
    main()
