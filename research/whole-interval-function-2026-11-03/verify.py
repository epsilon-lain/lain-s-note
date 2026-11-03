#!/usr/bin/env python3
"""Exact whole-interval certificates for the alternating-rectangle model.

Standard library only. The curve construction uses its affine fixed point;
verification also integrates the original rectangles independently.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import random

HERE = Path(__file__).resolve().parent


def polygon(b, depth):
    """Exact values on x=j/b**(2*depth), with denominator b**(3*depth)."""
    nums = [0, 0]
    grid = scale = 1
    for _ in range(depth):
        out = []
        for cell in range(b*b):
            ell = cell % b
            for i in range(grid):
                baseline = F(ell*(b-ell), b*b) + F((b-2*ell-1)*i, b*b*grid)
                val = baseline * (scale*b**3) + nums[i]
                assert val.denominator == 1
                out.append(val.numerator)
        nums = out + [0]
        grid *= b*b
        scale *= b**3
    return nums, grid, scale


def rectangle_interval(x, b, depth):
    """Integrate filled rectangles, then bound the remaining active box."""
    assert 0 <= x <= 1
    if x == 1:
        return F(0), F(0)
    r = x
    y = F(0)
    width = height = F(1)
    area = moment = F(0)
    for step in range(depth):
        r *= b
        digit = r.__floor__()
        r -= digit
        if step % 2 == 0:
            filled_width = width * digit / b
            piece_area = filled_width * height
            moment += piece_area * (y + height/2)
            area += piece_area
            width /= b
        else:
            filled_height = height * digit / b
            piece_area = width * filled_height
            moment += piece_area * (y + filled_height/2)
            area += piece_area
            y += filled_height
            height /= b
    tail_area = x - area
    assert 0 <= tail_area <= width*height
    prefix = b*(area-2*moment)
    bounds = [b*tail_area*(1-2*v) for v in (y, y+height)]
    return prefix+min(bounds), prefix+max(bounds)


def rational_value(x, b):
    """Independent exact summation using the finite two-digit residue orbit."""
    if x == 1:
        return F(0)
    a, q = x.numerator, x.denominator
    seen, terms = {}, []
    s = a
    psi = lambda t: t*(1-t)
    while s not in seen:
        seen[s] = len(terms)
        t = (b*s) % q
        next_s = (b*b*s) % q
        terms.append(psi(F(t,q))-psi(F(next_s,q))/b**2)
        s = next_s
    pre = seen[s]
    period = len(terms)-pre
    z = F(1,b**3)
    return (sum((z**j*terms[j] for j in range(pre)), F(0))
            + z**pre*sum((z**j*terms[pre+j] for j in range(period)), F(0))
            / (1-z**period))


def maximum(b):
    return F(1,4) if b % 2 == 0 else F(b*(b+1),4*(b*b+b+1))


def mean(b):
    return F(b*(b+1),6*(b*b+b+1))


def main():
    rng = random.Random(20261103)
    checks = dict(exact_rectangle_node_checks=0, rational_uniform_band_checks=0,
                  odd_maximum_witness_checks=0, general_extremum_checks=0,
                  exact_polygon_mean_checks=0, exact_node_symmetry_checks=0)
    curves = []
    for b, depth in [(2,6),(3,4)]:
        nums, grid, scale = polygon(b, depth)
        assert nums == list(reversed(nums))
        checks['exact_node_symmetry_checks'] += len(nums)
        assert min(nums) == 0
        assert [j for j,v in enumerate(nums) if v == 0] == [k*grid//b for k in range(b+1)]
        lipschitz = max(abs(nums[j+1]-nums[j])*F(grid,scale) for j in range(grid))
        assert lipschitz == b*(1-F(1,b**depth))
        avg = F(2*sum(nums)-nums[0]-nums[-1],2*grid*scale)
        assert avg == mean(b)*(1-F(1,b**(3*depth)))
        checks['exact_polygon_mean_checks'] += 1
        indices = sorted({0,grid,*[rng.randrange(grid) for _ in range(260)]})
        for j in indices:
            lo,hi = rectangle_interval(F(j,grid),b,2*depth)
            assert lo == hi == F(nums[j],scale)
            checks['exact_rectangle_node_checks'] += 1
        # Rationals away from nodes: compare the entire certified band against
        # independently integrated rectangles, with exact rational arithmetic.
        for _ in range(160):
            x = F(rng.randrange(997),997)
            cell = (x*grid).__floor__()
            u = x*grid-cell
            lower = F(nums[cell],scale)*(1-u)+F(nums[cell+1],scale)*u
            exact = rational_value(x,b)
            assert lower <= exact <= lower+F(1,4*b**(3*depth))
            lo,hi = rectangle_interval(x,b,24)
            assert lo <= exact <= hi
            checks['rational_uniform_band_checks'] += 1
        if b % 2 == 0:
            assert max(nums)*F(1,scale) == maximum(b)
        curves.append({'base':b,'depth_pairs':depth,'grid_denominator':grid,
                       'value_denominator':scale,'value_numerators':nums,
                       'uniform_error_upper_exact':str(F(1,4*b**(3*depth))),
                       'maximum_exact':str(maximum(b)),'mean_exact':str(mean(b)),
                       'polygon_mean_exact':str(avg)})
    maxima = []
    b, depth = 3, 4
    for word in range(b**depth):
        odd_digits = []
        t = word
        for _ in range(depth):
            odd_digits.append(t % b)
            t //= b
        prefix = 0
        for digit in reversed(odd_digits):
            prefix = prefix*b+digit
            prefix = prefix*b+1
        x = F(2*prefix+1,2*b**(2*depth))
        assert rational_value(x,b) == maximum(b)
        lo,hi = rectangle_interval(x,b,24)
        assert lo <= maximum(b) <= hi
        checks['odd_maximum_witness_checks'] += 1
        maxima.append(str(x))
    for b in range(2,16):
        witness = F(1,2*b) if b % 2 == 0 else F(1,2*(b+1))
        assert rational_value(witness,b) == maximum(b)
        lo,hi = rectangle_interval(witness,b,24)
        assert lo <= maximum(b) <= hi
        checks['general_extremum_checks'] += 1
    output = {'scope':'Fixed integer base, all area/address values x in [0,1].',
              'certification':'Exact rational nodes and a proved uniform bound, not only samples.',
              'recursion':'F((k*b+ell+u)/b^2)=[ell*(b-ell)+(b-2*ell-1)*u]/b^2+F(u)/b^3',
              'checks':checks,'curves':curves,
              'ternary_maximum_witnesses_exact':maxima}
    (HERE/'results.json').write_text(json.dumps(output,separators=(',',':'))+'\n')
    print(json.dumps(checks,indent=2))


if __name__ == '__main__':
    main()
