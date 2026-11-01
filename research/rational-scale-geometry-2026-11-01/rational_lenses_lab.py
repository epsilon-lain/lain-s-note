#!/usr/bin/env python3
"""Exact finite checks of two distinct rational observation models.

No pi digits, floating point rank decisions, or endpoint claims are used.
All targets in this script are rational controls. General statements are
proved in rational_lenses.md; finite checks do not prove them.
"""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json


def floor(x):
    return x.numerator // x.denominator


def lattice_lens(x, b):
    x, b = F(x), F(b)
    assert b > 1
    k = floor(x/b+F(1, 2))
    point = k*b
    distance = abs(x-point)
    assert point.denominator == b.denominator//gcd(k, b.denominator)
    # Independent search verifies the nearest index, including ties.
    search = range(max(0, floor(x/b)-3), max(0, floor(x/b))+5)
    assert distance == min(abs(x-j*b) for j in search)
    return {'k': k, 'point': str(point), 'distance': str(distance),
            'reduced_point_denominator': point.denominator,
            'kappa': str(point.denominator**2*distance)}


def prefix(x, b, L):
    x, b = F(x), F(b)
    assert x >= 0 and b > 1 and L >= 0
    a0 = floor(x)
    remainder = x-a0
    result = F(a0)
    digits = []
    for k in range(1, L+1):
        scaled = b*remainder
        digit = floor(scaled)
        remainder = scaled-digit
        assert 0 <= remainder < 1
        assert 0 <= digit < b
        digits.append(digit)
        result += digit/b**k
    m, n = b.numerator, b.denominator
    numerator = a0*m**L+sum(d*n**k*m**(L-k)
                              for k, d in enumerate(digits, 1))
    independent = F(numerator, m**L)
    assert result == independent
    assert x-result == remainder/b**L
    assert 0 <= x-result < 1/b**L
    assert m**L % result.denominator == 0
    return {'x': str(x), 'b': str(b), 'L': L, 'digits': digits,
            'prefix': str(result), 'remainder': str(remainder),
            'nominal_denominator': m**L,
            'reduced_denominator': result.denominator,
            'error': str(x-result),
            'kappa': str(result.denominator**2*(x-result))}


def main():
    bases = sorted({F(m, n) for n in range(1, 11)
                    for m in range(n+1, 16)})
    targets = [F(3), F(22, 7), F(355, 113), F(1, 3)]
    lens_count = 0
    for target in targets:
        for base in bases:
            lattice_lens(target, base)
            lens_count += 1
    aligned = lattice_lens(F(22, 7), F(22, 21))
    assert aligned['k'] == 3 and aligned['reduced_point_denominator'] == 7
    assert aligned['distance'] == '0'
    prefixes = []
    prefix_count = 0
    for target in targets:
        for base in [F(3, 2), F(7, 2), F(10), F(22, 7)]:
            for L in range(1, 13):
                record = prefix(target, base, L)
                if F(record['prefix']) != target:
                    D = record['reduced_denominator']
                    assert F(record['error']) >= F(1, target.denominator*D)
                    assert F(record['kappa']) >= F(D, target.denominator)
                prefix_count += 1
                if target == F(22, 7) and base == F(3, 2):
                    prefixes.append(record)
    never_terminates_control = []
    for L in range(1, 21):
        record = prefix(F(1, 3), F(3, 2), L)
        assert F(record['remainder']).denominator == 2**L
        assert F(record['error']) > 0
        never_terminates_control.append({k: record[k] for k in
            ['L', 'prefix', 'remainder', 'reduced_denominator']})
    # Same rational number, very different roles of m and n.
    b = F(3, 2)
    lattice = lattice_lens(F(22, 7), b)
    expansion = prefix(F(22, 7), b, 8)
    assert lattice['reduced_point_denominator'] == 1
    assert expansion['nominal_denominator'] == 3**8
    basin_values = [lattice_lens(F(3), b) for b in [F(3), F(3,2), F(9,4)]]
    assert [F(record['distance']) for record in basin_values] == [F(0), F(0), F(3,4)]
    output = dict(status='passed',
        scope=__doc__, nearest_lattice_checks=lens_count,
        greedy_prefix_checks=prefix_count,
        same_base_different_models=dict(base=str(b), lattice=lattice,
                                         digit_prefix=expansion),
        exact_rational_alignment=dict(target='22/7', base='22/21', **aligned),
        nonconvex_sublevel_control=dict(target='3', bases=['3','3/2','9/4'],
                                        distances=['0','0','3/4']),
        rational_prefix_series=prefixes,
        nontermination_control=never_terminates_control,
        nontermination_general_proof_in='rational_lenses.md',
        pi_claim='none; no pi numerical data in this script')
    destination = Path(__file__).with_name('rational_lenses_results.json')
    destination.write_text(json.dumps(output, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(dict(status='passed', lattice_checks=lens_count,
                         prefix_checks=prefix_count, output=destination.name)))


if __name__ == '__main__':
    main()
