#!/usr/bin/env python3
"""Exhaustive decimal-denominator scan with independent, exact pi bounds.

Every terminating decimal in lowest terms has denominator 2**a * 5**c.
For each such 1 < q <= 10**120, the nearest coprime numerator on each
side of q*pi is found. Their exact interval comparison certifies the minimum
of q**2 * abs(pi - p/q) over all integers p coprime to q. No floating-point
arithmetic, continued fractions, or imports from the other experiments.

Pi is bounded by 4*(atan(1/2)+atan(1/3)), with termwise outward rounding
to a dyadic grid and an alternating-series remainder bound.
"""

from math import gcd
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json


def dyadic_arctan(inv, bits=1536, terms=900):
    scale = 1 << bits
    low = high = 0
    power = inv
    for k in range(terms):
        denominator = (2*k+1)*power
        numerator = scale if k % 2 == 0 else -scale
        low += numerator // denominator
        high += -((-numerator)//denominator)
        power *= inv*inv
    denominator = (2*terms+1)*power
    numerator = scale if terms % 2 == 0 else -scale
    low += min(0, numerator//denominator)
    high += max(0, -((-numerator)//denominator))
    return low, high


def pi_bounds(bits=1536, terms=900):
    low2, high2 = dyadic_arctan(2, bits, terms)
    low3, high3 = dyadic_arctan(3, bits, terms)
    return 4*(low2+low3), 4*(high2+high3), 1 << bits


def decimal_bounds(low, high, scale, digits=30):
    factor = 10**digits
    lower = low*factor//scale
    upper = -((-high*factor)//scale)
    def fmt(n):
        sign = '-' if n < 0 else ''
        n = abs(n)
        return f'{sign}{n//factor}.{n%factor:0{digits}d}'
    return [fmt(lower), fmt(upper)]


def all_denominators(bound):
    values = []
    power2, a = 1, 0
    while power2 <= bound:
        q, c = power2, 0
        while q <= bound:
            if q > 1:
                values.append((q,a,c))
            q *= 5
            c += 1
        power2 *= 2
        a += 1
    values.sort()
    assert len({q for q,a,c in values}) == len(values)
    return values


def scan(decimal_exponent=120):
    limit = 10**decimal_exponent
    low_pi, high_pi, scale = pi_bounds()
    assert 3*scale < low_pi < high_pi < 4*scale
    denominator_data = all_denominators(limit)
    digest = hashlib.sha256()
    history = []
    nearest_records = []
    counts = {'denominators': len(denominator_data),
              'coprime_side_candidates': 0,
              'nearest_coprime_comparisons': 0,
              'below_one': 0,
              'ambiguous_pi_floors': 0,
              'max_downward_coprime_steps': 0,
              'max_upward_coprime_steps': 0}
    best = None
    for q,a,c in denominator_data:
        integer_part = q*low_pi//scale
        if integer_part != q*high_pi//scale:
            counts['ambiguous_pi_floors'] += 1
            raise AssertionError('Insufficient pi precision for floor',q)
        p_down = integer_part
        while gcd(p_down,q) != 1:
            p_down -= 1
        p_up = integer_part+1
        while gcd(p_up,q) != 1:
            p_up += 1
        counts['max_downward_coprime_steps'] = max(counts['max_downward_coprime_steps'],integer_part-p_down)
        counts['max_upward_coprime_steps'] = max(counts['max_upward_coprime_steps'],p_up-integer_part-1)
        down = (q*(q*low_pi-p_down*scale),q*(q*high_pi-p_down*scale),p_down)
        up = (q*(p_up*scale-q*high_pi),q*(p_up*scale-q*low_pi),p_up)
        assert down[0] > 0 and up[0] > 0
        if down[1] < up[0]:
            lower,upper,p = down
        elif up[1] < down[0]:
            lower,upper,p = up
        else:
            raise AssertionError('Insufficient precision for nearest comparison',q)
        assert gcd(p,q)==1
        counts['coprime_side_candidates'] += 2
        counts['nearest_coprime_comparisons'] += 1
        if upper < scale:
            counts['below_one'] += 1
        else:
            assert lower > scale, ('Ambiguous threshold 1',p,q)
        digest.update(f'{q}:{a}:{c}:{p_down}:{p_up}:{p}\n'.encode())
        record = {'p':p,'q':q,'a':a,'c':c,'lower':lower,'upper':upper}
        nearest_records.append(record)
        if best is None or upper < best['lower']:
            best = record
            history.append({'p':p,'q':q,'powers':[a,c],
                            'quality_interval_decimal':decimal_bounds(lower,upper,scale)})
        elif lower <= best['upper']:
            raise AssertionError('Insufficient precision for minimum comparison',p,q)

    assert best is not None
    competitors = [r for r in nearest_records if r['q'] != best['q']]
    competitor_lower = min(r['lower'] for r in competitors)
    assert competitor_lower > best['upper']
    nearest_records.sort(key=lambda r:r['upper'])
    # Verify all ordered reported top records are separated, not sorted by
    # an unresolved interval overlap.
    for first,second in zip(nearest_records[:10],nearest_records[1:11]):
        assert first['upper'] < second['lower']
    return {
        'status':'PASS',
        'scope':'All noninteger terminating decimals with reduced denominator 1 < q <= 10^120; finite statement only',
        'arithmetic':'Exact integers and rational/dyadic bounds; no floats',
        'denominator_bound':str(limit),
        'denominator_definition':'q = 2^a * 5^c, a,c nonnegative integers; each q appears exactly once',
        'numerator_completeness':'Check nearest coprime integer below q*pi and nearest coprime integer above q*pi; all other coprime numerators are farther on their respective side',
        'quality_definition':'q^2 * abs(pi - p/q), with gcd(p,q)=1',
        'pi_certificate':{'identity':'pi = 4*(atan(1/2)+atan(1/3))',
                          'bits':1536,'terms_per_arctan':900,
                          'lower':str(F(low_pi,scale)),
                          'upper':str(F(high_pi,scale)),
                          'decimal':decimal_bounds(low_pi,high_pi,scale)},
        'counts':counts,
        'minimum':{'numerator':best['p'],'denominator':best['q'],
                   'denominator_powers_2_5':[best['a'],best['c']],
                   'quality_lower':str(F(best['lower'],scale)),
                   'quality_upper':str(F(best['upper'],scale)),
                   'quality_interval_decimal':decimal_bounds(best['lower'],best['upper'],scale)},
        'other_denominators_lower_bound':str(F(competitor_lower,scale)),
        'other_denominators_lower_bound_decimal':decimal_bounds(competitor_lower,competitor_lower,scale),
        'minimum_history':history,
        'top_10': [{'numerator':r['p'],'denominator':r['q'],
                    'quality_interval_decimal':decimal_bounds(r['lower'],r['upper'],scale)}
                   for r in nearest_records[:10]],
        'ordered_candidate_digest_sha256':digest.hexdigest(),
        'program_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'limitations':['No assertion for denominators above the bound',
                       'No inference about all real-world measurements or pi normality',
                       'Integers excluded; both 3/1 and 4/1 have quality below 1',
                       'Independent pi identity and denominator enumeration, but no external formal proof checker'],
    }


if __name__=='__main__':
    result = scan()
    output = Path(__file__).with_name('smooth_results.json')
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'counts':result['counts'],
                      'minimum':{k:v for k,v in result['minimum'].items() if 'lower' not in k and 'upper' not in k},
                      'output':str(output)},ensure_ascii=False))
