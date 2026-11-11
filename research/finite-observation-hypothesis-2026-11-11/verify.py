#!/usr/bin/env python3
"""Independent audit of the finite rounding-gcd certificates.

This file deliberately does not import the experiment generator. Its pi
interval uses a different arctangent identity, and its continued fractions
come from ordinary Euclid algorithms applied separately to both endpoints.
"""
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
import json

BITS = 11520
UNIT = 1 << BITS


def reciprocal_angle_interval(den):
    lower = upper = 0
    power = den
    i = 0
    while True:
        divisor = (2*i+1)*power
        if divisor > UNIT * (1 << 32):
            break
        sign = 1 if i % 2 == 0 else -1
        numerator = sign*UNIT
        lower += numerator // divisor
        upper -= (-numerator) // divisor
        i += 1
        power *= den*den
    # Alternation bounds the whole omitted tail by its first term.
    numerator = (1 if i % 2 == 0 else -1)*UNIT
    lower += min(0, numerator // divisor)
    upper += max(0, -((-numerator)//divisor))
    return lower, upper, i


def independent_pi():
    a, b, n = reciprocal_angle_interval(2)
    c, d, m = reciprocal_angle_interval(3)
    return 4*(a+c), 4*(b+d), (n, m)


def euclid_expansion(a, b):
    result = []
    while b:
        digit, remainder = divmod(a, b)
        result.append(digit)
        a, b = b, remainder
    return result


def common_convergents(lo, hi, denominator_limit):
    left = euclid_expansion(lo, UNIT)
    right = euclid_expansion(hi, UNIT)
    prefix = []
    for a, b in zip(left, right):
        if a != b:
            break
        prefix.append(a)
    p0, p1, q0, q1 = 0, 1, 1, 0
    rows = []
    for index, a in enumerate(prefix):
        p, q = a*p1+p0, a*q1+q0
        if q > denominator_limit:
            return prefix, rows, (index, p, q)
        assert gcd(p, q) == 1
        rows.append((index, p, q))
        p0, p1, q0, q1 = p1, p, q1, q
    raise AssertionError("Independent interval cannot reach denominator cutoff")


def abs_linear_interval(p, q, lo, hi):
    a, b = q*lo-p*UNIT, q*hi-p*UNIT
    if a > 0:
        return a, b
    assert b < 0
    return -b, -a


def full_window(p, q, lo, hi):
    # k*|q*pi-p|<1/2, computed directly as a linear-form interval.
    e0, e1 = abs_linear_interval(p, q, lo, hi)
    k0 = UNIT // (2*e1)
    k1 = UNIT // (2*e0)
    assert k0 == k1, "Uncertified reciprocal threshold"
    assert 2*k0*e1 < UNIT
    assert 2*(k0+1)*e0 > UNIT
    return k0, (q*e0, q*e1)


def rounded_reading(B, lo, hi):
    a = (2*B*lo+UNIT)//(2*UNIT)
    b = (2*B*hi+UNIT)//(2*UNIT)
    assert a == b
    return a


def best_for_height(H, centers):
    # B=1, A=3 supplies baseline one; any larger value must have a CF center.
    best = (Fraction(1), 0, 3, 1, 1)
    for index, p, q, K in centers:
        if q*q > H:
            break
        k = min(H//q, K)
        value = Fraction(k, q)
        if value > best[0]:
            best = (value, index, p, q, k)
    return best


def fraction_inside_decimal(value, bounds):
    assert Fraction(bounds[0]) <= value <= Fraction(bounds[1])


def interval_inside_decimal(lo, hi, bounds):
    assert Fraction(bounds[0]) <= Fraction(lo, UNIT)
    assert Fraction(hi, UNIT) <= Fraction(bounds[1])


def main():
    folder = Path(__file__).parent
    source = json.loads((folder/'results.json').read_text())
    assert source['status'] == 'PASS'
    assert source['script_sha256'] == sha256((folder/'experiment.py').read_bytes()).hexdigest()
    H = int(source['resolution_limit'])
    root_H = isqrt(H)
    assert int(source['convergent_denominator_limit']) == root_H
    lo, hi, terms = independent_pi()
    original = source['pi_interval']
    old_unit = 1 << original['bits']
    old_lo, old_hi = int(original['lower_numerator']), int(original['upper_numerator'])
    assert old_lo*UNIT < lo*old_unit < hi*old_unit < old_hi*UNIT
    assert 3*UNIT < lo < hi < 4*UNIT

    prefix, raw_centers, next_center = common_convergents(lo, hi, root_H)
    claimed_digits = source['certified_cf_digits']
    assert prefix[:len(claimed_digits)] == claimed_digits
    assert len(claimed_digits) == len(raw_centers)+1
    assert source['convergents_within_root_limit'] == len(raw_centers)
    assert int(source['first_convergent_above_root_limit']['p']) == next_center[1]
    assert int(source['first_convergent_above_root_limit']['q']) == next_center[2]
    assert len(source['centers']) == len(raw_centers)
    centers = []
    kappas = {}
    for row, stated in zip(raw_centers, source['centers']):
        index, p, q = row
        assert (stated['index'], int(stated['p']), int(stated['q'])) == row
        K, quality = full_window(p, q, lo, hi)
        assert int(stated['full_rounding_multiplier_K']) == K
        interval_inside_decimal(*quality, stated['quality_decimal_interval'])
        side = 'below' if q*lo > p*UNIT else 'above'
        assert side == stated['side_of_fraction']
        if K:
            assert rounded_reading(K*q, lo, hi) == K*p
            assert gcd(K*p, K*q) == K
        assert rounded_reading((K+1)*q, lo, hi) != (K+1)*p
        centers.append((index, p, q, K))
        kappas[index] = quality

    def verify_winner(claim, height):
        value, index, p, q, k = best_for_height(height, centers)
        assert int(claim['H']) == height
        assert Fraction(*map(int, claim['G_ratio'])) == value
        assert int(claim['B']) == k*q
        assert int(claim['A']) == k*p
        assert int(claim['gcd']) == k == gcd(k*p, k*q)
        assert rounded_reading(k*q, lo, hi) == k*p
        if value > 1:
            assert claim['index'] == index
            assert int(claim['p']) == p and int(claim['q']) == q
            assert int(claim['full_K']) == centers[index][3]
            assert int(claim['capped_K']) == k
            fraction_inside_decimal(value, claim['G_decimal_interval'])
            interval_inside_decimal(*kappas[index], claim['quality_decimal_interval'])
        return value, index, p, q, k

    global_result = verify_winner(source['global_maximum'], H)
    assert global_result[0] > 1

    expected_changes = []
    previous = None
    for exponent in range(1, 1001):
        height = 10**exponent
        result = best_for_height(height, centers)
        key = result[0], result[1]
        if key != previous:
            expected_changes.append((exponent, height))
            previous = key
    assert len(expected_changes) == len(source['power_of_ten_changes'])
    for (exponent, height), claim in zip(expected_changes, source['power_of_ten_changes']):
        assert claim['decimal_exponent'] == exponent
        verify_winner(claim, height)

    for claim in source['examples']:
        index = claim['index']
        _, p, q, K = centers[index]
        assert (int(claim['p']), int(claim['q']), int(claim['full_K'])) == (p, q, K)
        assert int(claim['last_resolution']) == K*q
        at = rounded_reading(K*q, lo, hi)
        after = rounded_reading((K+1)*q, lo, hi)
        assert int(claim['last_readout']) == at == K*p
        assert int(claim['last_gcd']) == gcd(at, K*q) == K
        assert int(claim['next_multiple_resolution']) == (K+1)*q
        assert int(claim['next_actual_readout']) == after
        assert int(claim['next_center_numerator']) == (K+1)*p != after
        assert claim['next_multiple_exits_center'] is True
        fraction_inside_decimal(Fraction(K, q), claim['G_at_last_multiple'])
        interval_inside_decimal(*kappas[index], claim['quality_decimal_interval'])

    # Direct enumeration of a separate smaller range. Every record and every
    # prefix maximum is checked against the globally complete CF classifier.
    direct_limit = 10000
    direct_stream = sha256()
    direct_max = Fraction(0)
    direct_records = []
    cf_index = {(p,q): (i,K) for i,p,q,K in centers}
    above_one = 0
    for B in range(1, direct_limit+1):
        A = rounded_reading(B, lo, hi)
        g = gcd(A, B)
        value = Fraction(g*g, B)
        if value > 1:
            above_one += 1
            p, q = A//g, B//g
            index, full_K = cf_index[p,q]
            assert q*q < B and g <= full_K
            assert 2*kappas[index][1]*g < q*UNIT
        if value > direct_max:
            direct_max = value
            direct_records.append({'B':B,'A':A,'gcd':g,
                                   'G_ratio':[str(value.numerator),str(value.denominator)]})
        assert best_for_height(B, centers)[0] == direct_max
        direct_stream.update(f'{B},{A},{g}\n'.encode())
    claimed_records = [r for r in source['direct_crosscheck']['records'] if r['B']<=direct_limit]
    assert claimed_records == direct_records

    # Generate LCM clocks independently by multiplying at prime-power events,
    # rather than using the generator's successive gcd/lcm updates.
    smallest = list(range(1001))
    for p in range(2, isqrt(1000)+1):
        if smallest[p] == p:
            for j in range(p*p, 1001, p):
                if smallest[j] == j:
                    smallest[j] = p
    D = 1
    lcm_hash = sha256()
    lcm_changes = []
    last_index = None
    divisor_candidates = 0
    for n in range(1,1001):
        if n >= 2:
            p, rest = smallest[n], n
            while rest % p == 0:
                rest //= p
            if rest == 1:
                D *= p
        eligible = [i for i,p,q,K in centers if q<=D and D%q==0]
        divisor_candidates += len(eligible)
        best_i = min(eligible, key=lambda i:sum(kappas[i]))
        # The midpoint ordering is accepted only after proving disjointness
        # of the winning interval and every competing interval.
        for i in eligible:
            if i != best_i:
                assert kappas[best_i][1] < kappas[i][0]
        _, p, q, _ = centers[best_i]
        lcm_hash.update(f'{n},{D},{best_i},{p},{q}\n'.encode())
        if best_i != last_index:
            lcm_changes.append((n,D,best_i,p,q))
            last_index = best_i
    assert all(D % n == 0 for n in range(1,1001))
    lcm_source = source['lcm_nested_grid']
    assert divisor_candidates == lcm_source['divisibility_candidates']
    assert lcm_hash.hexdigest() == lcm_source['all_layers_sha256']
    assert len(lcm_changes) == len(lcm_source['record_changes'])
    for actual, claim in zip(lcm_changes, lcm_source['record_changes']):
        n,d,i,p,q = actual
        assert (claim['first_n'],int(claim['D_n']),claim['index'],int(claim['p']),int(claim['q'])) == actual
        assert int(claim['grid_numerator']) == p*(d//q)
        assert int(claim['grid_gcd']) == gcd(p*(d//q),d) == d//q
        interval_inside_decimal(*kappas[i],claim['quality_decimal_interval'])
    final_lcm = lcm_source['final']
    assert final_lcm['n'] == 1000 and int(final_lcm['D_n']) == D
    assert final_lcm['D_n_decimal_digits'] == len(str(D))
    assert (final_lcm['index'],int(final_lcm['p']),int(final_lcm['q'])) == (best_i,p,q)
    interval_inside_decimal(*kappas[best_i],final_lcm['quality_decimal_interval'])

    report = {
        'status':'PASS',
        'scope':'Independent exact finite maximum audit; no infinite pi conclusion',
        'independent_pi':{
            'method':'4*(atan(1/2)+atan(1/3)), directed integer term bounds and alternating tails',
            'bits':BITS,'terms':list(terms),'lower_numerator':str(lo),
            'upper_numerator':str(hi),'width_numerator':str(hi-lo),
            'strictly_inside_generator_interval':True},
        'continued_fractions':{
            'method':'separate ordinary Euclid expansions of rational endpoints',
            'common_prefix_digits':len(prefix),'required_digits':len(claimed_digits),
            'all_required_digits_match':True,'centers_verified':len(centers)},
        'rounding_windows_verified':len(centers),
        'window_exit_checks':len(centers),
        'global_maximum':{
            'resolution_limit':str(H),'winning_center_index':global_result[1],
            'winning_denominator_digits':len(str(global_result[3])),
            'G_exact_ratio':[str(global_result[0].numerator),str(global_result[0].denominator)],
            'G_decimal_interval':source['global_maximum']['G_decimal_interval'],
            'complete_above_baseline_one':True},
        'power_of_ten_envelopes_verified':1000,
        'power_of_ten_change_records_verified':len(expected_changes),
        'examples_verified':len(source['examples']),
        'lcm_nested_grid':{'clock_construction':'independent prime-power events',
                           'layers_verified':1000,
                           'record_layers':[r[0] for r in lcm_changes],
                           'final_clock_decimal_digits':len(str(D)),
                           'divisibility_candidates':divisor_candidates,
                           'all_layers_sha256':lcm_hash.hexdigest()},
        'direct_scan':{'resolution_range':[1,direct_limit],
                       'prefix_maxima_verified':direct_limit,
                       'G_above_one_count':above_one,
                       'record_count':len(direct_records),
                       'records_match_generator_prefix':True,
                       'sha256':direct_stream.hexdigest()},
        'inputs_sha256':{name:sha256((folder/name).read_bytes()).hexdigest()
                         for name in ('experiment.py','results.json')},
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    (folder/'verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','centers':len(centers),
                      'global_G':source['global_maximum']['G_decimal_interval'],
                      'power_changes':len(expected_changes),'direct_resolutions':direct_limit}))


if __name__ == '__main__':
    main()

