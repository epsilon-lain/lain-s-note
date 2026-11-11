#!/usr/bin/env python3
"""Exact finite checks for rounded-readout gcd properties.

Source center windows are reused as certified data, never imported as code.
Finite counts and frequencies do not determine the infinite pi question.
"""
from collections import Counter
from hashlib import sha256
from math import gcd
from pathlib import Path
import json

BITS = 2048
SCALE = 1 << BITS
DIRECT_LIMIT = 10**6
THRESHOLDS = (1, 10, 100, 1000)
EXPONENTS = list(range(1, 11)) + [20, 50, 100, 200, 300, 400, 500, 1000]
MEAN_BITS = 128


def atan_interval(inv):
    lo = hi = k = 0
    power = inv
    while True:
        den = (2*k+1)*power
        if den.bit_length() > BITS+40:
            break
        v = SCALE if k % 2 == 0 else -SCALE
        lo += v//den
        hi += -((-v)//den)
        power *= inv*inv
        k += 1
    v = SCALE if k % 2 == 0 else -SCALE
    lo += min(0, v//den)
    hi += max(0, -((-v)//den))
    return lo, hi, k


def pi_interval():
    l5,u5,n5 = atan_interval(5)
    l239,u239,n239 = atan_interval(239)
    return 16*l5-4*u239,16*u5-4*l239,[n5,n239]


def readout(B,lo,hi):
    a = (2*B*lo+SCALE)//(2*SCALE)
    assert a == (2*B*hi+SCALE)//(2*SCALE)
    g = gcd(a,B)
    e_lo,e_hi = B*lo-a*SCALE,B*hi-a*SCALE
    quartile = (4*e_lo+2*SCALE)//SCALE
    assert quartile == (4*e_hi+2*SCALE)//SCALE
    assert 0 <= quartile < 4
    return a,g,quartile


def decimal_bounds(lo,hi,den,places=24):
    s = 10**places
    a,b = lo*s//den,-((-hi*s)//den)
    def display(n):
        return str(n//s)+'.'+str(n%s).zfill(places)
    return [display(a),display(b)]


def load_source():
    parent = Path(__file__).resolve().parent.parent/'finite-observation-hypothesis-2026-11-11'
    data = (parent/'results.json').read_bytes()
    manifest = json.loads((parent/'manifest.json').read_text())
    digest = sha256(data).hexdigest()
    assert digest == manifest['files']['results.json']['sha256']
    source = json.loads(data)
    assert source['status'] == 'PASS'
    assert source['resolution_limit'] == str(10**1000)
    centers = [(int(r['p']),int(r['q']),int(r['full_rounding_multiplier_K']))
               for r in source['centers']]
    assert len(centers) == 969
    return centers,digest,sha256((parent/'manifest.json').read_bytes()).hexdigest()


def count_high(centers,N,T):
    return sum(max(0,min(K,N//q)-T*q) for p,q,K in centers if q*q < N//T+1)


def main():
    lo,hi,terms = pi_interval()
    centers,source_hash,source_manifest_hash = load_source()
    expected = {T:set() for T in THRESHOLDS}
    for p,q,K in centers:
        for T in THRESHOLDS:
            if T*q*q >= DIRECT_LIMIT:
                continue
            expected[T].update(k*q for k in range(T*q+1,min(K,DIRECT_LIMIT//q)+1))

    hist = Counter()
    divisibility = [0]*21
    bins = {str(d):[0]*4 for d in range(1,6)}
    bins['other'] = [0]*4
    actual = {T:[] for T in THRESHOLDS}
    event_hash = {T:sha256() for T in THRESHOLDS}
    read_hash,adj_hash = sha256(),sha256()
    mean_g_sum = mean_recip_lo = 0
    snapshots = {}
    adjacent_checks = 0
    previous = None

    def event_segments(N,T):
        answer = []
        for index,(p,q,K) in enumerate(centers):
            first,last = T*q+1,min(K,N//q)
            if first <= last:
                answer.append({'center_index':index,'q':str(q),
                               'k_first':str(first),'k_last':str(last)})
        return answer

    for B in range(1,DIRECT_LIMIT+1):
        A,g,quartile = readout(B,lo,hi)
        hist[g] += 1
        mean_g_sum += g
        mean_recip_lo += (1 << MEAN_BITS)//g
        for d in range(1,21):
            divisibility[d] += int(g%d == 0)
        bins[str(g) if g <= 5 else 'other'][quartile] += 1
        for T in THRESHOLDS:
            is_high = g*g > T*B
            assert is_high == (B in expected[T])
            if is_high:
                actual[T].append(B)
                event_hash[T].update(f'{B}\n'.encode())
        read_hash.update(f'{B},{A},{g},{quartile}\n'.encode())
        if previous and B >= 5:
            prev_B,prev_A,prev_g = previous
            m = prev_A-3*prev_B
            w = (A-3*B)-m
            assert w in (0,1)
            D = m if w == 0 else prev_B-m
            assert 0 < D < prev_B
            assert D%(prev_g*g) == 0
            assert 4*(prev_g*g)**2 < (prev_B+B)**2
            assert (prev_g*g)**2 <= D*D < prev_B*B
            adjacent_checks += 1
            adj_hash.update(f'{prev_B},{prev_g},{g},{m},{w},{D}\n'.encode())
        previous = B,A,g
        if B in (100000,DIRECT_LIMIT):
            assert all(len(actual[T]) == count_high(centers,B,T) for T in THRESHOLDS)
            den = B*(1 << MEAN_BITS)
            snapshots[str(B)] = {
                'count':B,'g_histogram_1_to_20':{str(d):hist[d] for d in range(1,21)},
                'g_above_20_count':sum(v for d,v in hist.items() if d>20),
                'divisibility_counts_1_to_20':{str(d):divisibility[d] for d in range(1,21)},
                'residual_quartiles_by_g':{k:list(v) for k,v in bins.items()},
                'g_sum':str(mean_g_sum),
                'mean_reciprocal_g_dyadic_bound':{
                    'bits':MEAN_BITS,'lower_sum_numerator':str(mean_recip_lo),
                    'upper_sum_numerator':str(mean_recip_lo+B),'denominator':str(den),
                    'decimal_interval':decimal_bounds(mean_recip_lo,mean_recip_lo+B,den)},
                'readout_stream_sha256':read_hash.hexdigest(),
                'adjacent_stream_sha256':adj_hash.hexdigest(),
                'adjacent_checks':adjacent_checks,
                'high_G_by_threshold':{
                    str(T):{'count':len(actual[T]),'bases_first_20':actual[T][:20],
                            'bases_last_20':actual[T][-20:],
                            'arithmetic_segments':event_segments(B,T),
                            'sha256':event_hash[T].hexdigest()} for T in THRESHOLDS}}

    large_counts = []
    for n in EXPONENTS:
        N = 10**n
        large_counts.append({'decimal_exponent':n,'N':str(N),
                             'counts':{str(T):str(count_high(centers,N,T)) for T in THRESHOLDS},
                             'squared_normalized_counts':{
                                 str(T):[str(count_high(centers,N,T)**2),str(N)]
                                 for T in THRESHOLDS}})

    checkpoints = []
    for index,(p,q,K) in enumerate(centers):
        cap = min(K,10**1000//q)
        for k in sorted({q+1,cap}):
            if q < k <= cap:
                checkpoints.append((k*q,k,index,p,q))
    checkpoints.sort()
    pair_hash = sha256()
    pair_checks = separated_checks = 0
    samples = []
    pairs = set()
    for i in range(len(checkpoints)):
        for gap in (1,7,31):
            if i+gap < len(checkpoints):
                pairs.add((i,i+gap))
    for i,j in sorted(pairs):
        B,g,index,p,q = checkpoints[i]
        C,h,other,r,s = checkpoints[j]
        if index == other:
            continue
        assert B < C
        D = abs((g*p)*C-(h*r)*B)
        assert D > 0 and D%(g*h) == 0
        assert g*h <= D and 2*D < B+C
        assert 4*g*g*h*h < (B+C)**2
        for T in THRESHOLDS:
            if g*g > T*B and h*h > T*C:
                assert (B+C)**2 > 4*T*T*B*C
                separated_checks += 1
        pair_hash.update(f'{i},{j},{B},{C},{g},{h},{index},{other}\n'.encode())
        pair_checks += 1
        if len(samples) < 5 and len(str(C)) <= 20:
            samples.append({'B':str(B),'C':str(C),'g_B':str(g),'g_C':str(h),
                            'center_indices':[index,other],'determinant':str(D)})

    examples = []
    for B in (1874217,1874218,1874219):
        A,g,quartile = readout(B,lo,hi)
        examples.append({'B':B,'A':A,'gcd':g,'reduced_p':A//g,'reduced_q':B//g,
                         'G_ratio':[str(g*g),str(B)],
                         'G_decimal_interval':decimal_bounds(g*g,g*g,B)})
    result = {
        'status':'PASS','scope':'Finite exact property checks and high-G counts; no infinite pi bound',
        'definitions':{'N_B':'round(B*pi)','g_B':'gcd(N_B,B)','G':'g_B^2/B',
                       'residual_quartiles':'floor(4*(B*pi-N_B+1/2)); bins 0,1,2,3',
                       'count_formula':'sum_centers max(0,min(fullK,N//q)-T*q), T integer >=1'},
        'reused_certificate':{'directory':'../finite-observation-hypothesis-2026-11-11',
                              'results_sha256':source_hash,'manifest_sha256':source_manifest_hash,
                              'center_count':len(centers),'N_max':str(10**1000)},
        'pi_interval':{'method':'16 atan(1/5)-4 atan(1/239), integer directed dyadic terms',
                       'bits':BITS,'terms':terms,'lower_numerator':str(lo),'upper_numerator':str(hi)},
        'direct_scan':{'resolutions':[1,DIRECT_LIMIT],
                       'every_high_event_matches_complete_classifier':True,'snapshots':snapshots},
        'large_counts':large_counts,
        'different_center_pair_checks':{'checkpoint_count':len(checkpoints),'pair_count':pair_checks,
                                        'threshold_separation_checks':separated_checks,
                                        'sha256':pair_hash.hexdigest(),'small_samples':samples},
        'neighboring_peak_example':examples,
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    target = Path(__file__).with_name('results.json')
    target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','direct_count':DIRECT_LIMIT,
                      'high_counts':{T:len(actual[T]) for T in THRESHOLDS},
                      'adjacent_checks':adjacent_checks,'pair_checks':pair_checks,
                      'bytes':target.stat().st_size}))


if __name__ == '__main__':
    main()
