#!/usr/bin/env python3
"""Exact CF reconstruction and exploratory, chronological window comparison.

Requires gmpy2. This is an exploratory analysis of a record-selected prefix,
not a proof of predictive structure or a formal significance test for pi.
"""
from array import array
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import gmpy2 as g

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'snapshot' / 'records'


def exact_cf(p, q):
    result = []
    while q:
        a, r = divmod(p, q)
        result.append(int(a))
        p, q = q, r
    return result


def asserted_representation(record):
    p = g.mpz(record['p_hex'], 16)
    q = g.mpz(record['q_hex'], 16)
    assert q > 0 and g.gcd(p, q) == 1
    result = exact_cf(p, q)
    wanted = record['convergent_index_zero_based'] + 1
    canonical_length = len(result)
    if canonical_length + 1 == wanted:
        assert result[-1] > 1
        result[-1] -= 1
        result.append(1)
    assert len(result) == wanted, (canonical_length, wanted)
    return result, canonical_length


def previous_bin(a):
    if a == 1:
        return '1'
    if a == 2:
        return '2'
    if a <= 4:
        return '3-4'
    if a <= 9:
        return '5-9'
    return '10+'


def cell():
    return {'n': 0, 'hits': 0, 'expected': 0.0, 'variance': 0.0}


def update(item, hit, probability):
    item['n'] += 1
    item['hits'] += int(hit)
    item['expected'] += probability
    item['variance'] += probability * (1.0 - probability)


def finish(item):
    out = dict(item)
    out['observed_rate'] = item['hits'] / item['n'] if item['n'] else None
    out['expected_rate'] = item['expected'] / item['n'] if item['n'] else None
    out['observed_over_expected'] = item['hits'] / item['expected'] if item['expected'] else None
    out['descriptive_standardized_residual'] = (item['hits'] - item['expected']) / math.sqrt(item['variance']) if item['variance'] else None
    return out


def analysis(coefficients, threshold):
    # Largest appended a at the endpoint was selected as a record: exclude it.
    end_exclusive = len(coefficients) - 1
    split = end_exclusive // 2
    bins = ['1', '2', '3-4', '5-9', '10+']
    groups = {part: {'all': cell(), 'previous_a': {b: cell() for b in bins},
                     'previous_one_run': {str(k): cell() for k in range(4)}}
              for part in ['train', 'test']}
    beta = 0.0
    one_run = 0
    test_rows = []
    for j in range(1, end_exclusive):
        previous = coefficients[j - 1]
        if previous == 1:
            one_run += 1
        else:
            one_run = 0
        # beta = q_(j-2)/q_(j-1), before consuming a_j.
        probability = (1.0 + beta) / (threshold + beta)
        if j >= 64:
            part = 'train' if j < split else 'test'
            hit = coefficients[j] >= threshold
            b = previous_bin(previous)
            update(groups[part]['all'], hit, probability)
            update(groups[part]['previous_a'][b], hit, probability)
            update(groups[part]['previous_one_run'][str(min(3, one_run))], hit, probability)
            if part == 'test':
                test_rows.append((b, hit, probability))
        beta = 1.0 / (coefficients[j] + beta)
    # Predefined bins, shrink training ratios toward the universal CF model.
    # 20 expected pseudo-events per cell; no tuning on the held-out half.
    factors = {b: (groups['train']['previous_a'][b]['hits'] + 20.0) /
               (groups['train']['previous_a'][b]['expected'] + 20.0) for b in bins}
    global_factor = (groups['train']['all']['hits'] + 20.0) / (groups['train']['all']['expected'] + 20.0)
    baseline_loss = modified_loss = global_calibrated_loss = 0.0
    for b, y, probability in test_rows:
        modified = min(1.0 - 1e-12, probability * factors[b])
        baseline_loss -= math.log(probability) if y else math.log1p(-probability)
        modified_loss -= math.log(modified) if y else math.log1p(-modified)
        calibrated = min(1.0 - 1e-12, probability * global_factor)
        global_calibrated_loss -= math.log(calibrated) if y else math.log1p(-calibrated)
    result = {'threshold': threshold, 'minimum_target_index': 64,
              'train_indices': [64, split - 1],
              'test_indices': [split, end_exclusive - 1],
              'selected_endpoint_excluded': end_exclusive,
              'groups': {}, 'training_previous_a_multipliers': factors,
              'held_out_total_log_loss_baseline': baseline_loss,
              'held_out_total_log_loss_modified': modified_loss,
              'held_out_log_loss_improvement': baseline_loss - modified_loss,
              'training_global_multiplier': global_factor,
              'held_out_total_log_loss_global_calibration_only': global_calibrated_loss,
              'held_out_log_loss_improvement_over_global_calibration': global_calibrated_loss - modified_loss}
    for part, data in groups.items():
        result['groups'][part] = {'all': finish(data['all']),
            'previous_a': {k: finish(v) for k, v in data['previous_a'].items()},
            'previous_one_run': {k: finish(v) for k, v in data['previous_one_run'].items()}}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--records-dir', type=Path, default=SOURCE,
                        help='Directory of exact uploaded record JSON files')
    args = parser.parse_args()
    started = time.monotonic()
    records = [json.loads(path.read_text()) for path in sorted(args.records_dir.glob('*.json'))]
    if not records:
        raise SystemExit('No record JSON files found; set --records-dir to the snapshot records directory.')
    largest = max(records, key=lambda record: record['convergent_index_zero_based'])
    coefficients, canonical_length = asserted_representation(largest)
    coefficients.append(int(largest['following_partial_quotient']))
    integer_checks = []
    for record in records:
        c, canonical = asserted_representation(record)
        assert c == coefficients[:len(c)]
        assert int(record['following_partial_quotient']) == coefficients[len(c)]
        integer_checks.append({'index': record['convergent_index_zero_based'],
            'canonical_cf_length': canonical, 'asserted_cf_length': len(c),
            'matches_reconstructed_prefix': True,
            'following_a_matches_prefix': True,
            'last_a': c[-1]})
    data = array('I', coefficients)
    assert data.itemsize == 4
    if sys.byteorder != 'little':
        data.byteswap()
    raw = data.tobytes()
    with (HERE / 'prefix.u32le.gz').open('wb') as stream:
        stream.write(gzip.compress(raw, mtime=0))
    reconstruction = {'terms': len(coefficients), 'canonical_rational_length': canonical_length,
        'last_term': coefficients[-1], 'sha256_u32le': hashlib.sha256(raw).hexdigest(),
        'record_integer_checks': integer_checks,
        'method': 'Exact ordinary Euclidean algorithm on each coprime uploaded p/q; finite CF ambiguity resolved against stated length; last following term supplied by source and independently verified by parent.'}
    (HERE / 'reconstruction.json').write_text(json.dumps(reconstruction, indent=2))
    print('PREFIX_READY', len(coefficients), reconstruction['sha256_u32le'], flush=True)
    primary = analysis(coefficients, 1000)
    sensitivity = [analysis(coefficients, k) for k in [100, 10000]]
    out = {'reconstruction': reconstruction, 'primary': primary, 'sensitivity': sensitivity,
        'elapsed_seconds': time.monotonic() - started,
        'cautions': [
            'The available prefix ends at a record-selected location; it is not a random sample, even after omitting the final large digit.',
            'The last half was held out chronologically for the specified bin multipliers, not for all possible pattern searches.',
            'Residual z values are descriptive and are not calibrated discovery p-values.',
            'Conditional baseline is (1+beta)/(K+beta) for the uniform-real prior conditional on the observed finite CF prefix; applying it to fixed pi is a model comparison.',
            'The 10,184,980-term coverage claim cannot be recovered from these nine fractions alone.'
        ]}
    (HERE / 'window_results.json').write_text(json.dumps(out, indent=2))
    print(json.dumps({'seconds': out['elapsed_seconds'], 'threshold': 1000,
        'train': primary['groups']['train']['all'], 'test': primary['groups']['test']['all'],
        'held_out_log_loss_improvement': primary['held_out_log_loss_improvement'],
        'test_previous_a': primary['groups']['test']['previous_a']}, indent=2))


if __name__ == '__main__':
    main()
