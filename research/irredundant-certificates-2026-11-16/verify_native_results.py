#!/usr/bin/env python3
"""Post-run external mathematical verification, using only Python's standard library.

This is newly written review code, not code executed by the original native run.
It never imports or executes the uploaded research programs. Default fixtures are
sanitized native outputs. General theorems require their proofs, not these tests.
"""

import argparse
from collections import Counter, defaultdict
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def subset_sums(parts):
    bits = 1
    for part in parts:
        bits |= bits << part
    return bits


def covers(bits, first, width):
    if first < 0 or width < 0:
        return False
    mask = (1 << width) - 1
    return (bits >> first) & mask == mask


def interval_critical(parts, half):
    if sum(parts) != 2 * half or (subset_sums(parts) >> half) & 1:
        return False
    return all(covers(subset_sums(parts[:i] + parts[i + 1:]),
                      half - part + 1, part - 1)
               for i, part in enumerate(parts))


def literal_critical(parts, half):
    if (subset_sums(parts) >> half) & 1:
        return False
    for i, part in enumerate(parts):
        others = parts[:i] + parts[i + 1:]
        for left in range(1, part):
            split = others + (left, part - left)
            if not (subset_sums(split) >> half) & 1:
                return False
    return True


def partitions(total, lower=1):
    if total == 0:
        yield ()
    for first in range(lower, total + 1):
        for rest in partitions(total - first, first):
            yield (first,) + rest


def verify_theory(data, exhaustive_to):
    spectra = []
    all_patterns = 0
    for row in data['theory_spectra']:
        half = row['m']
        stored = {tuple(a) for a in row['critical_patterns']}
        require(len(stored) == len(row['critical_patterns']), 'duplicate pattern')
        require(sorted({2 * half - len(a) for a in stored}) == row['spectrum'],
                'stored spectrum inconsistent with patterns')
        for parts in stored:
            require(interval_critical(parts, half), 'stored pattern not critical')
            k = len(parts)
            require(2 * half <= k * (2 ** (k - 1) + 1), 'component bound failed')
            all_patterns += 1
        if half <= exhaustive_to:
            candidates = list(partitions(2 * half))
            require(len(candidates) == row['partitions_examined'], 'partition count')
            actual = {a for a in candidates if literal_critical(a, half)}
            require(actual == stored, 'complete critical-pattern set mismatch')
            spectra.append({'m': half, 'spectrum': row['spectrum'],
                            'partitions_examined': len(candidates)})
    require(all_patterns == 349, 'expected 349 stored critical patterns')
    family_count = literal_count = 0
    parameter_pairs = set()
    for row in data['near_double_cases']:
        r, stage = row['r'], row['s']
        parameter_pairs.add((r, stage))
        parts = ((1,) * (r - 1) + (r + 1,) +
                 tuple(b for j in range(1, stage + 1) for b in [2 ** j * r + 1] * 2))
        half = (2 ** (stage + 1) - 1) * r + stage
        require(sorted(parts) == row['parts'] and half == row['m'], 'family formula')
        require(len(parts) == r + 2 * stage == row['component_count'], 'family k')
        require(2 * half - len(parts) == row['cost'], 'family cost')
        require(interval_critical(parts, half), 'family interval criticality')
        width = 2 ** stage * r
        bits = subset_sums(parts)
        require(covers(bits, half - width, width) and
                covers(bits, half + 1, width), 'near-half coverage')
        if half <= 512:
            require(literal_critical(parts, half), 'family literal criticality')
            literal_count += 1
        family_count += 1
    require(parameter_pairs == {(r, s) for r in range(1, 13)
                                for s in range(1, r + 2)} and family_count == 90,
            'family parameter coverage')
    return {'complete_literal_partition_checks': spectra,
            'all_stored_patterns_interval_and_bound_checked': all_patterns,
            'family_interval_and_coverage_cases': family_count,
            'family_literal_cases_half_at_most_512': literal_count,
            'limit': 'm above exhaustive_to: stored patterns checked for validity; '
                     'completeness of their lists is not independently re-enumerated'}


def inverse(word):
    return tuple(-x for x in reversed(word))


def reduce_word(word):
    stack = []
    for letter in word:
        if stack and stack[-1] + letter == 0:
            stack.pop()
        else:
            stack.append(letter)
    return tuple(stack)


def move(state, action):
    a, b = state
    if action == 0:
        return inverse(a), b
    if action == 1:
        return a, inverse(b)
    if action == 2:
        return reduce_word(a + b), b
    if action == 3:
        return reduce_word(a + inverse(b)), b
    if action == 4:
        return a, reduce_word(b + a)
    if action == 5:
        return a, reduce_word(b + inverse(a))
    index = (action - 6) // 4
    generator = (1, -1, 2, -2)[(action - 6) % 4]
    changed = reduce_word((generator,) + state[index] + (-generator,))
    return (changed, b) if index == 0 else (a, changed)


INVERSE_IDS = (0, 1, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12)


def ball(root):
    distance = {root: 0}
    parents = defaultdict(list)
    layers = [[root]]
    evaluations = 0
    for depth in range(4):
        following = []
        for state in layers[-1]:
            for action in range(14):
                target = move(state, action)
                require(move(target, INVERSE_IDS[action]) == state, 'inverse action')
                evaluations += 1
                if target not in distance:
                    distance[target] = depth + 1
                    following.append(target)
                if distance[target] == depth + 1:
                    parents[target].append((state, action))
        layers.append(following)
    return distance, parents, list(map(len, layers)), evaluations


def parent_paths(parents):
    @lru_cache(None)
    def visit(state):
        if not parents[state]:
            return ((),)
        return tuple(prefix + (action,) for previous, action in parents[state]
                     for prefix in visit(previous))
    return visit


def state_digest(states):
    return sha256(json.dumps(sorted(states), separators=(',', ':')).encode('ascii'))


def token_cost(path, macro):
    costs = [0] + [99] * len(path)
    for i in range(len(path)):
        costs[i + 1] = min(costs[i + 1], costs[i] + 1)
        end = i + len(macro)
        if path[i:end] == macro:
            costs[end] = min(costs[end], costs[i] + 1)
    return costs[-1]


def signature(lengths):
    return len(lengths) - 1, max(lengths), sum(lengths[1:])


def verify_ac(data):
    start = ((-2, -1), (-2, -2, -1, -1, 2, 2, 1))
    target = ((1,), (2,))
    front, fp, fl, fe = ball(start)
    back, bp, bl, be = ball(target)
    f3 = {state for state in front if front[state] <= 3}
    common = front.keys() & back.keys()
    require(not f3 & back.keys() and len(common) == 2, 'ball intersections')
    expectation = data['ac_ball_expectations']
    require(fl == expectation['forward_layer_sizes'] and
            bl == expectation['target_layer_sizes'], 'B4 layer counts')
    digests = {'forward3_sha256': state_digest(f3),
               'forward4_sha256': state_digest(front),
               'target4_sha256': state_digest(back)}
    require(all(digests[k] == expectation[k] for k in digests), 'ball digests')
    prefix, suffix = parent_paths(fp), parent_paths(bp)
    paths = sorted(p + tuple(INVERSE_IDS[i] for i in reversed(q))
                   for state in common for p in prefix(state) for q in suffix(state))
    require(len(paths) == len(set(paths)) == 8, 'all eight labelled paths')
    records, histogram = [], Counter()
    signature_cuts = 0
    for path in paths:
        state, lengths = start, [9]
        for action in path:
            state = move(state, action)
            lengths.append(sum(map(len, state)))
        require(state == target, 'path endpoint')
        peak, work = max(lengths), sum(lengths)
        histogram[peak, work] += 1
        records.append({'moves': list(path), 'lengths': lengths,
                        'peak': peak, 'work': work})
        for cut in range(9):
            left, right = signature(lengths[:cut + 1]), signature(lengths[cut:])
            combined = (left[0] + right[0], max(left[1], right[1]), left[2] + right[2])
            require(combined == signature(lengths), 'resource composition')
            require(sum(lengths[:cut + 1]) + sum(lengths[cut:]) - lengths[cut] == work,
                    'initial-inclusive work composition')
            signature_cuts += 1
    require(records == data['ac_paths'], 'native paths or resource trajectories')
    require(histogram == {(9, 48): 2, (9, 49): 4, (9, 51): 2}, 'resource distribution')
    macros = {p[i:i + k] for p in paths for k in (2, 3, 4) for i in range(9 - k)}
    costs = [(len(m) + sum(token_cost(p, m) for p in paths), m) for m in macros]
    best = min([64] + [cost for cost, _ in costs])
    winners = sorted(m for cost, m in costs if cost == best)
    singles = [min([8] + [len(m) + token_cost(p, m) for m in macros]) for p in paths]
    require(len(macros) == 55 and best == 55 and winners == [(2, 5, 5)] and
            singles == [8] * 8, 'finite flat-macro optimum')
    return {'forward_layer_sizes': fl, 'target_layer_sizes': bl,
            'generated_edges': fe + be, 'ball_digests': digests,
            'intersection_B3_B4': 0, 'intersection_B4_B4': len(common),
            'all_length8_paths': len(paths), 'exact_distance_in_frozen_table': 8,
            'resource_distribution': [{'peak': p, 'work': w, 'count': n}
                                      for (p, w), n in sorted(histogram.items())],
            'composition_cut_checks': signature_cuts,
            'flat_macro': {'candidate_count': len(macros), 'baseline': 64,
                           'best_paid_cost': best, 'optimal_macros': winners,
                           'single_path_paid_optima': singles},
            'official_core_or_cli_run': False, 'platform_upload': False}


def main():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=base / 'data/native-results.json')
    parser.add_argument('--output', type=Path, default=base / 'external-verification.json')
    parser.add_argument('--exhaustive-to', type=int, default=12, choices=range(3, 21),
                        help='literal all-partition checking; default 12, optional through 20')
    args = parser.parse_args()
    raw = args.data.read_bytes()
    data = json.loads(raw)
    result = {'format': 'native465-external-review/1',
              'verification_stage': 'external post-native-run review',
              'original_native_receipt': False,
              'supplied_research_code_executed': False,
              'checker_sha256': sha256(Path(__file__).read_bytes()),
              'fixture_sha256': sha256(raw), 'exhaustive_to': args.exhaustive_to}
    try:
        result['theory'] = verify_theory(data, args.exhaustive_to)
        result['application'] = verify_ac(data)
        result['status'] = 'passed'
    except Exception as error:
        result['status'] = 'failed'
        result['error'] = str(error)
        args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
        raise
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'status': result['status'], 'stage': result['verification_stage'],
                      'theory_patterns_checked': result['theory'][
                          'all_stored_patterns_interval_and_bound_checked'],
                      'family_cases': result['theory']['family_interval_and_coverage_cases'],
                      'ac_paths': result['application']['all_length8_paths']}, indent=2))


if __name__ == '__main__':
    main()
