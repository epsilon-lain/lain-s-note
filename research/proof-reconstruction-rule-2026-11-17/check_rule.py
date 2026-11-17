#!/usr/bin/env python3
"""Check a classical donor-conjugation identity against pinned AC semantics.

This checks synthetic rule instances and a separate frozen development scan,
not a newly solved competition instance.
It makes no native-host calls and does not establish a novel mathematical rule.
Only checks.json next to this script is written; bytecode output is disabled.
"""

import argparse
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import random
import sys
import types

sys.dont_write_bytecode = True

ROOT = Path('/workspace/scratch/7a6c9da7b31a/ac_official')
VERIFIER = ROOT / 'competition/tools/verifier'
EXPECTED_CORE_BLOB = '0a5da5fe6e3e42c96877734dad46b155a69c02b2'
PINNED_COMMIT = 'a0fd6e6f52d82c93ccc06ab91d81e8fa3678256e'
LETTERS = (1, -1, 2, -2)
CONJ_MOVE = {1: 6, -1: 7, 2: 8, -2: 9}
COMPILED = {1: [0, 3, 0, 3], -1: [0, 2, 0, 2]}
COUNTS = {'official_apply_move_calls': 0, 'official_verify_calls': 0,
          'route_checks': 0}


def git_blob(path):
    data = path.read_bytes()
    return hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()


def reduced(word):
    """Independent expected-value oracle: stack reduction in F(x,y)."""
    stack = []
    for letter in word:
        if stack and stack[-1] == -letter:
            stack.pop()
        else:
            stack.append(letter)
    return tuple(stack)


def inverse(word):
    return tuple(-letter for letter in reversed(word))


def load_official_core():
    # Refuse to execute the cached authoritative core before checking its blob.
    actual = git_blob(VERIFIER / 'core.py')
    if actual != EXPECTED_CORE_BLOB:
        raise RuntimeError(f'Pinned core Git blob mismatch: {actual}')
    # A minimal package loads only canon/core, not verifier.__init__ re-exports.
    package_name = '_pinned_donor_rule_verifier'
    package = types.ModuleType(package_name)
    package.__path__ = [str(VERIFIER)]
    sys.modules[package_name] = package
    loaded = {}
    for name in ('canon', 'core'):
        spec = importlib.util.spec_from_file_location(
            f'{package_name}.{name}', VERIFIER / f'{name}.py')
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        setattr(package, name, module)
        spec.loader.exec_module(module)
        loaded[name] = module
    core = loaded['core']
    original_apply = core.apply_move

    def counted_apply(state, move):
        COUNTS['official_apply_move_calls'] += 1
        return original_apply(state, move)

    core.apply_move = counted_apply
    return core


def expected_endpoint(u, v, sign):
    conjugator = v if sign == 1 else inverse(v)
    return (reduced(conjugator + u + inverse(conjugator)), v)


def literal_moves(v, sign):
    conjugator = v if sign == 1 else inverse(v)
    # Successive C_a compose in reverse execution order: execute reverse(w).
    return [CONJ_MOVE[letter] for letter in reversed(conjugator)]


def check_route(core, u, v, sign, moves, challenge_id):
    state = (u, v)
    expected = expected_endpoint(u, v, sign)
    lengths = [len(u) + len(v)]
    for move in moves:
        state = core.apply_move(state, move)
        assert state[1] == v, 'donor changed during replay'
        lengths.append(sum(map(len, state)))
    assert state == expected, (u, v, sign, moves, state, expected)
    challenge = {
        'challenge_id': challenge_id,
        'move_spec_version': core.MOVE_SPEC_VERSION,
        'initial_relators': [list(u), list(v)],
        'target_relators': [list(w) for w in expected],
    }
    COUNTS['official_verify_calls'] += 1
    result = core.verify(challenge, list(moves), core.MOVE_SPEC_VERSION, {
        'max_path_length': 100000,
        'max_total_relator_length': 10000,
        'max_work': 5000000,
    })
    assert result['ok'], result
    assert result['work'] == sum(lengths)
    assert result['peak_total_relator_length'] == max(lengths)
    COUNTS['route_checks'] += 1
    return {
        'moves': list(moves),
        'endpoint': [list(w) for w in state],
        'total_relator_lengths_including_initial': lengths,
        'official_verdict': result,
    }


def check_pair(core, u, v, prefix):
    for sign in (1, -1):
        suffix = 'positive' if sign == 1 else 'negative'
        literal = check_route(core, u, v, sign, literal_moves(v, sign),
                              f'{prefix}-{suffix}-literal')
        compiled = check_route(core, u, v, sign, COMPILED[sign],
                               f'{prefix}-{suffix}-compiled')
        assert literal['endpoint'] == compiled['endpoint']


def random_word(rng, length):
    word = []
    for _ in range(length):
        choices = [letter for letter in LETTERS
                   if not word or letter != -word[-1]]
        word.append(rng.choice(choices))
    return tuple(word)


def resource_guard_counterexample(core):
    """Endpoint equality does not imply feasibility under official caps."""
    u, v = (1,) * 9992, (1,) * 5
    challenge = {
        'challenge_id': 'synthetic-donor-resource-guard',
        'move_spec_version': core.MOVE_SPEC_VERSION,
        'initial_relators': [list(u), list(v)],
        'target_relators': [list(u), list(v)],
    }
    limits = {'max_path_length': 100000,
              'max_total_relator_length': 10000, 'max_work': 5000000}
    output = {}
    for name, moves in (('literal', [6] * 5), ('compiled', COMPILED[1])):
        state = (u, v)
        lengths = [sum(map(len, state))]
        for move in moves:
            state = core.apply_move(state, move)
            lengths.append(sum(map(len, state)))
        assert state == (u, v)
        COUNTS['official_verify_calls'] += 1
        verdict = core.verify(challenge, moves, core.MOVE_SPEC_VERSION, limits)
        output[name] = {
            'moves': moves,
            'endpoint_equals_initial': True,
            'total_relator_lengths_including_initial': lengths,
            'uncapped_replay_work': sum(lengths),
            'uncapped_replay_peak': max(lengths),
            'official_capped_verdict': verdict,
        }
    assert output['literal']['official_capped_verdict']['ok']
    assert output['literal']['official_capped_verdict']['work'] == 59982
    rejected = output['compiled']['official_capped_verdict']
    assert not rejected['ok'] and rejected['code'] == 'E_LENGTH_LIMIT'
    assert rejected['move_index'] == 1
    assert rejected['total_relator_length'] == 10002
    return {
        'description': '(x^9992, x^5); same endpoint, cheaper rule violates length cap',
        'initial_relator_lengths': [9992, 5], 'limits': limits,
        'literal': output['literal'], 'compiled': output['compiled'],
        'consequence': 'Always replay a proposed replacement with full official resource guards; an endpoint identity is not a universally legal optimization.',
    }


def frozen_dev_scan():
    """One fixed scan of maximal same-target conjugation blocks on four dev paths."""
    source = ROOT / 'competition/examples/stable_training_424.json'
    if not source.exists():
        return {'performed': False, 'reason': 'local source absent; no fetch'}
    stable_path = VERIFIER / 'stable_core.py'
    expected_blob = 'af48e32f766dbcd0cad0becd8c1c4edb4358efac'
    assert git_blob(stable_path) == expected_blob
    name = '_pinned_donor_rule_verifier.stable_core'
    spec = importlib.util.spec_from_file_location(name, stable_path)
    stable = importlib.util.module_from_spec(spec)
    sys.modules[name] = stable
    spec.loader.exec_module(stable)
    COUNTS['official_stable_apply_move_calls'] = 0
    COUNTS['official_stable_verify_calls'] = 0
    original_apply = stable.apply_move

    def counted_apply(state, move):
        COUNTS['official_stable_apply_move_calls'] += 1
        return original_apply(state, move)

    stable.apply_move = counted_apply
    dev_ids = ('ms-train-0379', 'ms-train-0160',
               'ms-train-0013', 'ms-train-0015')
    # Inspect moves only after exact membership in the frozen development set.
    records = {record['training_id']: record
               for record in json.loads(source.read_text())['instances']
               if record['training_id'] in dev_ids}
    assert set(records) == set(dev_ids)
    limits = {'max_path_length': 100000,
              'max_total_relator_length': 10000, 'max_work': 5000000}
    results = []
    hits = []
    for dev_id in dev_ids:
        record = records[dev_id]
        moves = record['moves']
        challenge = {
            'challenge_id': dev_id,
            'move_spec_version': record['move_spec_version'],
            'initial_relators': record['initial_relators'],
            'target_relators': record['target_relators'],
        }
        COUNTS['official_stable_verify_calls'] += 1
        original_verdict = stable.verify(
            challenge, moves, stable.MOVE_SPEC_VERSION, limits)
        assert original_verdict['ok'], original_verdict
        states = [tuple(tuple(w) for w in record['initial_relators'])]
        for move in moves:
            state, error = stable.apply_move(states[-1], move)
            assert error is None
            states.append(state)
        assert states[-1] == tuple(tuple(w) for w in record['target_relators'])
        blocks = []
        start = 0
        while start < len(moves):
            row = stable.MOVE_TABLE[moves[start]]
            if row['category'] != 'conjugation':
                start += 1
                continue
            target = row['relator']
            end = start
            conjugator = ()
            while end < len(moves):
                next_row = stable.MOVE_TABLE[moves[end]]
                if (next_row['category'] != 'conjugation'
                        or next_row['relator'] != target):
                    break
                conjugator = reduced((next_row['conjugator'],) + conjugator)
                end += 1
            if end - start > 4:
                blocks.append({'start': start, 'end_exclusive': end,
                               'length': end - start, 'target_relator': target})
                for donor, donor_word in enumerate(states[start]):
                    if donor == target:
                        continue
                    for sign in (1, -1):
                        required = donor_word if sign == 1 else inverse(donor_word)
                        if conjugator != required:
                            continue
                        invert_id = next(r['id'] for r in stable.MOVE_TABLE
                                         if r['category'] == 'inversion'
                                         and r['relator'] == target)
                        multiply_id = next(r['id'] for r in stable.MOVE_TABLE
                                           if r['category'] == 'multiplication'
                                           and r['relator'] == target
                                           and r['other'] == donor
                                           and r['invert_other'] == (sign == 1))
                        replacement = [invert_id, multiply_id,
                                       invert_id, multiply_id]
                        candidate = moves[:start] + replacement + moves[end:]
                        COUNTS['official_stable_verify_calls'] += 1
                        verdict = stable.verify(
                            challenge, candidate, stable.MOVE_SPEC_VERSION, limits)
                        hits.append({'dev_id': dev_id, 'start': start,
                                     'end_exclusive': end, 'donor_relator': donor,
                                     'sign': sign, 'replacement': replacement,
                                     'original_length': len(moves),
                                     'candidate_length': len(candidate),
                                     'whole_path_official_verdict': verdict})
            start = end
        results.append({'dev_id': dev_id, 'path_length': len(moves),
                        'original_whole_path_verdict': original_verdict,
                        'eligible_maximal_blocks': blocks})
    return {
        'performed': True, 'source': str(source),
        'stable_core_git_blob': expected_blob,
        'frozen_dev_ids': list(dev_ids), 'limits': limits,
        'eval_path_inspection': False,
        'scan_definition': 'Maximal consecutive conjugations of one target, length > 4; accumulated reduced conjugator must equal an unchanged other relator or its inverse.',
        'tuning_after_results': False,
        'records': results, 'hits': hits, 'hit_count': len(hits),
        'whole_path_validated_strict_shorter_count': sum(
            hit['whole_path_official_verdict']['ok'] for hit in hits),
    }


def main():
    global ROOT, VERIFIER
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--official-root', type=Path, default=ROOT,
                        help='Root of the cached pinned official checkout')
    parser.add_argument('--skip-dev-scan', action='store_true',
                        help='Run only the rule checks and resource guard fixture')
    args = parser.parse_args()
    ROOT = args.official_root.resolve()
    VERIFIER = ROOT / 'competition/tools/verifier'
    core = load_official_core()
    words = [()]
    for length in (1, 2):
        words.extend(word for word in itertools.product(LETTERS, repeat=length)
                     if reduced(word) == word)
    assert len(words) == 17
    for index, (u, v) in enumerate(itertools.product(words, repeat=2)):
        check_pair(core, u, v, f'synthetic-exhaustive-{index:03d}')

    seed = 20261117
    rng = random.Random(seed)
    random_cases = []
    for index in range(100):
        u = random_word(rng, rng.randint(3, 40))
        v = random_word(rng, rng.randint(3, 40))
        check_pair(core, u, v, f'synthetic-random-{index:03d}')
        random_cases.append({'u': list(u), 'v': list(v)})

    # Deliberately synthetic fixture: (x, x^11 y), not an official challenge.
    u, v = (1,), (1,) * 11 + (2,)
    literal = check_route(core, u, v, 1, [8] + [6] * 11,
                          'synthetic-donor-fixture-literal')
    compiled = check_route(core, u, v, 1, [0, 3, 0, 3],
                           'synthetic-donor-fixture-compiled')
    assert literal['endpoint'] == compiled['endpoint']
    assert literal['official_verdict']['work'] == 325
    assert compiled['official_verdict']['work'] == 113
    assert literal['official_verdict']['peak_total_relator_length'] == 37
    assert compiled['official_verdict']['peak_total_relator_length'] == 37
    guard_counterexample = resource_guard_counterexample(core)
    dev_scan = ({'performed': False, 'reason': 'explicit --skip-dev-scan'}
                if args.skip_dev_scan else frozen_dev_scan())

    report = {
        'status': 'passed',
        'scope': 'local pinned AC rule checks plus a separate frozen-four-development certificate scan',
        'native_host_calls': 0,
        'official_challenge_improvement_claimed': False,
        'mathematical_novelty_claimed': False,
        'provenance': {
            'expected_official_commit': PINNED_COMMIT,
            'core_path': str(VERIFIER / 'core.py'),
            'expected_core_git_blob': EXPECTED_CORE_BLOB,
            'actual_core_git_blob': git_blob(VERIFIER / 'core.py'),
            'core_git_blob_checked_before_execution': True,
            'canon_git_blob': git_blob(VERIFIER / 'canon.py'),
            'move_spec_version': core.MOVE_SPEC_VERSION,
        },
        'lemma': {
            'state': '(u,v), freely reduced words in F(x,y)',
            'positive_endpoint': '(free_reduce(v u v^-1), v)',
            'positive_moves': COMPILED[1],
            'negative_endpoint': '(free_reduce(v^-1 u v), v)',
            'negative_moves': COMPILED[-1],
            'proof_positive': [
                'u -> u^-1', '-> u^-1 v^-1', '-> v u', '-> v u v^-1'],
            'proof_negative': [
                'u -> u^-1', '-> u^-1 v', '-> v^-1 u', '-> v^-1 u v'],
            'donor_unchanged_at_every_step': True,
            'letterwise_route': 'map official conjugation IDs over reverse(v^sign)',
            'cost_caution': 'Compiled length 4 is shorter only when literal donor length > 4; work need not improve for every pair.',
        },
        'exhaustive': {
            'max_reduced_word_length': 2, 'words': 17, 'ordered_pairs': 289,
            'signs_per_pair': 2, 'semantic_cases': 578, 'route_checks': 1156,
            'empty_words_included': True, 'passed': True,
        },
        'random': {
            'seed': seed, 'ordered_pairs': 100, 'length_range_inclusive': [3, 40],
            'signs_per_pair': 2, 'semantic_cases': 200, 'route_checks': 400,
            'passed': True, 'cases': random_cases,
        },
        'fixture': {
            'synthetic': True, 'not_official_challenge': True,
            'initial_relators': [list(u), list(v)],
            'description': '(x, x^11 y), positive donor conjugation',
            'literal': literal, 'compiled': compiled,
            'move_reduction': 8, 'work_reduction': 212,
        },
        'counts': COUNTS,
        'resource_guard_counterexample': guard_counterexample,
        'frozen_dev_scan': dev_scan,
        'bridge_hypothesis': {
            'description': 'A certified proof may contain a letterwise conjugation whose conjugator equals the unchanged donor relator; reconstructing that segment in relator coordinates can lower expanded primitive cost.',
            'status': 'classical identity verified; frozen development scan reported separately',
            'real_path_scan_performed': dev_scan['performed'],
            'next_evidence_needed': 'Find a matching segment in a certified real path, replace it, and verify the entire official certificate and actual metric.',
        },
    }
    output = Path(__file__).with_name('checks.json')
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({
        'status': report['status'], 'output': str(output), 'counts': COUNTS,
        'fixture': {
            name: {key: route['official_verdict'][key]
                   for key in ('length', 'work', 'peak_total_relator_length')}
            for name, route in (('literal', literal), ('compiled', compiled))
        },
    }, indent=2))


if __name__ == '__main__':
    main()
