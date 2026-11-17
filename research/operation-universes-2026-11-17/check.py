#!/usr/bin/env python3
"""Exact small arithmetic-circuit checks; no network or third-party packages.

Free inputs: formal x, 0, 1. Binary gates cost one; fanout is free.
Finite searches supplement, and do not replace, the proofs in README.md.
"""
import json
from itertools import combinations_with_replacement
from pathlib import Path


def norm(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)


ZERO, ONE, X = (0,), (1,), (0, 1)
INITIAL = tuple(sorted((ZERO, ONE, X)))


def add(a, b):
    return norm([(a[i] if i < len(a) else 0) +
                 (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, av in enumerate(a):
        for j, bv in enumerate(b):
            out[i + j] += av * bv
    return norm(out)


OPS = {'+': add, '*': mul}


def run(gates):
    """A witness uses initial nodes [0, 1, x], then prior gate indices."""
    values = [ZERO, ONE, X]
    for op, i, j in gates:
        assert op in OPS and 0 <= i < len(values) and 0 <= j < len(values)
        values.append(OPS[op](values[i], values[j]))
    return values[-1]


def exact_search(target, allowed, max_gates):
    """BFS all useful value-set states within sound coefficient/degree bounds.

    In N[x], along any path to a nonzero output, degree and maximum
    coefficient cannot decrease under + or multiplication by a nonzero
    polynomial. Multiplication by zero erases that dependency. Thus a
    minimal nonzero-output circuit has no live intermediate exceeding
    either target bound. Duplicate values can be removed using free fanout.
    """
    target = norm(target)
    assert target != ZERO
    degree_bound, coefficient_bound = len(target) - 1, max(target)
    states = {INITIAL}
    layer_counts = [1]
    for depth in range(max_gates + 1):
        if any(target in state for state in states):
            return {'first_found': depth, 'excluded_through': depth - 1,
                    'layer_state_counts': layer_counts}
        if depth == max_gates:
            break
        next_states = set()
        for state in states:
            for a, b in combinations_with_replacement(state, 2):
                for op in allowed:
                    p = OPS[op](a, b)
                    if (p in state or len(p) - 1 > degree_bound or
                            max(p) > coefficient_bound):
                        continue
                    next_states.add(tuple(sorted((*state, p))))
        states = next_states
        layer_counts.append(len(states))
    return {'first_found': None, 'excluded_through': max_gates,
            'layer_state_counts': layer_counts}


def mixed_linear_witness(k):
    gates = [('+', 1, 1)]
    current = 3
    for _ in range(k):
        gates.append(('*', current, current))
        current += 1
    gates.append(('*', current, 2))
    return gates


def addition_linear_witness(k):
    gates, current = [], 2
    for _ in range(2 ** k):
        gates.append(('+', current, current))
        current = len(gates) + 2
    return gates


def shifted_power_witness(k):
    gates, current = [('+', 2, 1)], 3
    for _ in range(k):
        gates.append(('*', current, current))
        current += 1
    return gates


def main():
    searches = []
    cases = [
        ('mixed_x_squared_plus_x', (0, 1, 1), ('+', '*'), 2, 2),
        ('mixed_shifted_square', (1, 2, 1), ('+', '*'), 2, 2),
        ('mixed_shifted_fourth', (1, 4, 6, 4, 1), ('+', '*'), 3, 3),
        ('mixed_16x', (0, 16), ('+', '*'), 4, 4),
        ('mixed_256x_exclude_up_to_4', (0, 256), ('+', '*'), 4, None),
        ('plus_16x', (0, 16), ('+',), 4, 4),
    ]
    for name, target, allowed, limit, expected in cases:
        result = exact_search(target, allowed, limit)
        assert result['first_found'] == expected, (name, result)
        searches.append({'name': name, 'target_coefficients': target,
                         'allowed': allowed, 'max_gates': limit, **result})

    witnesses = []
    for k in range(2, 9):
        alpha = 2 ** (2 ** k)
        mixed = mixed_linear_witness(k)
        plus = addition_linear_witness(k)
        assert run(mixed) == run(plus) == (0, alpha)
        assert len(mixed) == k + 2 and len(plus) == 2 ** k
        witnesses.append({'k': k, 'alpha_decimal': str(alpha),
                          'mixed_gate_count': len(mixed),
                          'plus_gate_count': len(plus),
                          'mixed_witness': mixed,
                          'polynomial_identity_verified': True})
    shifted = []
    for k in range(1, 7):
        gates = shifted_power_witness(k)
        actual = run(gates)
        from math import comb
        expected = tuple(comb(2 ** k, j) for j in range(2 ** k + 1))
        assert actual == expected and len(gates) == k + 1
        shifted.append({'k': k, 'degree': 2 ** k, 'gates': len(gates),
                        'nonzero_expanded_coefficients': len(actual),
                        'identity_verified': True})
    assert run([('+', 2, 2), ('+', 3, 3)]) == (0, 4)
    result = {
        'status': 'passed',
        'model': 'N[x]; free 0,1,x; binary +,*; unit gate count; free reuse',
        'scope': 'finite exact BFS plus symbolic witness checks; general lower bounds are mathematical proofs',
        'searches': searches,
        'linear_witness_checks': witnesses,
        'shifted_power_witness_checks': shifted,
        'k_1_linear_exception_verified': True,
        'not_claimed': ['bit-complexity speedup', 'novel theorem priority',
                        'AC competition improvement', 'native Yohaku execution'],
    }
    path = Path(__file__).with_name('checks.json')
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'status': result['status'],
                      'exact_searches': searches,
                      'linear_witnesses': len(witnesses),
                      'shifted_power_witnesses': len(shifted)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
