"""New external small check; does not execute or import the native source.

Only Python standard library. All labelled simple graphs on 0..4 vertices.
Run: python check_pair_inverse.py --output external_check_results.json
"""
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import argparse
import hashlib
import json


def powerset(items):
    items = tuple(items)
    return tuple(frozenset(c) for k in range(len(items) + 1)
                 for c in combinations(items, k))


def observe(x, edges):
    return frozenset(e for e in edges if e <= x)


def direct_query_depth(candidates, target, vertices):
    """Exact deterministic minimax DP; target labels come from enumeration.

    A state is the entire remaining candidate set. All original coordinates
    are eligible, including those outside the target coordinate set. Constant
    coordinates cannot help, so only strictly splitting queries are explored.
    """
    @lru_cache(None)
    def solve(current):
        if len({target[x] for x in current}) <= 1:
            return 0
        costs = []
        for v in vertices:
            yes = frozenset(x for x in current if v in x)
            no = current - yes
            if yes and no:
                costs.append(1 + max(solve(yes), solve(no)))
        assert costs
        return min(costs)
    return solve(frozenset(candidates))


def run():
    summaries = []
    for n in range(5):
        vertices = frozenset(range(n))
        possible_edges = tuple(frozenset(e) for e in combinations(range(n), 2))
        subsets = powerset(vertices)
        row = dict(n=n, graphs=0, objects=0, output_subsets=0,
                   legal_observations=0, target_fiber_cases=0, query_dp_cases=0)
        for edges in powerset(possible_edges):
            row['graphs'] += 1
            fibers = defaultdict(set)
            for x in subsets:
                fibers[observe(x, edges)].add(x)
                row['objects'] += 1
            global_patterns = {s: 0 for s in subsets}
            global_weights = {s: 0 for s in subsets}
            for f in powerset(edges):
                row['output_subsets'] += 1
                p = frozenset().union(*f) if f else frozenset()
                z = frozenset(v for e in edges for v in e
                              if e & p and v not in p)
                r = vertices - p - z
                legal = observe(p, edges) == f
                actual = fibers.get(f, set())
                assert legal == bool(actual)
                if not legal:
                    continue
                row['legal_observations'] += 1
                independent_r = {i for i in powerset(r) if not observe(i, edges)}
                assert actual == {p | i for i in independent_r}
                for s in subsets:
                    row['target_fiber_cases'] += 1
                    a = r & s
                    independent_a = {j for j in powerset(a) if not observe(j, edges)}
                    patterns = {x & s for x in actual}
                    weights = {len(x & s) for x in actual}
                    alpha = max(map(len, independent_a))
                    assert patterns == {(p & s) | j for j in independent_a}
                    assert weights == set(range(len(p & s), len(p & s) + alpha + 1))
                    assert (len(patterns) == 1) == (not a)
                    assert (len(weights) == 1) == (not a)
                    global_patterns[s] = max(global_patterns[s], len(patterns))
                    global_weights[s] = max(global_weights[s], len(weights))
                    for label in ({x: x & s for x in actual},
                                  {x: len(x & s) for x in actual}):
                        depth = direct_query_depth(actual, label, vertices)
                        assert depth == len(a), (n, edges, f, s, depth)
                        row['query_dp_cases'] += 1
            for s in subsets:
                independent_s = {i for i in powerset(s) if not observe(i, edges)}
                assert global_patterns[s] == len(independent_s)
                assert global_weights[s] == max(map(len, independent_s)) + 1
                # Exact ceiling(log2(m)), also valid for m=1.
                assert (global_patterns[s] - 1).bit_length() == (len(independent_s) - 1).bit_length()
                assert (global_weights[s] - 1).bit_length() == max(map(len, independent_s)).bit_length()
        summaries.append(row)
    return {
        'format': 'external-pair-inverse-small-check/1',
        'status': 'passed',
        'source_commit': '43c3c2e9b026cb098f9c33d338a3bd2be34c24b5',
        'scope': 'all labelled finite simple graphs n=0..4; all objects, edge observations and target subsets; two deterministic coordinate-query DPs per legal F/S',
        'provenance': 'new external review implementation; not native execution, not a rerun of native n=5 artifact',
        'algorithm': 'frozenset object enumeration and complete candidate-set minimax decision trees',
        'summaries': summaries,
        'totals': {k: sum(row[k] for row in summaries) for k in summaries[0] if k != 'n'},
        'check_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'limitations': ['finite evidence does not replace arbitrary-n proofs',
                        'no independent n=5 rerun',
                        'no native output artifact hash or full artifact reproduction',
                        'no algorithmic efficiency or benchmark claim'],
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('external_check_results.json'))
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'totals': result['totals']}, ensure_ascii=False))
