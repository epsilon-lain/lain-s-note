"""Fresh small check by reviewer; does not import any uploaded research source."""
from itertools import combinations, permutations
from collections import deque
import json


def subsets(items):
    items = tuple(items)
    for bits in range(1 << len(items)):
        yield tuple(items[i] for i in range(len(items)) if bits >> i & 1)


def matching(edges):
    best = 0
    for chosen in subsets(edges):
        flat = [v for edge in chosen for v in edge]
        if len(set(flat)) == len(flat):
            best = max(best, len(chosen))
    return best


def cover_leaves(n, edges, target):
    leaves = [1 << v for v in range(n)]
    leaves += [sum(1 << v for v in e) for e in edges]
    leaves = [leaf for leaf in leaves if leaf & target == leaf]
    dist = {0: 0}
    queue = deque([0])
    while queue:
        old = queue.popleft()
        for leaf in leaves:
            new = old | leaf
            if new not in dist:
                dist[new] = dist[old] + 1
                queue.append(new)
    return dist[target]


def star_triangle(n, edges):
    adj = {v: set() for v in range(n)}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    unseen = set(range(n))
    while unseen:
        todo = [unseen.pop()]
        vertices = set(todo)
        while todo:
            for v in adj[todo.pop()]:
                if v not in vertices:
                    vertices.add(v)
                    unseen.remove(v)
                    todo.append(v)
        degs = sorted(len(adj[v]) for v in vertices)
        k = len(vertices)
        if k == 1:
            continue
        if degs == [1] * (k - 1) + [k - 1]:
            continue
        if k == 3 and degs == [2, 2, 2]:
            continue
        return False
    return True


rows = []
for n in range(1, 5):
    all_edges = tuple(combinations(range(n), 2))
    formulas = submods = optimizations = 0
    for h in subsets(all_edges):
        for target in range(1, 1 << n):
            inside = tuple(e for e in h if all(target >> v & 1 for v in e))
            actual = 2 * cover_leaves(n, h, target) - 1
            predicted = 2 * (target.bit_count() - matching(inside)) - 1
            assert actual == predicted
            formulas += 1
        vals = {frozenset(d): matching(d) for d in subsets(h)}
        is_submod = True
        for d, val in vals.items():
            missing = set(h) - d
            for e, f in combinations(missing, 2):
                if vals[d | {e, f}] - vals[d | {f}] > vals[d | {e}] - val:
                    is_submod = False
        assert is_submod == star_triangle(n, h)
        submods += is_submod
        for q in range(1, 5):
            for delta in range(1, 7):
                brute = max(2 * q * m - delta * len(d) for d, m in vals.items())
                assert brute == max(0, 2 * q - delta) * matching(h)
                optimizations += 1
    rows.append(dict(n=n, graphs=1 << len(all_edges), targets=formulas,
                     submod_graphs=submods, dictionary_opt_checks=optimizations))

p4 = ((0, 1), (1, 2), (2, 3))
p4_gains = []
for order in permutations(p4):
    chosen = ()
    for e in order:
        if 4 * (matching(chosen + (e,)) - matching(chosen)) - 3 > 0:
            chosen += (e,)
    p4_gains.append(4 * matching(chosen) - 3 * len(chosen))
assert sorted(p4_gains) == [1, 1, 2, 2, 2, 2]

hetero = ((0, 1), (1, 2))
hetero_costs = []
for d in subsets(hetero):
    cost = 3 * len(d)
    for target in (0b011, 0b110):
        cost += 2 * (2 * cover_leaves(3, d, target) - 1)
    hetero_costs.append(cost)
assert hetero_costs == [12, 11, 11, 10]
print(json.dumps(dict(rows=rows, p4_gains=p4_gains, heterogeneous_costs=hetero_costs,
                     all_assertions_passed=True), ensure_ascii=False, indent=2))
