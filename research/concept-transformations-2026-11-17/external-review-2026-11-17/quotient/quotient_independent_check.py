#!/usr/bin/env python3
"""External bounded audit, stdlib only; never imports/executes author source.

Run: python quotient_independent_check.py [output.json]
Different reference methods: normalized arbitrary maps for partitions, Floyd-
Warshall distances, Cartesian complete trajectories for route endpoints.
"""
from collections import Counter
from itertools import product
import hashlib
import json
from pathlib import Path
import sys

COMMIT = "43c3c2e9b026cb098f9c33d338a3bd2be34c24b5"


def partitions(n):
    ans = set()
    for f in product(range(n), repeat=n):
        blocks = tuple(sorted((tuple(i for i in range(n) if f[i] == j)
                               for j in set(f)), key=lambda b: b[0]))
        p = tuple(next(j for j, b in enumerate(blocks) if i in b)
                  for i in range(n))
        ans.add(p)
    return sorted(ans)


def floyd(n, edges):
    d = [[0 if i == j else float("inf") for j in range(n)]
         for i in range(n)]
    for i, j in edges:
        d[i][j] = min(d[i][j], 1)
    for k in range(n):
        for i in range(n):
            for j in range(n):
                d[i][j] = min(d[i][j], d[i][k] + d[k][j])
    return d


def distance(d, s, target):
    return min((d[s][t] for t in target), default=float("inf"))


def complete_walks(n, edges):
    return [w for length in range(5)
            for w in product(range(n), repeat=length + 1)
            if all((u, v) in edges for u, v in zip(w, w[1:]))]


def forward(n, edges, p, s, w):
    layer = {s}
    for y in w[1:]:
        layer = {v for v in range(n) if p[v] == y
                 and any((u, v) in edges for u in layer)}
    return layer


def backward(n, edges, p, w, target):
    layer = {v for v in target if p[v] == w[-1]}
    for y in reversed(w[:-1]):
        layer = {u for u in range(n) if p[u] == y
                 and any((u, v) in edges for v in layer)}
    return layer


def run():
    total = Counter()
    rows = []
    for n in range(1, 4):
        c = Counter()
        slots = [(u, v) for u in range(n) for v in range(n) if u != v]
        for mask in range(1 << len(slots)):
            edges = {e for i, e in enumerate(slots) if mask >> i & 1}
            c["graphs"] += 1
            raw = floyd(n, edges)
            concrete_walks = complete_walks(n, edges)
            for p in partitions(n):
                c["graph_partitions"] += 1
                k = max(p) + 1
                qedges = {(p[u], p[v]) for u, v in edges}
                qdist = floyd(k, qedges)
                local = all({p[v] for u, v in edges if u == x}
                            == {b for a, b in qedges if a == p[x]}
                            for x in range(n))
                c["local_graph_partitions"] += local
                for s in range(n):
                    for bmask in range(1 << k):
                        b = {i for i in range(k) if bmask >> i & 1}
                        target = {i for i in range(n) if p[i] in b}
                        dr = distance(raw, s, target)
                        dq = distance(qdist, p[s], b)
                        assert dq <= dr
                        c["distance_checks"] += 1
                        c["strict_optimistic_distances"] += dq < dr
                        if local:
                            assert dq == dr
                            c["local_distance_equalities"] += 1
                # Exhaustive Cartesian trajectory enumeration is the reference.
                c["raw_tree_walks"] += len(concrete_walks)
                lookup = {}
                for w in concrete_walks:
                    key = (w[0], tuple(p[x] for x in w))
                    lookup.setdefault(key, set()).add(w[-1])
                all_lift = True
                for w in complete_walks(k, qedges):
                    for s in range(n):
                        if p[s] != w[0]:
                            continue
                        ends = lookup.get((s, w), set())
                        assert forward(n, edges, p, s, w) == ends
                        c["route_endpoint_checks"] += 1
                        c["spurious_routes"] += not ends
                        all_lift = all_lift and bool(ends)
                        if local:
                            assert ends
                            c["local_routes"] += 1
                        for tmask in range(1 << n):
                            target = {i for i in range(n) if tmask >> i & 1}
                            good = bool(ends & target)
                            assert (s in backward(n, edges, p, w, target)) == good
                            c["backward_target_checks"] += 1
                # Includes length 1, so this checks both directions of local iff.
                assert all_lift == local
                c["bounded_lifting_equivalence_checks"] += 1
        rows.append({"n": n, **dict(c)})
        total.update(c)
    # Explicit pedagogical examples, including n=4 without exhaustive n=4.
    examples = [
        ("three_state_single_edge", 3, {(1, 2)}, (0, 0, 1), 0, (0, 1), {2}),
        ("four_state_incompatible_middle", 4, {(0, 1), (2, 3)},
         (0, 1, 1, 2), 0, (0, 1, 2), {3}),
    ]
    example_results = []
    for name, n, edges, p, s, w, target in examples:
        assert not forward(n, edges, p, s, w)
        assert all((a, b) in {(p[u], p[v]) for u, v in edges}
                   for a, b in zip(w, w[1:]))
        example_results.append({"name": name, "edges": sorted(edges),
                                "representation": p, "source": s,
                                "quotient_route": w, "target": sorted(target),
                                "quotient_route_exists": True,
                                "concrete_lift_exists": False})
    # Local lifting alone does not suffice for an unsaturated target.
    e = {(0, 1), (1, 0)}
    assert distance(floyd(2, e), 0, {1}) == 1
    assert distance(floyd(1, {(0, 0)}), 0, {0}) == 0
    # n=4 coverage formulas are arithmetic only, never graph execution.
    stirling = {1: [1], 2: [1, 1], 3: [1, 3, 1], 4: [1, 7, 6, 1]}
    formula_rows = []
    for n, values in stirling.items():
        g = 2 ** (n * (n - 1))
        weighted = sum(v * 2 ** k for k, v in enumerate(values, 1))
        formula_rows.append({"n": n, "graphs": g,
                             "graph_partitions": g * sum(values),
                             "distance_checks": n * g * weighted})
    expected_routes = {"raw_tree_walks": 5797, "route_endpoint_checks": 7713,
                       "local_routes": 2461, "backward_target_checks": 61498,
                       "spurious_routes": 2909}
    assert all(total[key] == value for key, value in expected_routes.items())
    return {
        "stage": "external_round3_bounded_quotient_review",
        "reviewed_commit": COMMIT,
        "scope": "Independent stdlib implementation; n=1..3 exhaustive loop-free directed graphs; all set partitions, sources, saturated targets; walks length 0..4 and all actual targets. Two explicit n=3/4 examples only. n=4 full counts are formulas, not a rerun.",
        "methods": ["Floyd-Warshall", "Cartesian complete concrete trajectories",
                    "partitions deduplicated from all arbitrary maps"],
        "passed": True,
        "counts_by_n": rows, "counts_n_le_3": dict(total),
        "formula_only_n_le_4": formula_rows,
        "matches_author_route_aggregates": expected_routes,
        "examples": example_results,
        "limitations": ["No execution of uploaded/frozen author source",
                        "No full n=4 independent rerun",
                        "No external task/provenance authentication",
                        "No labels, resource history, application benchmark or official CLI",
                        "Bounded experiment supports but does not prove general theorem"],
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name("quotient_check.json")
    result = run()
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"passed": result["passed"], "counts": result["counts_n_le_3"]}, ensure_ascii=False))
