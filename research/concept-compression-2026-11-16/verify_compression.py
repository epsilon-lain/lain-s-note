"""Exact finite concept-expression examples; Python 3.10+, standard library only.

Run normally to rebuild exact_results.json, or use --check to compare with it.
No model, network, floating-point arithmetic, or contest solution is used here.
"""
from collections import deque
from itertools import product
from pathlib import Path
import argparse
import json


ROOT = Path(__file__).resolve().parent


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def tree_dp(target, atoms):
    """All nonnegative vectors <= target; atoms map to (cost, display name).

    Costs count every atom occurrence and every binary addition.  Since all
    atoms are nonzero, nonnegative, each proper split has smaller total mass.
    No overshooting subexpression can occur in a valid target expression.
    """
    states = sorted((v for v in product(*(range(x + 1) for x in target))
                     if any(v)), key=lambda v: (sum(v), v))
    best, split_count = {}, 0
    for v in states:
        choices = []
        if v in atoms:
            cost, name = atoms[v]
            choices.append((cost, name))
        for a in product(*(range(x + 1) for x in v)):
            b = tuple(x - y for x, y in zip(v, a))
            if not any(a) or not any(b) or a > b:
                continue
            split_count += 1
            if a in best and b in best:
                choices.append((best[a][0] + best[b][0] + 1,
                                f"({best[a][1]}+{best[b][1]})"))
        if choices:
            best[v] = min(choices)
    return {
        "state_count": len(states),
        "unordered_splits_checked": split_count,
        "target_cost": best[target][0],
        "target_expression": best[target][1],
        "all_state_minima": [
            {"vector": list(v), "cost": best[v][0],
             "expression": best[v][1]} for v in states
        ],
    }


def dag_bfs(target, bases):
    """Exact fewest addition gates when input nodes are shared and free to use.

    A state is the set of already computed vectors. Each step adds the sum of
    two available vectors (the same vector may be used twice), within target.
    Already available results are not added twice. Every acyclic straight-line
    addition circuit has a topological instruction order, so BFS is exhaustive.
    """
    initial = frozenset(bases)
    queue = deque([initial])
    parent = {initial: None}
    transition = {}
    depth = {initial: 0}
    while queue:
        state = queue.popleft()
        if target in state:
            path, current = [], state
            while parent[current] is not None:
                a, b, v = transition[current]
                path.append({"left": list(a), "right": list(b),
                             "result": list(v)})
                current = parent[current]
            path.reverse()
            return {"addition_gates": len(path),
                    "input_nodes": len(initial),
                    "total_nodes": len(initial) + len(path),
                    "instructions": path,
                    "first_target_bfs_depth": depth[state]}
        values = sorted(state)
        for i, a in enumerate(values):
            for b in values[i:]:
                v = add(a, b)
                if any(x > limit for x, limit in zip(v, target)) or v in state:
                    continue
                following = state | {v}
                if following in parent:
                    continue
                parent[following] = state
                transition[following] = (a, b, v)
                depth[following] = depth[state] + 1
                queue.append(following)
    raise AssertionError("target unreachable")


def graph_distance(start, goal, edges):
    queue = deque([(start, 0)])
    seen = {start}
    while queue:
        vertex, distance = queue.popleft()
        if vertex == goal:
            return distance
        for a, b in edges:
            if a == vertex and b not in seen:
                seen.add(b)
                queue.append((b, distance + 1))
    return None


def compute():
    target = (1, 1, 1, 1, 1)
    bases = {tuple(int(i == j) for i in range(5)): (1, f"e{j + 1}")
             for j in range(5)}
    baseline = tree_dp(target, bases)
    b, c = (1, 1, 0, 0, 0), (0, 0, 1, 1, 1)
    preexisting = tree_dp(target, bases | {b: (1, "b"), c: (1, "c")})
    expanded = tree_dp(target, bases | {b: (3, "b"), c: (5, "c")})
    full_target_named = tree_dp(target, bases | {target: (1, "w")})
    assert [baseline["target_cost"], preexisting["target_cost"],
            expanded["target_cost"], full_target_named["target_cost"]] == [9, 3, 9, 1]
    assert baseline["state_count"] == 31
    assert baseline["unordered_splits_checked"] == 90
    assert all(row["cost"] == 2 * sum(row["vector"]) - 1
               for row in baseline["all_state_minima"])
    repeated = []
    for r in range(1, 7):
        repeated.append({"uses": r, "expand_each_time": 9 * r,
                         "define_blocks_once": 8 + 3 * r,
                         "define_target_once": 9 + r,
                         "blocks_beat_expansion": 8 + 3 * r < 9 * r})
    assert [row["uses"] for row in repeated if row["blocks_beat_expansion"]] == list(range(2, 7))
    two_bases = {(1, 0): (1, "e1"), (0, 1): (1, "e2")}
    repeated_tree = tree_dp((2, 2), two_bases)
    repeated_dag = dag_bfs((2, 2), two_bases)
    assert repeated_tree["target_cost"] == 7
    assert repeated_dag["total_nodes"] == 4
    edges = [(0, 1), (1, 0), (1, 2), (2, 1)]
    distance_one = graph_distance(1, 0, edges)
    distance_two = graph_distance(2, 0, edges)
    assert (distance_one, distance_two) == (1, 2)
    return {
        "format": "concept-compression-exact/1",
        "scope": "declared finite nonnegative-vector languages; elementary examples, not new universal theory",
        "cost_model": "tree atom occurrence 1 except expanded block references b=3/c=5; addition 1; DAG example counts each distinct input/gate node once",
        "target": list(target),
        "basis_only": baseline,
        "preexisting_blocks": preexisting,
        "expanded_block_reference_costs": expanded,
        "target_given_a_unit_cost_name": full_target_named,
        "definition_accounting": {
            "binding_headers_cost": 0,
            "block_definition_body_costs": [3, 5],
            "single_use_total": 11,
            "repeated_use_scheme_comparisons_not_global_optima": repeated,
        },
        "shared_subexpression_example": {
            "target": [2, 2], "tree": repeated_tree, "dag": repeated_dag,
            "bound": "two input nodes and at least two addition gates; construction attains four nodes",
        },
        "description_path_counterexample": {
            "objects": [0, 1, 2], "literal_cost_each": 1,
            "directed_edges": [list(e) for e in edges],
            "target": 0, "target_description_cost": 1,
            "distance_from_1": distance_one,
            "distance_from_2": distance_two,
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = compute()
    destination = ROOT / "exact_results.json"
    if args.check:
        assert json.loads(destination.read_text(encoding="utf-8")) == result
        action = "checked"
    else:
        destination.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                               encoding="utf-8")
        action = "wrote"
    print(json.dumps({"action": action, "basis_tree_cost": result["basis_only"]["target_cost"],
                      "preexisting_block_tree_cost": result["preexisting_blocks"]["target_cost"],
                      "expanded_block_tree_cost": result["expanded_block_reference_costs"]["target_cost"],
                      "finite_states_per_five_vector_language": 31,
                      "tree_nodes_for_2_2": result["shared_subexpression_example"]["tree"]["target_cost"],
                      "dag_nodes_for_2_2": result["shared_subexpression_example"]["dag"]["total_nodes"]},
                     ensure_ascii=False))
