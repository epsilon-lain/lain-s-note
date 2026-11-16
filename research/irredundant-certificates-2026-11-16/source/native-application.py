"""Native finite geodesic/resource/one-flat-macro audit, not official CLI.
All words and relator positions are exact; no quotient or length pruning.
Sources SEED/EXTERNAL_SOURCES are external_chatgpt_work, not native discoveries.
"""
import hashlib
import json
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
import time

SEED = json.loads(r'''{"official_commit":"a0fd6e6f52d82c93ccc06ab91d81e8fa3678256e","official_manifest_sha256":"dc8587d7dc6855edfb8aad087f27f4db59c558f66bdf838b664c1634a357ad37","challenge_id":"ac-08491","move_spec_version":"ac-r2-v1","instance_hash":"sha256:2938e3e2872ea496fe609883c1e0e58c0be982a39e7b0469e61c2fe3fb36b25e","initial_relators":[[-2,-1],[-2,-2,-1,-1,2,2,1]],"target_relators":[[1],[2]],"witness":[7,4,12,1,12,2,5,5],"official_witness_verification":{"ok":true,"length":8,"peak_total_relator_length":9,"work":49,"certificate_hash":"sha256:33b48c0c493aee87616d07731abb7a5b1cd25e392eb93f89f1fd5ee1e1e38526"},"forward_ball":{"depth":3,"frontier_sizes":[1,13,129,1186],"states":1329,"primitive_evaluations":2002,"official_crosscheck_evaluations":2002,"peak_total_length":41,"length_pruning":false,"states_sha256":"e1f6dd0863a66abdb4bae8809fa1e39c6e66ee274d95549822f75c0a906bedae"},"target_ball":{"depth":4,"frontier_sizes":[1,10,77,574,4453],"states":5115,"primitive_evaluations":9268,"peak_total_length":17,"length_pruning":false,"global_shortest_within_ball":true,"official_crosscheck_evaluations":9268,"states_sha256":"3737904cb75b8f52f58027c33bd8f5c22a8f5a3379ac844961cf6d6db3158bfc"},"intersection_count":0,"excluded_move_lengths":{"minimum":0,"maximum":7},"global_shortest_AC_length":8,"proof":"All moves have inverses. A path of length <=7 would meet the complete forward radius-3 ball and target radius-4 ball. Their intersection is empty; the verified witness has length 8.","scope":"One fixed AC presentation with the exact ordered target and frozen 14-move table; no optimality claim for Stable AC or other challenges.","separate_from_pilot_arm_budgets":true}''')
EXTERNAL_SOURCES = json.loads(r'''[{"local_path":"next_research/yohaku/research/concept-compression-ac-2026-11-16/shortest_certificate.json","sha256":"6946fda0afee248115f7555cde7bc77ed30cc5f71692f20142057a88a42d1f30","fixed_source_url":"https://github.com/epsilon-lain/yohaku/blob/9afa13603ab26d435d7b5729d057e9b7de9bf8bb/research/concept-compression-ac-2026-11-16/shortest_certificate.json","revision":"9afa13603ab26d435d7b5729d057e9b7de9bf8bb","remote_verification":"direct_git_blob_match"},{"local_path":"two_teams_sol_20261116/application/certify.py","sha256":"d7235090980d9db6429fb78ba0e50a34c50ecf37dcb0dedb8a026eb4487efb7f","fixed_source_url":"https://github.com/epsilon-lain/yohaku/blob/9afa13603ab26d435d7b5729d057e9b7de9bf8bb/research/theory-application-sol-2026-11-16/application/certify.py","revision":"9afa13603ab26d435d7b5729d057e9b7de9bf8bb","remote_verification":"direct_git_blob_match"},{"local_path":"next_research/math/research/concept-compression-2026-11-16/README.md","sha256":"874b660f4a1c9f6846a7f437e9a7f86e7ee11ff463a663feb2922223690a6d7d","fixed_source_url":"https://github.com/epsilon-lain/lain-s-note/blob/51f6c2a413b829568efed24d1167e049d2756220/research/concept-compression-2026-11-16/README.md","revision":"51f6c2a413b829568efed24d1167e049d2756220","remote_verification":"direct_git_blob_match"}]''')

START = tuple(map(tuple, SEED["initial_relators"]))
TARGET = tuple(map(tuple, SEED["target_relators"]))
OLD_WITNESS = tuple(SEED["witness"])
INVERSE = (0, 1, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12)
GENS = (1, -1, 2, -2)
ACTION_TABLE = [
    "r0 <- inverse(r0)", "r1 <- inverse(r1)",
    "r0 <- reduce(r0 r1)", "r0 <- reduce(r0 inverse(r1))",
    "r1 <- reduce(r1 r0)", "r1 <- reduce(r1 inverse(r0))",
] + [f"r{i} <- reduce(({g},) r{i} ({-g},))" for i in (0, 1) for g in GENS]

def inv(w):
    return tuple(-x for x in w[::-1])

def free_reduce(word):
    stack = []
    for x in word:
        if stack and stack[-1] == -x:
            stack.pop()
        else:
            stack.append(x)
    return tuple(stack)

def apply(state, move):
    """A: stack normal form, frozen seed table."""
    out = list(state)
    if move < 2:
        out[move] = inv(out[move])
    elif move < 6:
        i = (move - 2) // 2
        rhs = state[1-i] if move % 2 == 0 else inv(state[1-i])
        out[i] = free_reduce(state[i] + rhs)
    else:
        i = (move - 6) // 4
        g = GENS[(move - 6) % 4]
        out[i] = free_reduce((g,) + state[i] + (-g,))
    return tuple(out)

def erase_pairs(word):
    """B: repeatedly erase first adjacent inverse pair, not stack reduction."""
    a = list(word)
    while True:
        for i in range(len(a)-1):
            if a[i] + a[i+1] == 0:
                del a[i:i+2]
                break
        else:
            return tuple(a)

def reference(state, move):
    """Explicit 14 cases independent of A's index arithmetic."""
    a, b = state
    ia = tuple(-x for x in reversed(a))
    ib = tuple(-x for x in reversed(b))
    if move == 0: return ia, b
    if move == 1: return a, ib
    if move == 2: return erase_pairs(a+b), b
    if move == 3: return erase_pairs(a+ib), b
    if move == 4: return a, erase_pairs(b+a)
    if move == 5: return a, erase_pairs(b+ia)
    if move == 6: return erase_pairs((1,)+a+(-1,)), b
    if move == 7: return erase_pairs((-1,)+a+(1,)), b
    if move == 8: return erase_pairs((2,)+a+(-2,)), b
    if move == 9: return erase_pairs((-2,)+a+(2,)), b
    if move == 10: return a, erase_pairs((1,)+b+(-1,))
    if move == 11: return a, erase_pairs((-1,)+b+(1,))
    if move == 12: return a, erase_pairs((2,)+b+(-2,))
    if move == 13: return a, erase_pairs((-2,)+b+(2,))
    raise ValueError(move)

def size(s):
    return sum(map(len, s))

def digest(states):
    data = json.dumps(sorted(states), separators=(",", ":")).encode("ascii")
    return hashlib.sha256(data).hexdigest()

def ball(source, radius):
    distances = {source: 0}
    parents = defaultdict(list)
    layers = [[source]]
    evals = 0
    max_size = size(source)
    for d in range(radius):
        next_layer = []
        for s in layers[-1]:
            for m in range(14):
                t = apply(s, m)
                assert t == reference(s, m), ("implementation_mismatch", s, m)
                assert apply(t, INVERSE[m]) == s, ("inverse_mismatch", s, m)
                assert reference(t, INVERSE[m]) == s
                evals += 1
                max_size = max(max_size, size(t))
                if t not in distances:
                    distances[t] = d + 1
                    next_layer.append(t)
                if distances[t] == d + 1:
                    parents[t].append((s, m))
        layers.append(sorted(next_layer))
    metadata = {
        "radius": radius, "states": len(distances),
        "layer_sizes": list(map(len, layers)),
        "layer_sha256": [digest(layer) for layer in layers],
        "all_state_sha256": digest(distances),
        "primitive_evaluations": evals,
        "second_implementation_crosschecks": evals,
        "inverse_checks_per_implementation": evals,
        "maximum_generated_total_length": max_size,
        "length_resource_pruning": False, "state_quotient": False,
    }
    return distances, parents, layers, metadata

def chains(source, parents):
    @lru_cache(None)
    def visit(state):
        if state == source:
            return ((),)
        return tuple(p + (m,) for prev, m in parents[state] for p in visit(prev))
    return visit

def chain_counts(source, parents):
    @lru_cache(None)
    def visit(state):
        return 1 if state == source else sum(visit(s) for s, m in parents[state])
    return visit

def trace(path, source=START):
    states = [source]
    for m in path:
        states.append(apply(states[-1], m))
    return states

def summary(states):
    return {
        "source": states[0], "target": states[-1],
        "expanded_length": len(states)-1,
        "peak": max(map(size, states)),
        "tail_work": sum(map(size, states[1:])),
    }

def compose(a, b):
    assert a["target"] == b["source"]
    return {
        "source": a["source"], "target": b["target"],
        "expanded_length": a["expanded_length"]+b["expanded_length"],
        "peak": max(a["peak"], b["peak"]),
        "tail_work": a["tail_work"]+b["tail_work"],
    }

def parse(path, macro):
    """Exact shortest tokenization with this one fixed flat macro."""
    n = len(path)
    dp = [None] * (n+1)
    dp[n] = (0, ())
    for i in range(n-1, -1, -1):
        options = [(1 + dp[i+1][0], (path[i],) + dp[i+1][1])]
        j = i + len(macro)
        if j <= n and path[i:j] == macro:
            options.append((1 + dp[j][0], ("M",) + dp[j][1]))
        dp[i] = min(options, key=lambda x: (x[0], repr(x[1])))
    return dp[0]

def expand(tokens, macro):
    return tuple(m for token in tokens for m in (macro if token == "M" else (token,)))

def compress(paths):
    candidates = sorted({p[i:i+k] for p in paths for k in (2, 3, 4)
                         for i in range(len(p)-k+1)})
    per_path = [{"path": p, "best_paid_cost": len(p), "macro": None,
                 "tokens": p} for p in paths]
    corpus_best = {"total_cost": sum(map(len, paths)), "macro": None,
                   "definition_cost": 0, "reference_tokens": sum(map(len, paths)),
                   "tokens": list(paths)}
    free_best = sum(map(len, paths))
    candidate_costs = []
    for macro in candidates:
        parses = [parse(p, macro) for p in paths]
        for i, (cost, tokens) in enumerate(parses):
            assert expand(tokens, macro) == paths[i]
            paid = len(macro) + cost
            if paid < per_path[i]["best_paid_cost"]:
                per_path[i] = {"path": paths[i], "best_paid_cost": paid,
                               "macro": macro, "tokens": tokens}
        token_sum = sum(cost for cost, tokens in parses)
        free_best = min(free_best, token_sum)
        paid = len(macro) + token_sum
        candidate_costs.append({"macro": macro, "definition_cost": len(macro),
                                "reference_tokens": token_sum, "total_cost": paid})
        if paid < corpus_best["total_cost"]:
            corpus_best = {
                "total_cost": paid, "macro": macro, "definition_cost": len(macro),
                "reference_tokens": token_sum,
                "tokens": [tokens for cost, tokens in parses],
            }
    macro = corpus_best["macro"]
    contexts = []
    for p, tokens in zip(paths, corpus_best["tokens"]):
        assert expand(tokens, macro or ()) == p
        s = START
        accumulated = summary([s])
        for token in tokens:
            body = macro if token == "M" else (token,)
            states = trace(body, s)
            segment = summary(states)
            accumulated = compose(accumulated, segment)
            if token == "M":
                contexts.append(segment)
            s = states[-1]
        assert accumulated == summary(trace(p))
        assert s == TARGET
    return {
        "language": "zero or one flat macro; body length 2..4; primitive/ref=1; definition=body length; header=0; no nesting",
        "selection_scope": "post-hoc finite corpus, not held-out search performance",
        "candidate_count": len(candidates),
        "all_candidate_costs": candidate_costs,
        "single_path_paid_cost_distribution": sorted(Counter(x["best_paid_cost"] for x in per_path).items()),
        "single_path_optima": per_path,
        "corpus_baseline_cost": sum(map(len, paths)),
        "corpus_best_paid": corpus_best,
        "corpus_best_if_definition_preexisting": free_best,
        "chosen_macro_contexts": contexts,
        "description_tokens_are_not_raw_moves": True,
        "expanded_path_and_resource_preservation_checks": len(paths),
    }

def write(name, value):
    Path("/output", name).write_text(json.dumps(value, indent=2) + "\n")

def main():
    started = time.monotonic()
    f, fp, fl, fm = ball(START, 4)
    b, bp, bl, bm = ball(TARGET, 4)
    f3 = {s for s, d in f.items() if d <= 3}
    common_34 = f3.intersection(b)
    common_44 = sorted(f.keys() & b.keys())
    assert len(f3) == SEED["forward_ball"]["states"]
    assert digest(f3) == SEED["forward_ball"]["states_sha256"]
    assert digest(b) == SEED["target_ball"]["states_sha256"]
    assert not common_34
    assert common_44
    assert all(f[s] == b[s] == 4 for s in common_44)
    cf, cb = chain_counts(START, fp), chain_counts(TARGET, bp)
    path_count = sum(cf(s) * cb(s) for s in common_44)
    if path_count > 1000000:
        raise RuntimeError("Enumeration guard: no path resource/compression completeness claim")
    vf, vb = chains(START, fp), chains(TARGET, bp)
    paths = sorted({p + tuple(INVERSE[m] for m in reversed(q))
                    for s in common_44 for p in vf(s) for q in vb(s)})
    assert len(paths) == path_count
    assert OLD_WITNESS in paths
    records, edges = [], set()
    states_by_layer = [set() for _ in range(9)]
    distribution = Counter()
    for p in paths:
        states = trace(p)
        assert len(p) == 8 and states[-1] == TARGET
        lengths = list(map(size, states))
        peak, work = max(lengths), sum(lengths)
        distribution[(peak, work)] += 1
        records.append({"moves": p, "states": states, "lengths": lengths,
                        "peak": peak, "work": work})
        for i, s in enumerate(states):
            states_by_layer[i].add(s)
        for s, m, t in zip(states, p, states[1:]):
            edges.add((s, m, t))
        assert compose(summary(states[:5]), summary(states[4:])) == summary(states)
        assert sum(map(size, states[:5])) + sum(map(size, states[4:])) - size(states[4]) == work
    pairs = sorted(distribution)
    frontier = [r for r in pairs if not any(q != r and q[0] <= r[0] and q[1] <= r[1] for q in pairs)]
    old = next(r for r in records if tuple(r["moves"]) == OLD_WITNESS)
    assert old["peak"] == SEED["official_witness_verification"]["peak_total_relator_length"]
    assert old["work"] == SEED["official_witness_verification"]["work"]
    best = min(records, key=lambda r: (r["peak"], r["work"], r["moves"]))
    compression = compress(paths)
    dag_counts = {START: 1}
    for i in range(8):
        next_counts = defaultdict(int)
        for s, m, t in edges:
            if s in states_by_layer[i] and t in states_by_layer[i+1]:
                next_counts[t] += dag_counts.get(s, 0)
        dag_counts = dict(next_counts)
    assert dag_counts == {TARGET: path_count}
    comparison = {
        "old_witness": old, "lexicographic_peak_work_best": best,
        "all_geodesics_old_peak_work_equal": len(pairs) == 1,
        "shorter_or_stable_or_official_cli_claim": False,
    }
    results = {
        "version": "native-application-geodesic-resource-v1", "status": "complete",
        "challenge_id": "ac-08491", "start": START, "target": TARGET,
        "action_table": ACTION_TABLE, "inverse_ids": INVERSE,
        "resource_definitions": {"peak": "max total length including initial",
                                 "work": "sum total lengths including initial and final"},
        "external_seed_sources": EXTERNAL_SOURCES,
        "adopted_math": [{"id": "external/ac-08491-distance-8", "revision": 1},
                         {"id": "external/macro-expansion-interface", "revision": 1}],
        "official_validation": {"new_official_cli_run": False, "platform_upload": False,
                                "archived_official_witness_verdict": SEED["official_witness_verification"]},
        "forward_ball": fm, "target_ball": bm,
        "forward_radius3_sha256": digest(f3), "intersection_34_count": len(common_34),
        "intersection_44_count": len(common_44), "exact_distance_in_frozen_table": 8,
        "all_8_step_labelled_paths": path_count,
        "geodesic_dag": {"layer_sizes": list(map(len, states_by_layer)), "labelled_edges": len(edges)},
        "resource_distribution": [{"peak": p, "work": w, "paths": distribution[(p,w)]} for p,w in pairs],
        "pareto_frontier": [{"peak": p, "work": w, "paths": distribution[(p,w)]} for p,w in frontier],
        "comparison": comparison, "compression": compression,
        "proof": {
            "lower_bound": "Every <=7 path has a vertex within 3 of start and 4 of target; those complete balls are disjoint. Inverse table checked on every generated edge.",
            "all_paths": "Every 8-step path is geodesic and has midpoint in both radius-4 balls. Every prefix and reverse suffix is a shortest parent chain, with all labelled parents retained.",
            "macro_optimum": "Any useful 2..4-ID flat macro occurs as a contiguous substring of an output path. Every such substring and every segmentation is considered by suffix DP; optional no-macro baseline included.",
        },
        "scope": "Frozen ordinary rank-2 AC action graph, exact ordered words. Complete only at length 8; no optimum over longer paths, Stable AC, general compression grammars, competition ranking, or cost claims.",
        "elapsed_seconds": time.monotonic()-started,
    }
    write("results.json", results)
    write("paths.json", {"records": records, "dag_edges": sorted(edges)})
    write("balls.json", {"forward": [{"state": s, "distance": f[s]} for s in sorted(f)],
                        "target": [{"state": s, "distance": b[s]} for s in sorted(b)]})
    Path("/output/submission-internal.txt").write_text("ac-08491: " + json.dumps(best["moves"]) + "\n")
    print(json.dumps({k: results[k] for k in [
        "status", "forward_ball", "target_ball", "intersection_34_count",
        "intersection_44_count", "all_8_step_labelled_paths", "geodesic_dag",
        "resource_distribution", "pareto_frontier", "elapsed_seconds"
    ]}), flush=True)
    print(json.dumps({"compression_candidates": compression["candidate_count"],
        "single_path_paid_cost_distribution": compression["single_path_paid_cost_distribution"],
        "corpus_baseline": compression["corpus_baseline_cost"],
        "corpus_best_paid_cost": compression["corpus_best_paid"]["total_cost"],
        "corpus_best_macro": compression["corpus_best_paid"]["macro"],
        "best_path": best["moves"], "best_peak": best["peak"], "best_work": best["work"],
        "new_official_cli_run": False}), flush=True)

if __name__ == "__main__":
    main()
