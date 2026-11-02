#!/usr/bin/env python3
"""Exact batch root-family search, without continued fractions as search input.

Target: all primitive p/q with 1 <= q <= 10**100 and
q*q*abs(pi-p/q) < 1/250. Dyadic denominator blocks are rejected by exact
integer counts in widened affine strips. Positive counts trigger bisection.
"""

from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path
import json

from build_cf_certificate import atan_bounds, pack, unpack, quality_bounds

ROOT = Path(__file__).resolve().parent


def floor_sum(n, m, a, b):
    """Return sum(floor((a*i+b)/m), i=0..n-1), including negative a,b."""
    assert n >= 0 and m > 0
    qa, a = divmod(a, m)
    qb, b = divmod(b, m)
    answer = qa * n * (n - 1) // 2 + qb * n
    steps = 0
    while True:
        steps += 1
        qa, a = divmod(a, m)
        qb, b = divmod(b, m)
        answer += qa * n * (n - 1) // 2 + qb * n
        y = a * n + b
        if y < m:
            return answer, steps
        n, b = divmod(y, m)
        m, a = a, m


def affine_floor_sum(A, B, slope, offset):
    """Exact sum floor(slope*q+offset), q=A..B."""
    denominator = lcm(slope.denominator, offset.denominator)
    a = slope.numerator * (denominator // slope.denominator)
    b = offset.numerator * (denominator // offset.denominator)
    return floor_sum(B - A + 1, denominator, a, b + A * a)


def strip_count(A, B, L, U, threshold):
    # A true hit obeys |q*pi-p| < threshold/q <= threshold/A.
    h = threshold / A
    lower_floor, steps1 = affine_floor_sum(A, B, L, -h)
    neg_upper_floor, steps2 = affine_floor_sum(A, B, -U, -h)
    count = -neg_upper_floor - lower_floor - (B - A + 1)
    assert count >= 0
    return count, lower_floor, neg_upper_floor, steps1 + steps2


def brute_strip_count(A, B, L, U, threshold):
    h = threshold / A
    answer = 0
    for q in range(A, B + 1):
        lower, upper = q * L - h, q * U + h
        first = lower.numerator // lower.denominator + 1
        last = -((-upper).numerator // (-upper).denominator) - 1
        answer += max(0, last - first + 1)
    return answer


def leaf_candidates(q, L, U, threshold):
    lower, upper = q * L - threshold / q, q * U + threshold / q
    first = lower.numerator // lower.denominator + 1
    last = -((-upper).numerator // (-upper).denominator) - 1
    rows = []
    for p in range(first, last + 1):
        d = gcd(p, q)
        try:
            quality_lower, quality_upper, side = quality_bounds(p, q, L, U)
        except AssertionError:
            value = F(p, q)
            quality_lower = F(0)
            quality_upper = q * q * max(value - L, U - value)
            side = "inside_pi_interval"
        if quality_upper < threshold:
            status = "retained_strictly_below_threshold"
        elif quality_lower >= threshold:
            status = "discarded_at_or_above_threshold"
        else:
            status = "unresolved_by_pi_interval"
        if d != 1:
            primitive_status = "nonprimitive_root_label"
        else:
            primitive_status = "primitive_candidate"
        rows.append({
            "p": str(p), "q": str(q), "gcd": str(d),
            "root_index_p_minus_3q": str(p - 3 * q),
            "reduced_p": str(p // d), "reduced_q": str(q // d),
            "primitive_status": primitive_status, "side": side,
            "raw_quality_lower": pack(quality_lower),
            "raw_quality_upper": pack(quality_upper),
            "reduced_quality_lower": pack(quality_lower / (d * d)),
            "reduced_quality_upper": pack(quality_upper / (d * d)),
            "status": status,
        })
    return rows


def run_search(N, L, U, threshold):
    # This assertion also guarantees that a leaf's uncertainty strip contains
    # at most one integer numerator. It is a guard, not a search assumption.
    assert N * (U - L) + 2 * threshold < 1
    nodes, roots, candidates = [], [], []
    total_steps = 0

    def visit(A, B, parent=None):
        nonlocal total_steps
        node_id = len(nodes)
        count, lo_floor, neg_hi_floor, steps = strip_count(A, B, L, U, threshold)
        total_steps += steps
        node = {"id": node_id, "parent": parent, "A": str(A), "B": str(B),
                "count": str(count), "lower_floor_sum": str(lo_floor),
                "negative_upper_floor_sum": str(neg_hi_floor),
                "floor_sum_euclidean_steps": steps}
        nodes.append(node)
        if count == 0:
            node["status"] = "zero_count_block_rejected"
        elif A == B:
            node["status"] = "leaf_scored"
            rows = leaf_candidates(A, L, U, threshold)
            assert len(rows) == count
            node["candidate_indices"] = list(range(len(candidates), len(candidates) + len(rows)))
            candidates.extend(rows)
        else:
            node["status"] = "bisected"
            middle = (A + B) // 2
            left_id = visit(A, middle, node_id)
            right_id = visit(middle + 1, B, node_id)
            node["children"] = [left_id, right_id]
        return node_id

    A = 1
    while A <= N:
        B = min(2 * A - 1, N)
        roots.append(visit(A, B))
        A *= 2

    retained = [r for r in candidates if r["status"] == "retained_strictly_below_threshold"]
    discarded = [r for r in candidates if r["status"] == "discarded_at_or_above_threshold"]
    unresolved = [r for r in candidates if r["status"] == "unresolved_by_pi_interval"]
    return {
        "denominator_cap": str(N), "threshold": pack(threshold),
        "pi_lower": pack(L), "pi_upper": pack(U),
        "search_uses_continued_fractions": False,
        "strict_strip": "q*L-threshold/A < p < q*U+threshold/A, A<=q<=B",
        "count_formula": "-sum_floor(-q*U-threshold/A)-sum_floor(q*L-threshold/A)-(B-A+1)",
        "coverage": "Every primitive fraction with denominator <= cap and quality < threshold lies in one of the leaf strips. Zero-count nodes reject all integer numerators in their whole denominator block.",
        "deduplication": "All q labels searched; gcd and reduced coordinates saved. Completeness claim concerns primitive labels. A nonprimitive label does not create a new rational candidate.",
        "dyadic_root_node_ids": roots, "nodes": nodes, "leaf_candidates": candidates,
        "retained": retained, "discarded_leaf_candidates": discarded, "unresolved": unresolved,
        "statistics": {
            "dyadic_blocks": len(roots), "node_count": len(nodes),
            "zero_count_blocks": sum(n["status"] == "zero_count_block_rejected" for n in nodes),
            "bisected_blocks": sum(n["status"] == "bisected" for n in nodes),
            "leaf_blocks": sum(n["status"] == "leaf_scored" for n in nodes),
            "leaf_denominators": [n["A"] for n in nodes if n["status"] == "leaf_scored"],
            "floor_sum_calls": 2 * len(nodes), "floor_sum_euclidean_steps": total_steps,
            "integer_numerator_candidates_scored": len(candidates),
            "retained_count": len(retained), "discarded_leaf_count": len(discarded),
            "unresolved_count": len(unresolved),
        },
    }


def validate_small(L, U, threshold):
    floor_cases = 0
    for n in range(11):
        for m in range(1, 9):
            for a in range(-9, 10):
                for b in range(-9, 10):
                    got, _ = floor_sum(n, m, a, b)
                    expected = sum((a * i + b) // m for i in range(n))
                    assert got == expected
                    floor_cases += 1
    strip_cases = 0
    for A in range(1, 9):
        for length in range(1, 9):
            B = A + length - 1
            for denominator in range(1, 8):
                l = F(3 * denominator + 1, denominator)
                u = l + F(1, denominator + 1)
                t = F(length, 3 * denominator + 1)
                got, _, _, _ = strip_count(A, B, l, u, t)
                assert got == brute_strip_count(A, B, l, u, t)
                strip_cases += 1
    # Directly check the actual pi strip count in each dyadic shell through 10000.
    A, pi_blocks, checked_q = 1, 0, 0
    while A <= 10000:
        B = min(2 * A - 1, 10000)
        got, _, _, _ = strip_count(A, B, L, U, threshold)
        assert got == brute_strip_count(A, B, L, U, threshold)
        pi_blocks += 1
        checked_q += B - A + 1
        A *= 2
    search = run_search(10000, L, U, threshold)
    brute_rows = []
    for q in range(1, 10001):
        brute_rows.extend(leaf_candidates(q, L, U, threshold))
    assert {(r["p"], r["q"], r["status"]) for r in search["leaf_candidates"]} == {
        (r["p"], r["q"], r["status"]) for r in brute_rows}
    return {
        "exhaustive_integer_floor_sum_cases": floor_cases,
        "small_rational_strip_cases": strip_cases,
        "pi_dyadic_strip_counts_checked_by_direct_summation": pi_blocks,
        "pi_denominators_directly_checked": checked_q,
        "pi_search_vs_direct_leaf_scan_cap": 10000,
        "all_checks_passed": True,
        "small_search_statistics": search["statistics"],
        "small_direct_candidate_pairs": [[r["p"], r["q"]] for r in brute_rows],
    }


def replay_tree(data):
    N, threshold = int(data["denominator_cap"]), unpack(data["threshold"])
    L, U = unpack(data["pi_lower"]), unpack(data["pi_upper"])
    nodes = data["nodes"]
    visited = set()

    def replay(node_id, A, B, parent):
        assert node_id not in visited
        visited.add(node_id)
        node = nodes[node_id]
        assert (int(node["A"]), int(node["B"]), node["parent"]) == (A, B, parent)
        count, lo_floor, neg_hi_floor, steps = strip_count(A, B, L, U, threshold)
        assert (count, lo_floor, neg_hi_floor, steps) == (
            int(node["count"]), int(node["lower_floor_sum"]),
            int(node["negative_upper_floor_sum"]), node["floor_sum_euclidean_steps"])
        if node["status"] == "zero_count_block_rejected":
            assert count == 0
        elif node["status"] == "leaf_scored":
            assert A == B and count > 0
            expected = leaf_candidates(A, L, U, threshold)
            got = [data["leaf_candidates"][i] for i in node["candidate_indices"]]
            assert got == expected
        else:
            assert node["status"] == "bisected" and A < B and count > 0
            middle = (A + B) // 2
            replay(node["children"][0], A, middle, node_id)
            replay(node["children"][1], middle + 1, B, node_id)

    A = 1
    for root_id in data["dyadic_root_node_ids"]:
        B = min(2 * A - 1, N)
        replay(root_id, A, B, None)
        A *= 2
    assert A > N and len(visited) == len(nodes)
    return {"all_tree_nodes_replayed": True, "nodes": len(nodes),
            "whole_denominator_interval_covered": True,
            "all_leaf_scores_rebuilt": True}


def postcheck_cf(data):
    # Called only after independent search, output write, and tree replay.
    cf_path = ROOT / "cf_all_denominators_1e100.json"
    if not cf_path.exists():
        return {"performed": False}
    cf = json.loads(cf_path.read_text())
    threshold = unpack(data["threshold"])
    certified_hits = {(r["p"], r["q"]) for r in data["retained"] if r["gcd"] == "1"}
    cf_hits = set()
    for row in cf["convergents"]:
        lo, hi = unpack(row["quality_lower"]), unpack(row["quality_upper"])
        if hi < threshold:
            cf_hits.add((row["p"], row["q"]))
        else:
            assert lo >= threshold
    assert certified_hits == cf_hits
    return {"performed": True, "search_input_used_cf": False,
            "primitive_hit_sets_identical": True,
            "primitive_hits": sorted([list(pair) for pair in certified_hits])}


if __name__ == "__main__":
    a, b = atan_bounds(5, 200)
    c, d = atan_bounds(239, 70)
    L, U, threshold = 16 * a - 4 * d, 16 * b - 4 * c, F(1, 250)
    small = validate_small(L, U, threshold)
    (ROOT / "batch_floor_small_validation.json").write_text(json.dumps(small, indent=2) + "\n")
    result = run_search(10 ** 100, L, U, threshold)
    out = ROOT / "batch_floor_search_1e100.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    replay = replay_tree(json.loads(out.read_text()))
    postcheck = postcheck_cf(result)
    verification = {"small_validation": small, "tree_replay": replay,
                    "continued_fraction_postcheck_only": postcheck}
    (ROOT / "batch_floor_verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    print(json.dumps({"result_path": str(out), "result_bytes": out.stat().st_size,
                      "statistics": result["statistics"], "small_validation": small,
                      "tree_replay": replay, "postcheck": postcheck}, indent=2))
