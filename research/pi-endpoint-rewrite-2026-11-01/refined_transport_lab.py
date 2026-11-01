#!/usr/bin/env python3
"""Exact finite relaxation with j=0, selected-column support and log gains.

All supplied gains are rational inputs. The caller must independently
prove that they are valid lower bounds in the analytic determinant.
The examples are finite rational toys, not measured or certified pi data.
"""
from fractions import Fraction
from itertools import product
from math import lcm
from pathlib import Path
import json


def fraction_text(x):
    x = Fraction(x)
    return f"{x.numerator}/{x.denominator}"


def target_domain(candidates, selected_alphas, enforce_selected_domain):
    return [a for a in candidates if not enforce_selected_domain
            or any(all(x <= y for x, y in zip(a, alpha)) for alpha in selected_alphas)]


def allowed_targets(j, beta, targets, enforce_j_zero):
    return [a for a in targets if all(x >= y for x, y in zip(a, beta))
            and (not enforce_j_zero or j != 0 or a == beta)]


def min_saving(demand, candidates, selected_alphas, lambda_by_j,
               d=Fraction(0), c0=Fraction(1, 6),
               enforce_j_zero=True, enforce_selected_domain=True):
    """Minimize c0*sum(n_a^2)+sum lambda_ji*(a_i-beta_i)+d*sum|a|.

    demand maps (j,beta_tuple) to a nonnegative integer multiplicity.
    candidates is a finite list of nonnegative integer multi-indices.
    selected_alphas contains the multi-indices of the SAME minor's columns.
    lambda_by_j[j] is a vector of rational certified lower bounds for
    -ln|epsilon_(j,i)|; lambda_by_j[0] is unused when j=0 is enforced.
    """
    c0, d = Fraction(c0), Fraction(d)
    assert c0 > 0 and d >= 0
    candidates = sorted(set(candidates))
    dimension = len(candidates[0])
    assert all(len(a) == dimension and all(isinstance(x, int) and x >= 0 for x in a)
               for a in candidates)
    assert all(len(alpha) == dimension for alpha in selected_alphas)
    demand = {group: count for group, count in demand.items() if count}
    assert all(isinstance(count, int) and count > 0 for count in demand.values())
    groups = list(demand)
    gains = {j: tuple(Fraction(x) for x in lambda_by_j[j]) for j, _ in groups}
    assert all(len(gain) == dimension for gain in gains.values())
    targets = target_domain(candidates, selected_alphas, enforce_selected_domain)
    options = {group: allowed_targets(*group, targets, enforce_j_zero) for group in groups}
    if any(not values for values in options.values()):
        raise ValueError("A supplied row group has no admissible target.")
    M = sum(demand.values())
    scale = lcm(c0.denominator, d.denominator,
                *(value.denominator for gain in gains.values() for value in gain))
    source = 0
    group_start = 1
    target_start = group_start + len(groups)
    sink = target_start + len(targets)
    graph = [[] for _ in range(sink+1)]

    def add(u, v, capacity, cost):
        cost = Fraction(cost)*scale
        assert cost.denominator == 1
        edge = [v, len(graph[v]), capacity, int(cost)]
        graph[u].append(edge)
        graph[v].append([u, len(graph[u])-1, 0, -int(cost)])
        return edge

    allocations = {}
    for index, (j, beta) in enumerate(groups):
        node = group_start+index
        add(source, node, demand[j, beta], 0)
        for target_index, a in enumerate(targets):
            if a not in options[j, beta]:
                continue
            linear = sum(gain*(x-y) for gain, x, y in zip(gains[j], a, beta)) + d*sum(a)
            allocations[j, beta, a] = (add(node, target_start+target_index, M, linear), M)
    occupancy_edges = {}
    for index, a in enumerate(targets):
        occupancy_edges[a] = [add(target_start+index, sink, 1, c0*(2*k-1))
                              for k in range(1, M+1)]

    total_scaled = 0
    # Bellman-Ford shortest augmenting paths handle even negative certified
    # lambda lower bounds. All capacities and scaled costs remain integers.
    for _ in range(M):
        distance = [None]*len(graph)
        previous = [None]*len(graph)
        distance[source] = 0
        for _pass in range(len(graph)-1):
            changed = False
            for u, edges in enumerate(graph):
                if distance[u] is None:
                    continue
                for edge_index, (v, _, capacity, cost) in enumerate(edges):
                    if not capacity:
                        continue
                    candidate = distance[u]+cost
                    if distance[v] is None or candidate < distance[v]:
                        distance[v] = candidate
                        previous[v] = (u, edge_index)
                        changed = True
            if not changed:
                break
        assert distance[sink] is not None
        for u, edges in enumerate(graph):
            if distance[u] is not None:
                assert all(not capacity or distance[v] <= distance[u]+cost
                           for v, _, capacity, cost in edges)
        v = sink
        while v != source:
            u, edge_index = previous[v]
            edge = graph[u][edge_index]
            total_scaled += edge[3]
            edge[2] -= 1
            graph[v][edge[1]][2] += 1
            v = u

    occupancy = {a: sum(edge[2] == 0 for edge in edges) for a, edges in occupancy_edges.items()}
    allocation = {key: capacity-edge[2] for key, (edge, capacity) in allocations.items()}
    assert sum(occupancy.values()) == M
    for j, beta in groups:
        assert sum(n for (jj, bb, _), n in allocation.items() if (jj, bb) == (j, beta)) == demand[j, beta]
    for a in targets:
        assert sum(n for (_, _, aa), n in allocation.items() if aa == a) == occupancy[a]
    direct = c0*sum(n*n for n in occupancy.values())
    direct += sum(n*(sum(gain*(x-y) for gain, x, y in zip(gains[j], a, beta))+d*sum(a))
                  for (j, beta, a), n in allocation.items())
    result = Fraction(total_scaled, scale)
    assert direct == result
    return result, occupancy, allocation


def brute_saving(demand, candidates, selected_alphas, lambda_by_j,
                 d=Fraction(0), c0=Fraction(1, 6),
                 enforce_j_zero=True, enforce_selected_domain=True):
    """Independent exhaustive full row assignments for tiny instances."""
    targets = target_domain(candidates, selected_alphas, enforce_selected_domain)
    rows = [group for group, count in demand.items() for _ in range(count)]
    options = [allowed_targets(*group, targets, enforce_j_zero) for group in rows]
    best = None
    count = 0
    for assignment in product(*options):
        occupancy = {a: assignment.count(a) for a in set(assignment)}
        score = Fraction(c0)*sum(n*n for n in occupancy.values())
        for (j, beta), a in zip(rows, assignment):
            score += sum(Fraction(gain)*(x-y) for gain, x, y in zip(lambda_by_j[j], a, beta))
            score += Fraction(d)*sum(a)
        best = score if best is None else min(best, score)
        count += 1
    return best, count


def compare_budget(saving, same_minor_arithmetic_budget,
                   selected_column_growth, test_cost, residual_cost=Fraction(0)):
    """Return saving minus all input costs, using exact rational arithmetic.

    Every cost must be independently justified for the SAME nonzero minor.
    A finite positive toy margin does not establish its applicability to pi.
    """
    return Fraction(saving)-sum(map(Fraction, [same_minor_arithmetic_budget,
        selected_column_growth, test_cost, residual_cost]))


def main():
    demand = {(j, (0,)): 2 for j in range(3)}
    candidates, selected = [(0,), (1,), (2,)], [(1,)]
    gains = {0: (Fraction(0),), 1: (Fraction(0),), 2: (Fraction(1, 2),)}
    d, c0 = Fraction(1, 60), Fraction(1, 6)
    variants = [
        ("all_information", gains, True, True),
        ("delete_j_zero_constraint", gains, False, True),
        ("delete_selected_column_domain", gains, True, False),
        ("replace_heterogeneous_lambda_by_coordinate_minimum",
         {j: (Fraction(0),) for j in range(3)}, True, True)]
    results = []
    for name, gain, j_zero, selected_domain in variants:
        value, occupancy, allocation = min_saving(demand, candidates, selected, gain,
            d=d, c0=c0, enforce_j_zero=j_zero, enforce_selected_domain=selected_domain)
        brute, number = brute_saving(demand, candidates, selected, gain,
            d=d, c0=c0, enforce_j_zero=j_zero, enforce_selected_domain=selected_domain)
        assert value == brute
        results.append({"name": name, "minimum_saving": fraction_text(value),
            "bruteforce_assignments": number,
            "occupancies": [{"a": a, "n": n} for a, n in occupancy.items() if n],
            "allocation": [{"j": j, "beta": beta, "a": a, "n": n}
                           for (j, beta, a), n in allocation.items() if n]})
    baseline = Fraction(results[0]["minimum_saving"])
    for result in results:
        result["saving_lost_against_all_information"] = fraction_text(baseline-Fraction(result["minimum_saving"]))
    assert all(Fraction(result["saving_lost_against_all_information"]) > 0 for result in results[1:])

    # A second, genuinely two-coordinate test includes a nonzero beta and
    # a negative lambda lower bound, and exercises integer cost scaling.
    second = dict(demand={(0, (0, 0)): 1, (1, (1, 0)): 1, (2, (0, 0)): 1},
                  candidates=[(0, 0), (1, 0), (0, 1), (1, 1), (2, 0)],
                  selected_alphas=[(1, 1), (2, 0)],
                  lambda_by_j={0: (Fraction(0), Fraction(0)),
                               1: (Fraction(2, 7), Fraction(-1, 11)),
                               2: (Fraction(3, 5), Fraction(1, 13))},
                  d=Fraction(1, 17), c0=Fraction(1, 6))
    value, _, _ = min_saving(**second)
    brute, number = brute_saving(**second)
    assert value == brute
    comparison = compare_budget(baseline, Fraction(1, 3), Fraction(2, 5), Fraction(1, 7), Fraction(1, 11))
    assert comparison == baseline-Fraction(1, 3)-Fraction(2, 5)-Fraction(1, 7)-Fraction(1, 11)
    payload = {"kind": "rational_finite_transport_toys_not_pi_data", "c0": fraction_text(c0),
               "d": fraction_text(d), "selected_alphas": selected,
               "candidate_targets": candidates,
               "demand": [{"j": j, "beta": beta, "multiplicity": n} for (j, beta), n in demand.items()],
               "lambda_by_j": {str(j): [fraction_text(x) for x in gain] for j, gain in gains.items()},
               "ablation_results": results,
               "two_coordinate_crosscheck": {"minimum_saving": fraction_text(value), "enumerated_assignments": number},
               "same_minor_budget_interface_arithmetic_crosscheck": fraction_text(comparison),
               "bruteforce_crosschecks_passed": 5,
               "boundaries": ["Every lambda,d,c0 is an input requiring independent analytic certification.",
                   "All supplies and targets in this file are finite toys, not data from pi approximants.",
                   "Selected-domain admissibility and the j=0 zero rule are retained separately.",
                   "Actual coefficient determinants can still vanish within the relaxed feasible set.",
                   "A positive toy comparison does not certify a nonzero minor, analytic bounds or pi property."]}
    destination = Path(__file__).with_name("refined_transport_results.json")
    destination.write_text(json.dumps(payload, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps({"ablation_variants": len(results), "bruteforce_crosschecks": 5,
                      "destination": str(destination)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
