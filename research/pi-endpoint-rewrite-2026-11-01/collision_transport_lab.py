#!/usr/bin/env python3
"""Exact finite relaxation of the pi017 translated-row saving.

This audits a bound, not pi or the interpolation-rank theorem.  The
transverse indices are over-approximated: every a >= beta with wa <= H
is allowed, even if its actual scalar or coefficient determinant vanishes.
The certified collision constant is 1/6 < log(2)/4.
"""
from fractions import Fraction
from itertools import product
from heapq import heappop, heappush
from pathlib import Path
import json


def indices(weights, budget):
    return [a for a in product(*(range(budget // w + 1) for w in weights))
            if sum(w*x for w, x in zip(weights, a)) <= budget]


def weight(weights, a):
    return sum(w*x for w, x in zip(weights, a))


def supplies(weights, H, theta, v0, K):
    out = {}
    for beta in indices(weights, H):
        remainder = Fraction(H) - Fraction(weight(weights, beta), theta)
        if remainder <= 0:
            continue
        # s >= 0 and v0*s < remainder: number = ceil(remainder/v0).
        ratio = remainder / v0
        out[beta] = K * ((ratio.numerator + ratio.denominator - 1)//ratio.denominator)
    return out


def min_saving(weights, H, demand, nu=Fraction(2)):
    """Return min [sum n_a^2/6 + nu sum_rows w(a-beta)].

    Min-cost flow is exact: all costs multiplied by 6*denominator(nu).
    Convex target occupancy has marginal costs 2*k-1.  An edge from
    beta to a exists exactly when a dominates beta componentwise.
    """
    targets = indices(weights, H)
    betas = list(demand)
    M = sum(demand.values())
    scale = 6*nu.denominator
    collision_unit = nu.denominator
    scalar_unit = 6*nu.numerator
    source = 0
    beta_start = 1
    a_start = beta_start + len(betas)
    sink = a_start + len(targets)
    graph = [[] for _ in range(sink+1)]

    def add(u, v, cap, cost):
        edge = [v, len(graph[v]), cap, cost]
        graph[u].append(edge)
        graph[v].append([u, len(graph[u])-1, 0, -cost])
        return edge

    for j, beta in enumerate(betas):
        add(source, beta_start+j, demand[beta], 0)
        for k, a in enumerate(targets):
            if all(x >= y for x, y in zip(a, beta)):
                add(beta_start+j, a_start+k, M, scalar_unit*weight(weights, a))
    target_edges = []
    for k, a in enumerate(targets):
        edges = []
        for occupancy in range(1, M+1):
            edges.append(add(a_start+k, sink, 1, collision_unit*(2*occupancy-1)))
        target_edges.append(edges)

    potentials = [0]*len(graph)
    total_cost = 0
    for _ in range(M):
        dist = [None]*len(graph)
        previous = [None]*len(graph)
        dist[source] = 0
        queue = [(0, source)]
        while queue:
            cur, u = heappop(queue)
            if cur != dist[u]:
                continue
            for edge_index, edge in enumerate(graph[u]):
                v, _, cap, cost = edge
                if not cap:
                    continue
                reduced = cost + potentials[u] - potentials[v]
                assert reduced >= 0
                candidate = cur + reduced
                if dist[v] is None or candidate < dist[v]:
                    dist[v] = candidate
                    previous[v] = (u, edge_index)
                    heappush(queue, (candidate, v))
        assert dist[sink] is not None
        for v, d in enumerate(dist):
            if d is not None:
                potentials[v] += d
        v = sink
        while v != source:
            u, edge_index = previous[v]
            edge = graph[u][edge_index]
            total_cost += edge[3]
            edge[2] -= 1
            graph[v][edge[1]][2] += 1
            v = u
    beta_cost = scalar_unit * sum(count*weight(weights, beta) for beta, count in demand.items())
    result = Fraction(total_cost-beta_cost, scale)
    occupancies = {a: sum(edge[2] == 0 for edge in edges)
                   for a, edges in zip(targets, target_edges)}
    assert sum(occupancies.values()) == M
    direct = sum(Fraction(n*n, 6)+nu*n*weight(weights, a)
                 for a, n in occupancies.items())
    direct -= nu*sum(count*weight(weights, beta) for beta, count in demand.items())
    assert direct == result
    return result, occupancies


def brute(weights, H, demand):
    rows = [beta for beta, multiplicity in demand.items() for _ in range(multiplicity)]
    targets = indices(weights, H)
    allowed = [[a for a in targets if all(x >= y for x, y in zip(a, beta))] for beta in rows]
    best = None
    for assignment in product(*allowed):
        counts = {a: assignment.count(a) for a in set(assignment)}
        value = sum(Fraction(n*n, 6) for n in counts.values())
        value += 2*sum(weight(weights, a)-weight(weights, beta) for a, beta in zip(assignment, rows))
        best = value if best is None else min(best, value)
    return best


def serial(x):
    return f"{x.numerator}/{x.denominator}" if isinstance(x, Fraction) else x


def compare_budget(saving, arithmetic_budget, error_budget=Fraction(0)):
    """Positive means this finite envelope beats the supplied full budget.

    The caller, not this function, must establish arithmetic_budget for
    the same nonzero minor and error_budget for its analytic translation.
    """
    return saving-Fraction(arithmetic_budget)-Fraction(error_budget)


def main():
    checks = []
    for H, demand in [(2, {(0,): 2, (1,): 1}),
                      (3, {(0,): 2, (1,): 1, (2,): 1}),
                      (3, {(0, 0): 2, (1, 0): 1})]:
        weights = (1,) if len(next(iter(demand))) == 1 else (1, 2)
        value, _ = min_saving(weights, H, demand)
        assert value == brute(weights, H, demand)
        checks.append({"weights": weights, "H": H, "minimum": serial(value)})
    examples = []
    for H in [4, 6, 8, 10]:
        weights, theta, v0, K = (1, 3), Fraction(1, 2), Fraction(1), 2
        demand = supplies(weights, H, theta, v0, K)
        M = sum(demand.values())
        weighted_beta = sum(n*weight(weights, beta) for beta, n in demand.items())
        value, occupancy = min_saving(weights, H, demand)
        # Main arithmetic budget with selected-column cost replaced by MH.
        # A caller may replace this by its exact SAME-minor arithmetic budget.
        arithmetic_budget = M*H-weighted_beta
        examples.append({"H": H, "M": M, "weighted_beta_sum": weighted_beta,
                         "minimum_joint_saving": serial(value),
                         "normalized_joint_saving": serial(value/(M*H)),
                         "main_arithmetic_budget": arithmetic_budget,
                         "gap_before_all_errors": serial(compare_budget(value, arithmetic_budget)),
                         "occupancies": [{"a": a, "n": n} for a, n in occupancy.items() if n]})
    # Exact algebraic upper bound L<1/9 under the old two-alternative scheme,
    # using e<3 and eta^2(1-eta)^m <= 4/27.  log(2)/4<1/4.
    exact_bounds = {"L_upper_strict": "1/9", "cL_upper_strict": "1/36",
                    "cL_divided_by_exact_main_arithmetic_cost_upper_strict": "1/8"}
    payload = {"kind": "finite_relaxation_not_a_pi_proof", "collision_constant": "1/6",
               "bruteforce_crosschecks": checks, "examples": examples,
               "two_alternative_no_go_bounds": exact_bounds,
               "scope": ["All actual nonzero translated terms satisfy these relaxed index constraints.",
                         "Some relaxed assignments can be zero determinants or absent actual scalars.",
                         "Finite margins omit translation, analytic, arithmetic and limiting errors.",
                         "Interpolation separation and full rank are not asserted for toy weights.",
                         "Supply an independently valid same-minor arithmetic budget before comparing."]}
    destination = Path(__file__).with_name("collision_transport_results.json")
    destination.write_text(json.dumps(payload, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps({"bruteforce_crosschecks": len(checks), "examples": len(examples),
                      "destination": str(destination)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
