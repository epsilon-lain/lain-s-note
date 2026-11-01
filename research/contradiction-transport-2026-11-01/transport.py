"""Exact finite convex transport with coordinatewise upward moves.

Standard library only. This optimizes a finite relaxation, not a pi determinant.
"""
from collections import Counter
from fractions import Fraction as F
from heapq import heappop, heappush
from itertools import product
from math import lcm


def leq(origin, destination):
    return all(x <= y for x, y in zip(origin, destination))


def dot(weights, point):
    return sum((x * y for x, y in zip(weights, point)), F(0))


def energy(counts, c, weights):
    return sum((c * n * n + 2 * dot(weights, a) * n
                for a, n in counts.items()), F(0))


def solve(supply, destinations, c=F(1, 6), weights=None):
    """Return exact minimum, destination counts and nonzero transport flows.

    `supply` maps nonnegative integer tuples to nonnegative integer counts.
    Destinations are an explicitly finite set. No unbounded-domain claim is made.
    Convex destination cost is expanded into unit edges; integer min-cost flow
    uses Dijkstra potentials, with all rational costs scaled to integers.
    """
    c = F(c)
    destinations = sorted(set(destinations))
    origins = sorted(a for a, n in supply.items() if n)
    if not origins:
        return F(0), {}, []
    dimension = len(origins[0])
    weights = tuple(map(F, weights or [1] * dimension))
    if dimension == 0 or c <= 0 or len(weights) != dimension or any(w <= 0 for w in weights):
        raise ValueError('c and all weights must be positive')
    for a in origins + destinations:
        if len(a) != dimension or any(not isinstance(x, int) or x < 0 for x in a):
            raise ValueError('points must have the same positive dimension')
    if any(not isinstance(n, int) or n < 0 for n in supply.values()):
        raise ValueError('supplies must be nonnegative integers')
    mass = sum(supply.values())
    scale = lcm(c.denominator, *(w.denominator for w in weights))
    source, sink = 0, 1 + len(origins) + len(destinations)
    graph = [[] for _ in range(sink + 1)]

    def edge(a, b, cap, cost):
        forward = [b, len(graph[b]), cap, cost]
        reverse = [a, len(graph[a]), 0, -cost]
        graph[a].append(forward)
        graph[b].append(reverse)
        return forward

    flow_edges = []
    for i, beta in enumerate(origins, 1):
        edge(source, i, supply[beta], 0)
        for j, a in enumerate(destinations, 1 + len(origins)):
            if leq(beta, a):
                e = edge(i, j, mass, 0)
                flow_edges.append((beta, a, e))
    destination_edges = []
    for j, a in enumerate(destinations, 1 + len(origins)):
        edges = []
        for k in range(1, mass + 1):
            cost = scale * (c * (2 * k - 1) + 2 * dot(weights, a))
            assert cost.denominator == 1
            edges.append(edge(j, sink, 1, int(cost)))
        destination_edges.append((a, edges))

    potentials = [0] * len(graph)
    total_cost = 0
    for _ in range(mass):
        distances = [None] * len(graph)
        previous = [None] * len(graph)
        distances[source] = 0
        queue = [(0, source)]
        while queue:
            d, a = heappop(queue)
            if distances[a] != d:
                continue
            for idx, e in enumerate(graph[a]):
                b, _, cap, cost = e
                if cap <= 0:
                    continue
                reduced = cost + potentials[a] - potentials[b]
                assert reduced >= 0
                candidate = d + reduced
                if distances[b] is None or candidate < distances[b]:
                    distances[b] = candidate
                    previous[b] = (a, idx)
                    heappush(queue, (candidate, b))
        if distances[sink] is None:
            raise ValueError('finite destination set cannot accept all supply')
        for a, d in enumerate(distances):
            if d is not None:
                potentials[a] += d
        b = sink
        while b != source:
            a, idx = previous[b]
            e = graph[a][idx]
            total_cost += e[3]
            e[2] -= 1
            graph[b][e[1]][2] += 1
            b = a
    counts = {a: sum(1 - e[2] for e in edges)
              for a, edges in destination_edges}
    counts = {a: n for a, n in counts.items() if n}
    flows = [(beta, a, mass - e[2]) for beta, a, e in flow_edges if e[2] < mass]
    assert sum(counts.values()) == mass
    assert energy(counts, c, weights) == F(total_cost, scale)
    assert Counter({b: sum(n for origin, _, n in flows if origin == b)
                    for b in origins}) == Counter(supply)
    return F(total_cost, scale), counts, flows


def exhaustive(supply, destinations, c, weights):
    """Independent tiny-instance reference: enumerate individual assignments."""
    choices = []
    for beta, n in supply.items():
        choices.extend([[a for a in destinations if leq(beta, a)]] * n)
    return min(energy(Counter(assignment), c, weights)
               for assignment in product(*choices))


def identity_terms(supply, counts, flows, c, weights):
    """Exact quadratic decomposition from README theorem 1."""
    points = set(supply) | set(counts)
    squared = c * sum((counts.get(a, 0) - supply.get(a, 0)) ** 2 for a in points)
    potential = lambda a: c * supply.get(a, 0) + dot(weights, a)
    drift = 2 * sum(n * (potential(a) - potential(beta)) for beta, a, n in flows)
    gap = energy(counts, c, weights) - energy(supply, c, weights)
    assert gap == squared + drift
    return gap, squared, drift
