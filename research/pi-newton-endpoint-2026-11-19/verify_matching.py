"""Check the sharp diagonal bound by two independent exact computations."""
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb, gcd, prod


def matching_bound(caps, weights, degree):
    """Choose an unused slot, then try leaving it free or pairing it."""
    @lru_cache(None)
    def solve(remaining, edge_budget):
        if edge_budget == 0 or not any(remaining):
            return 0
        first = next(i for i, count in enumerate(remaining) if count)
        fewer = list(remaining)
        fewer[first] -= 1
        best = solve(tuple(fewer), edge_budget)
        for other in range(first + 1, len(caps)):
            if fewer[other]:
                fewer[other] -= 1
                value = min(weights[first], weights[other])
                best = max(best, value + solve(tuple(fewer), edge_budget - 1))
                fewer[other] += 1
        return best

    return solve(caps, min(degree, sum(caps) - degree))


def add_row(pivots, row):
    """Integer echelon reduction, exact over Q; no modular rank assumptions."""
    for column in range(len(row)):
        if row[column] == 0:
            continue
        if column not in pivots:
            pivots[column] = row
            return
        previous = pivots[column]
        a, b = previous[column], row[column]
        row = [a * x - b * y for x, y in zip(row, previous)]
        common = gcd(*row)
        if common > 1:
            row = [v // common for v in row]


def jet_rank_bound(caps, weights, degree):
    """First jet-weight layer at which a homogeneous polynomial is forced zero."""
    exponent_box = list(product(*(range(cap + 1) for cap in caps)))
    columns = [alpha for alpha in exponent_box if sum(alpha) == degree]
    jets = [u for u in exponent_box if sum(u) <= degree]
    jets.sort(key=lambda u: sum(w * k for w, k in zip(weights, u)))
    pivots = {}
    for u in jets:
        # [e^u] (t+e)^alpha / t^(degree-sum(u)), identically in t.
        row = [prod(comb(a, k) for a, k in zip(alpha, u)) for alpha in columns]
        add_row(pivots, row)
        if len(pivots) == len(columns):
            return sum(w * k for w, k in zip(weights, u))
    raise AssertionError("all monomial Taylor jets should have full rank")


def pi_upper_bound():
    # Alternating arctan sums and Machin's formula give a rational upper bound.
    def partial(inverse, terms):
        return sum(Fraction((-1) ** k, (2 * k + 1) * inverse ** (2 * k + 1))
                   for k in range(terms))
    return 16 * partial(5, 13) - 4 * partial(239, 12)


def same_side_checks():
    upper = pi_upper_bound()
    families = [[(22, 7), (355, 113), (104348, 33215)],
                [(22, 7), (45, 14), (67, 21)]]
    caps = (4, 3, 2)
    total = sum(caps)
    cases = 0
    for family in families:
        denominators = [q for _, q in family]
        errors = [Fraction(p, q) - upper for p, q in family]
        assert all(0 < error for error in errors)
        assert all(Fraction(p, q) < 4 for p, q in family)  # π > 3
        common = gcd(*denominators)
        raw = prod(q ** d for q, d in zip(denominators, caps))
        pairs = [(0, 1), (0, 2), (1, 2)]
        for numbers in product(*(range(min(caps[i], caps[j]) + 1)
                                 for i, j in pairs)):
            usage = [0, 0, 0]
            edge_factor = Fraction(1)
            for (i, j), number in zip(pairs, numbers):
                usage[i] += number
                usage[j] += number
                edge_factor *= max(errors[i], errors[j]) ** number
            if any(n > d for n, d in zip(usage, caps)):
                continue
            edges = sum(numbers)
            for degree in range(edges, total - edges + 1):
                cleared = Fraction(raw, common ** (total - degree)) * edge_factor
                assert cleared >= common ** (degree - edges)
                cases += 1
    return cases


def verify():
    configurations = []
    for count in (2, 3, 4):
        for caps in product(range(1, 3), repeat=count):
            for weights in ((1,) * count, tuple(range(1, count + 1)),
                            tuple(2 ** (count - 1 - i) for i in range(count))):
                configurations.append((caps, weights))
    configurations.extend([((4, 3, 2), (1, 4, 11)),
                           ((2, 0, 3), (7, 7, 2)),
                           ((0,), (3,)), ((3,), (3,)),
                           ((32, 10), (1, 4))])
    checked = 0
    for caps, weights in configurations:
        values = []
        for degree in range(sum(caps) + 1):
            bound = matching_bound(caps, weights, degree)
            actual = jet_rank_bound(caps, weights, degree)
            if bound != actual:
                raise AssertionError((caps, weights, degree, bound, actual))
            values.append(actual)
            checked += 1
        if caps in ((4, 3, 2), (32, 10)):
            print(f"caps={caps}, weights={weights}, maximum={max(values)}")
    print(f"Matched {checked} homogeneous jet spaces in "
          f"{len(configurations)} capacity/weight configurations.")
    print(f"Verified {same_side_checks()} same-side bounds with rational π enclosures.")


if __name__ == "__main__":
    verify()
