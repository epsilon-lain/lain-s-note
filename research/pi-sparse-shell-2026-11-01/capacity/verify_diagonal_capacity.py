"""Exact finite checks for DIAGONAL_QUOTIENT_CAPACITY.md."""
from fractions import Fraction
from math import comb
from itertools import product
import json
from pathlib import Path


def ceil_fraction(x):
    return -((-x.numerator) // x.denominator)


def strict_count(budget, unit):
    return max(0, ceil_fraction(budget / unit))


def triangle_count(H, xweight, yweight, strict):
    result = 0
    xmax = int(H // xweight)
    for x in range(xmax + 1):
        budget = H - xweight * x
        if strict:
            result += strict_count(budget, yweight)
        elif budget >= 0:
            result += int(budget // yweight) + 1
    return result


def weak_compositions(k, m):
    if m == 1:
        yield (k,)
    else:
        for first in range(k + 1):
            for rest in weak_compositions(k - first, m - 1):
                yield (first,) + rest


H = Fraction(100)
K = 2
theta = Fraction(1, 4)
w0 = Fraction(1)
v0 = Fraction(1, 5)
w = (Fraction(1), Fraction(2))
wmin, wmax = min(w), max(w)

L = triangle_count(H, wmax / theta, v0, strict=True)
D = triangle_count(H, w0, wmin, strict=False)
assert (L, D, K * L) == (3380, 5151, 6760)
assert K * L > D

s0_rows = K * strict_count(theta * H, wmax)
s0_capacity = int(H // wmin) + 1
assert (s0_rows, s0_capacity) == (26, 101)
assert s0_rows <= s0_capacity

type_checks = []
for a in product(*(range(int(theta * H // wi) + 1) for wi in w)):
    wa = sum(wi * ai for wi, ai in zip(w, a))
    if wa >= theta * H:
        continue
    n = K * strict_count(H - wa / theta, v0)
    capacity = triangle_count(H - wa, w0, wmin, strict=False)
    assert n <= capacity
    type_checks.append({"a": a, "n": n, "capacity": capacity})

zero_type = next(item for item in type_checks if item["a"] == (0, 0))
assert (zero_type["n"], zero_type["capacity"]) == (1000, 5151)
assert K * theta ** len(w) < 1
assert K * w0 * theta ** len(w) / v0 < 1

# For a monomial Y^h X^alpha, the common time-dependent factor in
# sum_{|beta|=k} T_{j,s,beta} is multiplied by the binomial identity
# below. This verifies the complete row-sum identity, since that factor
# does not depend on the composition beta of k.
row_identity_checks = 0
for m in (1, 2, 3):
    for alpha in product(range(5), repeat=m):
        for k in range(sum(alpha) + 2):
            lhs = 0
            for beta in weak_compositions(k, m):
                term = 1
                for ai, bi in zip(alpha, beta):
                    term *= comb(ai, bi) if bi <= ai else 0
                lhs += term
            rhs = comb(sum(alpha), k) if k <= sum(alpha) else 0
            assert lhs == rhs
            row_identity_checks += 1

result = {
    "all_checks_passed": True,
    "L_H": L,
    "D_H": D,
    "K_L_H": K * L,
    "forced_rank_defect_at_least": K * L - D,
    "single_time_rows": s0_rows,
    "single_time_capacity": s0_capacity,
    "individual_types_checked": len(type_checks),
    "all_individual_type_capacities_pass": True,
    "zero_type": zero_type,
    "row_sum_binomial_identities_checked": row_identity_checks,
    "source_numerical_volume_conditions_pass": True,
    "source_separated_weight_thresholds_asserted": False,
}
output = Path(__file__).with_name("diagonal_capacity_verification.json")
output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
