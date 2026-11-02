"""Independent exact check of the coordinatewise relaxed-witness mechanism.

This constructs only multiplicities and column sets. It does not construct
time labels, a transport minor, a nonzero period minor, or an approximation
tuple. The limiting critical-family assertion remains the written proof.
"""

from fractions import Fraction as F
from math import ceil
from pathlib import Path
import json


def main():
    weights = (1, 2)
    height = 60
    theta = F(1, 4)
    v0 = F(1, 8)
    target_upper = F(7, 25)
    types = [(a, b) for a in range(height + 1)
             for b in range(height // 2 + 1) if a + 2 * b <= height]
    columns = [(h, a, b) for a, b in types
               for h in range(height - a - 2 * b + 1)]
    source = [(a, b, ceil(F(height - 4 * a - 8 * b) / v0))
              for a, b in types if 4 * a + 8 * b < height]
    mass = sum(n for a, b, n in source)
    assert v0 == 2 * theta ** 2
    assert mass <= len(columns)
    selected = sorted(columns, key=lambda c: (c[1] + 2 * c[2], c))[:mass]

    # m=2 gives the low-mean radial cutoff A=(m+2)t/m=2t.
    cutoff = 2 * target_upper
    values = [max(F(0), cutoff - F(a + 2 * b, height)) for a, b in types]
    normalization = F(mass, sum(values))
    real_profile = [normalization * value for value in values]
    profile = [int(value) for value in real_profile]
    remainder = mass - sum(profile)
    assert 0 <= remainder < len(types)
    cheapest = sorted(range(len(types)),
                      key=lambda i: (types[i][0] + 2 * types[i][1], i))
    for i in cheapest[:remainder]:
        profile[i] += 1
    assert sum(profile) == mass
    assert all(abs(F(n) - r) <= 1 for n, r in zip(profile, real_profile))

    source_totals = [sum(n * (a, b)[i] for a, b, n in source) for i in range(2)]
    target_totals = [sum(n * a[i] for a, n in zip(types, profile)) for i in range(2)]
    column_totals = [sum(c[i + 1] for c in selected) for i in range(2)]
    gamma = [target_totals[i] - source_totals[i] for i in range(2)]
    residual = [column_totals[i] - target_totals[i] for i in range(2)]
    gap = [column_totals[i] - source_totals[i] for i in range(2)]
    assert all(x > 0 for x in gamma + residual)
    assert gap == [gamma[i] + residual[i] for i in range(2)]

    def mean(totals):
        return F(sum(w * a for w, a in zip(weights, totals)), mass * height)

    assert mean(source_totals) < mean(target_totals) < mean(column_totals)
    assert mean(target_totals) <= target_upper
    result = {
        "status": "passed", "m": 2, "H": height, "weights": weights,
        "theta": str(theta), "v0": str(v0), "K": 1, "w0": 1,
        "source_rows": mass, "all_raw_columns": len(columns),
        "source_mean": str(mean(source_totals)),
        "target_mean": str(mean(target_totals)),
        "column_mean": str(mean(column_totals)),
        "Gamma": gamma, "E": residual, "D": gap,
        "collision_over_MH": str(F(sum(n * n for n in profile), mass * height)),
        "scope": "Exact feasible mass, mean, coordinate degrees, and selected raw columns "
                 "only. No row/time assignment or determinant nonvanishing is asserted.",
    }
    Path(__file__).with_name("independent_relaxed_witness.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
