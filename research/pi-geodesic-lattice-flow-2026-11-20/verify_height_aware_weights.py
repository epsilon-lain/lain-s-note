#!/usr/bin/env python3
"""Exact verification of height-aware integer formal-group weights.

All polynomial and moment calculations use Python integers only.
The candidates are improvements over the binomial weights, not claimed
global optima.
"""
from math import comb


def trim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a, b):
    out = [0] * max(len(a), len(b))
    for i, v in enumerate(a):
        out[i] += v
    for i, v in enumerate(b):
        out[i] += v
    return trim(out)


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def power(a, n):
    out = [1]
    for _ in range(n):
        out = mul(out, a)
    return out


def shift_x_minus_one_power(n):
    # (x-1)^n, coefficients in ascending order.
    return [(-1) ** (n-k) * comb(n, k) for k in range(n + 1)]


def divide_by_x(poly):
    assert poly[0] == 0
    return trim(poly[1:])


def weights_from_Q(m, Q):
    numerator = add([1], mul(shift_x_minus_one_power(m + 1), Q))
    L = divide_by_x(numerator)
    return {j + 2: c for j, c in enumerate(L) if c}


def binomial_weights(m):
    return {j: (-1) ** j * comb(m + 1, j - 1)
            for j in range(2, m + 3)}


def moment_check(weights, m):
    return [
        sum(lam * (j ** r) for j, lam in weights.items())
        for r in range(m + 1)
    ]


def costs(weights):
    c0 = sum(abs(lam) for lam in weights.values())
    c2 = sum(abs(lam) * j * j for j, lam in weights.items())
    return c0, c2


def main():
    # Q coefficients in ascending powers of x.
    candidates_Q = {
        4: [1, 1],  # x+1
        5: [-1, -2, -2, -1],  # -(x+1)(x^2+x+1)
        6: [1, 3, 4, 3, 1],  # (x+1)^2(x^2+x+1)
    }
    results = []
    for m, Q in candidates_Q.items():
        alt = weights_from_Q(m, Q)
        base = binomial_weights(m)
        moments = moment_check(alt, m)
        assert moments == [1] * (m + 1), (m, alt, moments)
        c0_alt, c2_alt = costs(alt)
        c0_base, c2_base = costs(base)
        results.append({
            "m": m,
            "alternative_weights": alt,
            "moments": moments,
            "alternative_l1": c0_alt,
            "binomial_l1": c0_base,
            "alternative_quadratic_cost": c2_alt,
            "binomial_quadratic_cost": c2_base,
            "quadratic_cost_reduction_percent":
                100 * (c2_base - c2_alt) / c2_base,
        })

    expected = {
        4: {2: 4, 3: -5, 5: 5, 6: -4, 7: 1},
        5: {2: 4, 3: -5, 4: 1, 5: 1, 6: 1, 7: 1, 8: -5, 9: 4, 10: -1},
        6: {2: 4, 3: -4, 4: -3, 5: 6, 8: -6, 9: 3, 10: 4, 11: -4, 12: 1},
    }
    for row in results:
        assert row["alternative_weights"] == expected[row["m"]]

    print("PASS: all exact moment identities")
    print("PASS: all alternative weights are integers")
    for row in results:
        print(row)


if __name__ == "__main__":
    main()
