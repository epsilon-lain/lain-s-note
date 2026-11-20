#!/usr/bin/env python3
"""Exact integer checks for high-order formal-group interpolation weights."""
from math import comb


def weights(m):
    return {j: (-1) ** j * comb(m + 1, j - 1)
            for j in range(2, m + 3)}


def check_interpolation(m):
    lam = weights(m)
    # For f(k)=k^d, f(1)=1 for every d>=0.
    for d in range(m + 1):
        lhs = sum(w * (j ** d) for j, w in lam.items())
        assert lhs == 1, (m, d, lhs)
    assert sum(w * j for j, w in lam.items()) == 1


def multiplicative_series_coefficients(m):
    lam = weights(m)
    # ((1+T)^j - 1)/T has coefficient binom(j,r+1) at T^r.
    coeffs = []
    for r in range(m + 2):
        coeffs.append(sum(w * comb(j, r + 1)
                          for j, w in lam.items()
                          if r + 1 <= j))
    return coeffs


def main():
    rows = []
    for m in range(1, 31):
        check_interpolation(m)
        coeffs = multiplicative_series_coefficients(m)
        assert coeffs[0] == 1, (m, coeffs)
        assert all(coeffs[r] == 0 for r in range(1, m)), (m, coeffs)
        assert coeffs[m] == (-1) ** m, (m, coeffs)
        rows.append({
            "m": m,
            "weight_l1_norm": sum(abs(w) for w in weights(m).values()),
            "constant_coefficient": coeffs[0],
            "first_nonzero_error_coefficient": coeffs[m],
            "lower_error_coefficients_zero": True,
        })

    print("PASS: interpolation identities for m=1,...,30")
    print("PASS: multiplicative formal-group cancellation through order m-1")
    for row in rows[:8]:
        print(row)


if __name__ == "__main__":
    main()
