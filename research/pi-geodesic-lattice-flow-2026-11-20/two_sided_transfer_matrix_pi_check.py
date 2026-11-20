#!/usr/bin/env python3
"""Check the two-sided transfer-matrix formulas against direct convergents.

Requires mpmath. This is a finite numerical regression, not a formal interval proof.
"""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 100
ALPHA = mp.pi
COUNT = 45
CANDIDATES = 25

I2 = ((mp.mpf(1), mp.mpf(0)), (mp.mpf(0), mp.mpf(1)))


def mat(a, signed=False):
    return ((mp.mpf(-a if signed else a), mp.mpf(1)),
            (mp.mpf(1), mp.mpf(0)))


def mmul(A, B):
    return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def mvec(A, v):
    return tuple(sum(A[i][j] * v[j] for j in range(2)) for i in range(2))


def product_right(indices, signed=False):
    P = I2
    for a in indices:
        P = mmul(P, mat(a, signed=signed))
    return P


def product_left(indices, signed=False):
    P = I2
    for a in indices:
        P = mmul(mat(a, signed=signed), P)
    return P


def convergents(alpha, count):
    z = mp.mpf(alpha)
    p0, p1 = 0, 1
    q0, q1 = 1, 0
    rows = []
    for n in range(count):
        a = int(mp.floor(z))
        p = a * p1 + p0
        q = a * q1 + q0
        E = abs(mp.mpf(q) * alpha - p)
        rows.append({"n": n, "a": a, "p": p, "q": q, "E": E})
        p0, p1 = p1, p
        q0, q1 = q1, q
        frac = z - a
        if frac == 0:
            break
        z = 1 / frac
    return rows


def lookahead_bound(rows, n):
    a2 = rows[n + 2]["a"]
    if a2 >= 2:
        return 2
    A = rows[n + 3]["a"]
    fib_prev, fib = 1, 1  # F_1,F_2
    k = 2
    B = mp.mpf(A + 2) / mp.sqrt(2 * A + 3)
    while fib < B:
        fib_prev, fib = fib, fib + fib_prev
        k += 1
    return k - 1


def matrix_ratios(rows, n, K):
    qn, En = mp.mpf(rows[n]["q"]), rows[n]["E"]
    r = mp.mpf(rows[n - 1]["q"]) / qn
    z = mp.mpf(rows[n + 1]["E"]) / En

    # Future denominator ratios Rhat_k.
    future_R = {}
    for k in range(1, K + 1):
        indices = [rows[j]["a"] for j in range(n + 1, n + k + 1)]
        P = product_left(indices, signed=False)
        future_R[k] = mvec(P, (mp.mpf(1), r))[0]

    # Future error ratios Zhat_k.
    future_Z = {1: z}
    for k in range(2, K + 1):
        indices = [rows[j]["a"] for j in range(n + 2, n + k + 1)]
        P = product_left(indices, signed=True)
        future_Z[k] = mvec(P, (z, mp.mpf(1)))[0]

    # Past denominator and error ratios.
    past = {}
    for h in range(1, n + 1):
        if h == 1:
            R = r
            Z = mp.mpf(rows[n - 1]["E"]) / En
        else:
            indices = [rows[j]["a"] for j in range(n - h + 2, n + 1)]
            Pq = product_right(indices, signed=True)
            R = mvec(Pq, (r, mp.mpf(1)))[0]
            Pe = product_right(indices, signed=False)
            Z = mvec(Pe, (mp.mpf(rows[n + 1]["a"]) + z, mp.mpf(1)))[0]
        past[h] = (R, Z)

    return future_R, future_Z, past, r, z


def main():
    rows = convergents(ALPHA, COUNT)
    checks = []
    for n in range(1, min(CANDIDATES, len(rows) - 4) + 1):
        K = lookahead_bound(rows, n)
        future_R, future_Z, past, r, z = matrix_ratios(rows, n, K)
        qn, En = mp.mpf(rows[n]["q"]), rows[n]["E"]
        rho2 = (qn / En) ** 2

        # Compare every transfer-matrix ratio with direct convergent ratios.
        max_error = mp.mpf(0)
        for k in range(1, K + 1):
            direct_R = mp.mpf(rows[n + k]["q"]) / qn
            direct_Z = mp.mpf(rows[n + k]["E"]) / En
            max_error = max(max_error, abs(future_R[k] - direct_R),
                            abs(future_Z[k] - direct_Z))
        for h, (R, Z) in past.items():
            direct_R = mp.mpf(rows[n - h]["q"]) / qn
            direct_Z = mp.mpf(rows[n - h]["E"]) / En
            max_error = max(max_error, abs(R - direct_R), abs(Z - direct_Z))

        P_norm = max(1 / rho2, *[
            (1 - R**2) / (Z**2 - 1) for R, Z in past.values()
        ])
        F_norm = min(
            (future_R[k]**2 - 1) / (1 - future_Z[k]**2)
            for k in range(1, K + 1)
        )

        P_direct = max(
            mp.mpf(1),
            *[
                mp.mpf(qn**2 - rows[i]["q"]**2)
                / (rows[i]["E"]**2 - En**2)
                for i in range(n)
            ],
        )
        F_direct = min(
            mp.mpf(rows[n + k]["q"]**2 - qn**2)
            / (En**2 - rows[n + k]["E"]**2)
            for k in range(1, K + 1)
        )
        assert abs(rho2 * P_norm - P_direct) < mp.mpf("1e-60")
        assert abs(rho2 * F_norm - F_direct) < mp.mpf("1e-60")
        assert max_error < mp.mpf("1e-60")

        checks.append({
            "q": rows[n]["q"],
            "lookahead": K,
            "max_matrix_ratio_error": mp.nstr(max_error, 8),
            "past_threshold": mp.nstr(P_direct, 25),
            "future_threshold": mp.nstr(F_direct, 25),
            "active_by_selector": bool(P_direct < F_direct),
        })

    result = {
        "alpha": "pi",
        "precision_decimal_digits": mp.mp.dps,
        "convergents_generated": len(rows),
        "candidates_checked": len(checks),
        "checks": checks,
        "all_transfer_formulas_match_direct": True,
        "scope": "Finite high-precision regression. Symbolic transfer identities and lookahead theorem are proved separately.",
    }
    out = Path(__file__).with_name("two_sided_transfer_matrix_pi_check.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "candidates_checked": len(checks),
        "all_transfer_formulas_match_direct": True,
        "output": str(out),
    }, indent=2))


if __name__ == "__main__":
    main()
