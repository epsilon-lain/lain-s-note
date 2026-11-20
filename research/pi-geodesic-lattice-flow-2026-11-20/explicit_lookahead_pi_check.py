#!/usr/bin/env python3
"""Verify the explicit partial-quotient future-lookahead bound for pi.

Requires mpmath. The theorem is symbolic; this script is a finite numerical
regression of its predicted lookahead against the stopping certificate.
"""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 100
ALPHA = mp.pi
COUNT = 45
CANDIDATES = 20


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


def fibonacci_bound(a_n2, a_n3):
    if a_n2 >= 2:
        return 3, None, 2
    A = a_n3
    fib = [0, 1, 1]
    while fib[-1] ** 2 * (2 * A + 3) < (A + 2) ** 2:
        fib.append(fib[-1] + fib[-2])
    k = len(fib) - 1
    return k, A, k - 1


def actual_stopping(rows, n):
    qn = rows[n]["q"]
    En = rows[n]["E"]
    U = mp.inf
    witness = None
    for m in range(n + 1, len(rows) - 1):
        qj = rows[m]["q"]
        Ej = rows[m]["E"]
        crossing = mp.mpf(qj * qj - qn * qn) / (En * En - Ej * Ej)
        if crossing < U:
            U, witness = crossing, m

        qnext = rows[m + 1]["q"]
        unseen_lower_bound = mp.mpf(qnext * qnext - qn * qn) / (En * En)
        if unseen_lower_bound >= U:
            return {
                "actual_lookahead": m - n,
                "future_witness_q": rows[witness]["q"],
                "future_threshold": mp.nstr(U, 30),
                "stop_lower_bound": mp.nstr(unseen_lower_bound, 30),
            }
    raise RuntimeError("Increase COUNT: stopping certificate not reached")


def main():
    rows = convergents(ALPHA, COUNT)
    results = []
    for n in range(1, CANDIDATES + 1):
        a_n2 = rows[n + 2]["a"]
        a_n3 = rows[n + 3]["a"]
        k, A, bound = fibonacci_bound(a_n2, a_n3)
        actual = actual_stopping(rows, n)
        assert actual["actual_lookahead"] <= bound, (n, bound, actual)
        results.append({
            "cf_index_zero_based": n,
            "candidate_q": rows[n]["q"],
            "a_n_plus_2": a_n2,
            "a_n_plus_3": A,
            "fibonacci_k": k,
            "proved_lookahead_bound": bound,
            **actual,
        })

    result = {
        "alpha": "pi",
        "precision_decimal_digits": mp.mp.dps,
        "convergents_generated": len(rows),
        "candidates_checked": CANDIDATES,
        "results": results,
        "all_actual_lookaheads_within_bound": True,
        "scope": "Finite numerical regression only; the lookahead theorem is proved independently.",
    }
    path = Path(__file__).with_name("explicit_lookahead_pi_check.json")
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "candidates_checked": CANDIDATES,
        "all_within_bound": True,
        "rows": [
            (r["candidate_q"], r["proved_lookahead_bound"], r["actual_lookahead"])
            for r in results
        ],
        "output": str(path),
    }, indent=2))


if __name__ == "__main__":
    main()
