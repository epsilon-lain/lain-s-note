#!/usr/bin/env python3
"""Numerical diagnostic for the exact local partial-quotient threshold."""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 100
ALPHA = mp.pi
QMAX = 1_000_000


def convergents(alpha, qmax):
    z = mp.mpf(alpha)
    p0, p1 = 0, 1
    q0, q1 = 1, 0
    out = []
    n = 0
    while True:
        a = int(mp.floor(z))
        p = a * p1 + p0
        q = a * q1 + q0
        if q > qmax:
            break
        out.append({"n": n, "a": a, "p": p, "q": q})
        p0, p1 = p1, p
        q0, q1 = q1, q
        frac = z - a
        if frac == 0:
            break
        z = 1 / frac
        n += 1
    return out


def main():
    c = convergents(ALPHA, QMAX)
    rows = []
    for n in range(1, len(c) - 2):
        if c[n + 1]["a"] != 1:
            continue
        beta = mp.mpf(c[n]["q"]) / c[n - 1]["q"]
        alpha_forward = mp.mpf(0)
        # Recompute the complete quotient alpha_{n+2} from the CF tail.
        # It equals a_{n+2} + 1/(a_{n+3}+...); derive it backwards
        # from a high-precision pi CF tail built in the convergent routine.
        # Store the tail by restarting the standard CF recurrence.
        z = mp.mpf(ALPHA)
        for _ in range(n + 2):
            z = 1 / (z - mp.floor(z))
        alpha_forward = z
        threshold = (beta + 2) / (beta - 1)
        passed = alpha_forward < threshold
        rows.append({
            "q": c[n]["q"],
            "a_next": c[n + 1]["a"],
            "beta_backward": mp.nstr(beta, 30),
            "alpha_forward": mp.nstr(alpha_forward, 30),
            "threshold": mp.nstr(threshold, 30),
            "local_convex": bool(passed),
        })

    expected = {106: False, 33102: False, 33215: True, 66317: True, 265381: False}
    observed = {r["q"]: r["local_convex"] for r in rows}
    for q, status in expected.items():
        assert observed.get(q) == status, (q, observed.get(q), status)

    result = {
        "alpha": "pi",
        "precision_decimal_digits": mp.mp.dps,
        "convergent_denominator_cutoff": QMAX,
        "partial_quotient_one_cases": rows,
        "scope": "Finite high-precision diagnostic; the threshold equivalence is proved symbolically.",
    }
    path = Path(__file__).with_name("partial_quotient_one_threshold_check.json")
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
