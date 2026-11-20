#!/usr/bin/env python3
"""Finite future-veto stopping experiment for selected pi convergents.

Requires mpmath. Decimal outputs are diagnostics, not interval certificates.
"""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 100
ALPHA = mp.pi
QMAX = 10**9


def convergents(alpha, qmax):
    z = mp.mpf(alpha)
    p0, p1 = 0, 1
    q0, q1 = 1, 0
    rows = []
    n = 0
    while True:
        a = int(mp.floor(z))
        p = a * p1 + p0
        q = a * q1 + q0
        if q > qmax:
            break
        E = abs(mp.mpf(q) * alpha - p)
        rows.append({"n": n, "a": a, "p": p, "q": q, "E": E})
        p0, p1 = p1, p
        q0, q1 = q1, q
        frac = z - a
        if frac == 0:
            break
        z = 1 / frac
        n += 1
    return rows


def thresholds(rows, n):
    qn, En = rows[n]["q"], rows[n]["E"]
    past = []
    for i in range(n):
        qi, Ei = rows[i]["q"], rows[i]["E"]
        past.append(((mp.mpf(qn*qn-qi*qi))/(Ei*Ei-En*En), i))
    P, past_i = max(past, default=(mp.mpf(1), None))
    P = max(mp.mpf(1), P)

    U = mp.inf
    arg = None
    trace = []
    for m in range(n + 1, len(rows) - 1):
        qj, Ej = rows[m]["q"], rows[m]["E"]
        x = mp.mpf(qj*qj-qn*qn)/(En*En-Ej*Ej)
        if x < U:
            U, arg = x, m
        qnext = rows[m+1]["q"]
        lower_future = mp.mpf(qnext*qnext-qn*qn)/(En*En)
        trace.append({
            "examined_index": m,
            "examined_q": qj,
            "current_min_x": mp.nstr(U, 35),
            "minimizer_index": arg,
            "minimizer_q": rows[arg]["q"],
            "unseen_future_lower_bound": mp.nstr(lower_future, 35),
            "stop": bool(lower_future >= U),
        })
        if lower_future >= U:
            return {
                "n": n, "q": qn, "p": rows[n]["p"],
                "past_threshold": mp.nstr(P, 35),
                "past_witness_q": rows[past_i]["q"] if past_i is not None else None,
                "future_threshold": mp.nstr(U, 35),
                "future_witness_q": rows[arg]["q"],
                "active": bool(P < U),
                "lookahead_convergents": m - n,
                "stopping_trace": trace,
            }
    raise RuntimeError("Increase QMAX: stopping certificate not reached")


def main():
    rows = convergents(ALPHA, QMAX)
    results = [thresholds(rows, n) for n in (7, 8)]
    assert results[0]["q"] == 99532 and results[0]["active"]
    assert results[1]["q"] == 265381 and not results[1]["active"]
    assert results[1]["future_witness_q"] == 364913

    out = {
        "alpha": "pi",
        "precision_decimal_digits": mp.mp.dps,
        "convergent_denominator_cutoff": QMAX,
        "selected_candidates": results,
        "scope": "High-precision numerical illustration. The finite stopping theorem is exact; decimal diagnostics are not interval certificates.",
    }
    path = Path(__file__).with_name("future_veto_stopping_experiment.json")
    path.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
