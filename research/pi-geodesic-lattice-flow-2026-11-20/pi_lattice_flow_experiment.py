#!/usr/bin/env python3
"""High-precision finite lower-envelope experiment for the pi lattice flow.

Requires mpmath. Numerical output is exploratory, not a proof for unbounded q.
"""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 80
PI = mp.pi
QMAX = 200_000


def nearest_integer(x):
    return int(mp.floor(x + mp.mpf("0.5")))


def cf_convergents(x, count=12):
    y = mp.mpf(x)
    p0, p1 = 0, 1
    q0, q1 = 1, 0
    out = []
    for _ in range(count):
        a = int(mp.floor(y))
        p = a * p1 + p0
        q = a * q1 + q0
        out.append((p, q))
        p0, p1 = p1, p
        q0, q1 = q1, q
        frac = y - a
        if frac == 0:
            break
        y = 1 / frac
    return out


def lower_hull(qmax):
    # Each line is b*x+a, with a=q^2 and b=||q*pi||^2.
    lines = [(mp.mpf(1), 0, 0, 1, mp.mpf(1))]  # vertical vector (0,1)
    for q in range(1, qmax + 1):
        p = nearest_integer(PI * q)
        e = abs(mp.mpf(p) - PI * q)
        lines.append((e * e, q * q, q, p, e))

    # For x increasing, process slopes in descending order.
    lines.sort(key=lambda z: (-z[0], z[1]))
    hull, starts = [], []
    for b, a, q, p, e in lines:
        x_start = mp.ninf
        while hull:
            bp, ap = hull[-1][0], hull[-1][1]
            x_start = (a - ap) / (bp - b)
            if x_start <= starts[-1]:
                hull.pop()
                starts.pop()
            else:
                break
        if not hull:
            x_start = mp.ninf
        hull.append((b, a, q, p, e))
        starts.append(x_start)

    active = []
    for i, line in enumerate(hull):
        lo = starts[i]
        hi = starts[i + 1] if i + 1 < len(starts) else mp.inf
        if hi >= 1:
            q, p, e = line[2], line[3], line[4]
            t_enter = mp.log(lo) / 4 if lo > 0 and lo != mp.ninf else mp.mpf(0)
            t_balance = (mp.log(mp.mpf(q) / e) / 2) if q else mp.mpf(0)
            active.append({
                "q": q, "p": p,
                "error": mp.nstr(e, 30),
                "q_times_error": mp.nstr(q * e, 30) if q else None,
                "x_enter": mp.nstr(lo, 30) if lo != mp.ninf else None,
                "t_enter": mp.nstr(t_enter, 30),
                "t_balance": mp.nstr(t_balance, 30) if q else None,
                "x_exit": mp.nstr(hi, 30) if hi != mp.inf else None,
            })
    return active


def main():
    active = lower_hull(QMAX)
    active_q = {row["q"] for row in active}
    convergents = []
    for p, q in cf_convergents(PI, 12):
        e = abs(mp.mpf(p) - PI * q)
        convergents.append({
            "p": p, "q": q,
            "q_error": mp.nstr(q * e, 30),
            "active_in_qmax_window": q in active_q,
        })

    result = {
        "alpha": "pi",
        "precision_decimal_digits": mp.mp.dps,
        "q_max": QMAX,
        "method": "lower convex hull of lines q^2 + ||q*pi||^2*x, x=exp(4t)",
        "active_vectors": active,
        "continued_fraction_convergents": convergents,
        "scope": "Finite high-precision numerical experiment only; no infinite itinerary claim.",
    }
    out = Path(__file__).with_name("pi_lattice_flow_experiment.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "q_max": QMAX,
        "active_vectors": [(x["q"], x["p"]) for x in active],
        "convergents": [(x["p"], x["q"], x["active_in_qmax_window"]) for x in convergents],
        "output": str(out),
    }, indent=2))


if __name__ == "__main__":
    main()
