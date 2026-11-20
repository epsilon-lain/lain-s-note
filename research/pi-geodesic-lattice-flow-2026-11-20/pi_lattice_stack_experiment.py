#!/usr/bin/env python3
"""Streaming lower-hull/PAVA experiment for the pi lattice flow.

Requires mpmath. The finite-prefix hull is exact as a combinatorial algorithm;
its coordinates are evaluated at high precision. It is not an infinite proof.
"""
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
        error = abs(mp.mpf(q) * alpha - p)
        out.append({"n": n, "a": a, "p": p, "q": q, "E": error})
        p0, p1 = p1, p
        q0, q1 = q1, q
        fraction = z - a
        if fraction == 0:
            break
        z = 1 / fraction
        n += 1
    return out


def drop(a, b):
    """Positive drop slope between points with increasing X and decreasing Y."""
    return (a["Y"] - b["Y"]) / (b["X"] - a["X"])


def streaming_hull(alpha, qmax):
    conv = convergents(alpha, qmax)
    points = [{
        "X": mp.mpf(0), "Y": mp.mpf(1),
        "q": 0, "p": 1, "n": None,
    }]
    points.extend({
        "X": mp.mpf(row["q"]) ** 2,
        "Y": row["E"] ** 2,
        "q": row["q"], "p": row["p"], "n": row["n"],
    } for row in conv)

    stack = []
    popped = []
    for point in points:
        while len(stack) >= 2 and drop(stack[-2], stack[-1]) <= drop(stack[-1], point):
            popped.append(stack.pop())
        stack.append(point)

    active = []
    for i, point in enumerate(stack):
        left_drop = mp.mpf(0) if i == 0 else drop(stack[i - 1], point)
        right_drop = mp.mpf(0) if i == len(stack) - 1 else drop(point, stack[i + 1])
        x_start = max(mp.mpf(1), 1 / left_drop) if left_drop > 0 else mp.mpf(1)
        x_end = mp.inf if right_drop == 0 else 1 / right_drop
        if x_start < x_end:
            active.append({
                "q": point["q"], "p": point["p"], "cf_index": point["n"],
                "x_start": mp.nstr(x_start, 40),
                "x_end": None if x_end == mp.inf else mp.nstr(x_end, 40),
                "t_start": mp.nstr(mp.log(x_start) / 4, 35),
                "t_end": None if x_end == mp.inf else mp.nstr(mp.log(x_end) / 4, 35),
            })

    return conv, stack, popped, active


def main():
    conv, stack, popped, active = streaming_hull(ALPHA, QMAX)
    hull_pairs = [(p["q"], p["p"]) for p in stack]
    popped_pairs = [(p["q"], p["p"]) for p in popped]

    expected = [
        (0, 1), (1, 3), (7, 22), (113, 355),
        (33215, 104348), (66317, 208341),
        (99532, 312689), (364913, 1146408),
    ]
    assert hull_pairs == expected
    assert (106, 333) in popped_pairs
    assert (33102, 103993) in popped_pairs
    assert (265381, 833719) in popped_pairs

    result = {
        "alpha": "pi",
        "precision_decimal_digits": mp.mp.dps,
        "convergent_denominator_cutoff": QMAX,
        "convergents_processed": len(conv),
        "hull_vertices_including_vertical": hull_pairs,
        "popped_convergents": popped_pairs,
        "active_intervals": active,
        "algorithm": "append convergent; while previous_drop <= new_drop, pop middle; append",
        "complexity": {
            "time": "O(number of convergents)",
            "space": "O(number of hull vertices)",
            "each_point_pushed_once": True,
            "each_point_popped_at_most_once": True,
        },
        "scope": "Finite-prefix numerical experiment. The stack theorem is exact; the displayed pi data are cutoff-dependent.",
    }
    out = Path(__file__).with_name("pi_lattice_stack_experiment.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "convergents_processed": len(conv),
        "hull_vertices": hull_pairs,
        "popped": popped_pairs,
        "output": str(out),
    }, indent=2))


if __name__ == "__main__":
    main()
