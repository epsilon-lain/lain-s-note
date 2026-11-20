#!/usr/bin/env python3
"""Compare the full denominator hull with the convergent-only hull for pi."""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 100
PI = mp.pi


def convergents(alpha, qmax):
    z = mp.mpf(alpha)
    p0, p1 = 0, 1
    q0, q1 = 1, 0
    out = []
    while True:
        a = int(mp.floor(z))
        p = a * p1 + p0
        q = a * q1 + q0
        if q > qmax:
            break
        e = abs(mp.mpf(p) - alpha * q)
        out.append((q, p, e))
        p0, p1 = p1, p
        q0, q1 = q1, q
        frac = z - a
        if frac == 0:
            break
        z = 1 / frac
    return out


def hull(lines):
    # line data: (slope, intercept, q, p, error)
    lines = sorted(lines, key=lambda z: (-z[0], z[1]))
    hs, starts = [], []
    for line in lines:
        b, a = line[0], line[1]
        x0 = mp.ninf
        while hs:
            x0 = (a - hs[-1][1]) / (hs[-1][0] - b)
            if x0 <= starts[-1]:
                hs.pop()
                starts.pop()
            else:
                break
        if not hs:
            x0 = mp.ninf
        hs.append(line)
        starts.append(x0)
    active = []
    for i, line in enumerate(hs):
        right = starts[i + 1] if i + 1 < len(starts) else mp.inf
        if right >= 1:
            active.append({
                "q": line[2], "p": line[3],
                "error": mp.nstr(line[4], 35),
                "q_error": mp.nstr(line[2] * line[4], 35),
                "x_enter": mp.nstr(starts[i], 35),
                "t_enter": mp.nstr(max(mp.mpf(0), mp.log(starts[i]) / 4), 35)
                if starts[i] > 0 and starts[i] != mp.ninf else "0",
            })
    return active


def convergent_hull(alpha, qmax):
    conv = convergents(alpha, qmax)
    lines = [(e * e, q * q, q, p, e) for q, p, e in conv]
    active = hull(lines)
    return conv, active


def full_hull(alpha, qmax):
    lines = [(mp.mpf(1), 0, 0, 1, mp.mpf(1))]
    for q in range(1, qmax + 1):
        p = int(mp.floor(alpha * q + mp.mpf("0.5")))
        e = abs(mp.mpf(p) - alpha * q)
        lines.append((e * e, q * q, q, p, e))
    return hull(lines)


def main():
    qmax = 200_000
    conv, cf_active = convergent_hull(PI, qmax)
    full_active = full_hull(PI, qmax)
    cf_pairs = [(x["q"], x["p"]) for x in cf_active]
    full_pairs = [(x["q"], x["p"]) for x in full_active]
    assert cf_pairs == full_pairs, (cf_pairs, full_pairs)

    result = {
        "q_max": qmax,
        "precision": mp.mp.dps,
        "convergents_count": len(conv),
        "full_scan_active": full_pairs,
        "convergent_only_active": cf_pairs,
        "identical": True,
        "scope": "Finite numerical regression; infinite subset theorem is proved separately.",
    }
    path = Path(__file__).with_name("continued_fraction_hull_check.json")
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
