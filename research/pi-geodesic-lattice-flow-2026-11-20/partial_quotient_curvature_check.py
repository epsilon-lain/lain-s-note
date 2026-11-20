#!/usr/bin/env python3
"""Check local partial-quotient curvature against finite lower-hull activity."""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 100
ALPHA = mp.pi
QMAX = 200_000


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
        e = abs(mp.mpf(q) * alpha - p)
        out.append({"n": n, "a": a, "p": p, "q": q, "E": e})
        p0, p1 = p1, p
        q0, q1 = q1, q
        frac = z - a
        if frac == 0:
            break
        z = 1 / frac
        n += 1
    return out


def local_test(c, n):
    prev, cur, nxt = c[n-1], c[n], c[n+1]
    a = nxt["a"]  # a_{n+1} in q_{n+1}=a_{n+1}q_n+q_{n-1}
    qprev, q, qnext = prev["q"], cur["q"], nxt["q"]
    Eprev, E, Enext = prev["E"], cur["E"], nxt["E"]
    Dleft = (Eprev**2 - E**2) / (q**2 - qprev**2)
    Dright = (E**2 - Enext**2) / (qnext**2 - q**2)

    r = mp.mpf(qprev) / q
    z = Enext / E
    lhs = ((a + z)**2 - 1) / (1 - r**2)
    rhs = (1 - z**2) / ((a + r)**2 - 1)
    assert abs(Dleft / (E/q)**2 - lhs) < mp.mpf('1e-70')
    assert abs(Dright / (E/q)**2 - rhs) < mp.mpf('1e-70')
    return {
        "n": n, "p": cur["p"], "q": q, "a_n": a,
        "D_left": mp.nstr(Dleft, 35),
        "D_right": mp.nstr(Dright, 35),
        "local_convex": bool(Dleft > Dright),
        "partial_quotient_test": bool(lhs > rhs),
    }


def finite_hull(c):
    lines = [(mp.mpf(1), 0, 0, 1)]  # vertical vector
    for row in c:
        q, E = row["q"], row["E"]
        lines.append((E**2, q**2, q, row["p"]))
    lines.sort(key=lambda z: (-z[0], z[1]))
    h, starts = [], []
    for line in lines:
        b, a = line[0], line[1]
        x0 = mp.ninf
        while h:
            x0 = (a - h[-1][1]) / (h[-1][0] - b)
            if x0 <= starts[-1]:
                h.pop()
                starts.pop()
            else:
                break
        if not h:
            x0 = mp.ninf
        h.append(line)
        starts.append(x0)
    active = set()
    for i, line in enumerate(h):
        right = starts[i+1] if i+1 < len(starts) else mp.inf
        if right >= 1:
            active.add(line[2])
    return active


def main():
    c = convergents(ALPHA, QMAX)
    active = finite_hull(c)
    tests = []
    for n in range(1, len(c)-1):
        row = local_test(c, n)
        row["active_in_finite_hull"] = row["q"] in active
        tests.append(row)

    # The two well-known skipped convergents must fail the local test.
    skipped = {106, 33102}
    assert all(not x["local_convex"] for x in tests if x["q"] in skipped)
    # Local convexity is necessary but not asserted sufficient globally.
    assert all(not x["active_in_finite_hull"] or x["local_convex"] for x in tests)

    result = {
        "alpha": "pi",
        "q_max": QMAX,
        "precision_decimal_digits": mp.mp.dps,
        "convergents_count": len(c),
        "local_tests": tests,
        "scope": "Finite numerical audit of the derived local formula; global activity is determined by the global interval criterion.",
    }
    out = Path(__file__).with_name("partial_quotient_curvature_check.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("convergents:", len(c))
    print("skipped 106 and 33102 rejected:", True)
    print("active implies local convexity in finite hull:", True)
    print("saved:", out)


if __name__ == "__main__":
    main()
