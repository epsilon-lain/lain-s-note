#!/usr/bin/env python3
"""Exact rational angular quadrature, cancellation, and denominator accounting.

Candidate generation uses no approximation to pi.  A Machin interval is used only
after the finite candidate grid has been generated.  Standard library only.
"""

from decimal import Decimal, localcontext
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path
import json


CS = [F(-1, 2), F(-1, 4), F(0), F(1, 4), F(1, 2)]
TRAP_NS = [1, 2, 4, 8, 16]
SIMPSON_NS = [2, 4, 8, 16]
AITKEN_NS = [1, 2, 4]


def kernel(t, c):
    """4 z'(t)/(1+z(t)^2), z=t+c*t*(1-t), endpoint z=0,1."""
    z = t + c * t * (1 - t)
    dz = 1 + c - 2 * c * t
    return 4 * dz / (1 + z * z)


def reduce_pair(p, q):
    assert q != 0
    if q < 0:
        p, q = -p, -q
    g = gcd(abs(p), q)
    r = F(p, q)
    assert r.numerator == p // g and r.denominator == q // g
    return r, {
        "raw_numerator": str(p),
        "raw_denominator": str(q),
        "actual_gcd": str(g),
        "raw_numerator_bits": abs(p).bit_length(),
        "raw_denominator_bits": q.bit_length(),
    }


def trapezoid(n, c):
    values = [kernel(F(j, n), c) for j in range(n + 1)]
    d = lcm(*(v.denominator for v in values))
    p = sum((1 if j in (0, n) else 2) * v.numerator * (d // v.denominator)
            for j, v in enumerate(values))
    r, cost = reduce_pair(p, 2 * n * d)
    cost["kernel_denominator_max_bits"] = max(v.denominator.bit_length() for v in values)
    cost["kernel_common_denominator_bits"] = d.bit_length()
    return r, cost


def simpson(n, c):
    assert n % 2 == 0
    values = [kernel(F(j, n), c) for j in range(n + 1)]
    d = lcm(*(v.denominator for v in values))
    p = sum((1 if j in (0, n) else (4 if j % 2 else 2))
            * v.numerator * (d // v.denominator)
            for j, v in enumerate(values))
    r, cost = reduce_pair(p, 3 * n * d)
    # This equivalence is an exact algebraic check, rather than a floating test.
    tn, _ = trapezoid(n, c)
    th, _ = trapezoid(n // 2, c)
    assert r == (4 * tn - th) / 3
    cost["kernel_denominator_max_bits"] = max(v.denominator.bit_length() for v in values)
    cost["kernel_common_denominator_bits"] = d.bit_length()
    return r, cost


def aitken(n, c):
    a, b, e = [trapezoid(k * n, c)[0] for k in (1, 2, 4)]
    d = lcm(a.denominator, b.denominator, e.denominator)
    aa, bb, ee = [r.numerator * (d // r.denominator) for r in (a, b, e)]
    delta = aa - 2 * bb + ee
    assert delta != 0
    r, cost = reduce_pair(aa * ee - bb * bb, d * delta)
    assert r == (a * e - b * b) / (a - 2 * b + e)
    cost["source_fractions"] = [str(a), str(b), str(e)]
    cost["second_difference"] = str(a - 2 * b + e)
    values = [kernel(F(j, 4 * n), c) for j in range(4 * n + 1)]
    cost["kernel_denominator_max_bits"] = max(v.denominator.bit_length() for v in values)
    cost["kernel_common_denominator_bits"] = lcm(*(v.denominator for v in values)).bit_length()
    return r, cost


def atan_interval(k, terms):
    # Alternating series with the next term gives strict enclosing endpoints.
    s = sum((F((-1) ** j, (2 * j + 1) * k ** (2 * j + 1))
             for j in range(terms)), F())
    e = s + F((-1) ** terms, (2 * terms + 1) * k ** (2 * terms + 1))
    return min(s, e), max(s, e)


def pi_interval():
    # tan(4 atan(1/5))=120/119 and tan(theta-atan(1/239))=1.
    u = 2 * F(1, 5) / (1 - F(1, 5) ** 2)
    v = 2 * u / (1 - u * u)
    assert v == F(120, 119)
    assert (v - F(1, 239)) / (1 + v * F(1, 239)) == 1
    a, b = atan_interval(5, 70), atan_interval(239, 20)
    return 16 * a[0] - 4 * b[1], 16 * a[1] - 4 * b[0]


def decimal_display(x, digits=16):
    with localcontext() as ctx:
        ctx.prec = digits
        return str(Decimal(x.numerator) / Decimal(x.denominator))


def absolute_error_interval(r, pi_bounds):
    lo, hi = pi_bounds
    lower = 0 if lo <= r <= hi else min(abs(lo - r), abs(hi - r))
    return F(lower), max(abs(lo - r), abs(hi - r))


def score(r, pi_bounds):
    e0, e1 = absolute_error_interval(r, pi_bounds)
    q = r.denominator
    return {
        "absolute_error_interval": [str(e0), str(e1)],
        "quality_interval": [str(q * q * e0), str(q * q * e1)],
        "absolute_error_display": decimal_display((e0 + e1) / 2),
        "quality_display": decimal_display(q * q * (e0 + e1) / 2),
    }


def common_multiplier_checks():
    cases = 0
    # C=1+5t-3t², C>=1 on [0,1].  Compare the scalar-multiplied
    # spinor formula directly, including both derivative product terms.
    for c in CS:
        for n in TRAP_NS:
            for j in range(n + 1):
                t = F(j, n)
                p, pp = F(1), F(0)
                q, qp = t + c * t * (1 - t), 1 + c - 2 * c * t
                cc, cp = 1 + 5 * t - 3 * t * t, 5 - 6 * t
                assert cc > 0
                p2, q2, p2p, q2p = cc * p, cc * q, cp * p + cc * pp, cp * q + cc * qp
                direct = 4 * (p2 * q2p - q2 * p2p) / (p2 * p2 + q2 * q2)
                assert direct == kernel(t, c)
                cases += 1
    return cases


def error_identity_checks():
    cases = 0
    for alpha in [F(-7, 3), F(0), F(11, 4)]:
        for e0 in [F(-1, 3), F(2, 7)]:
            for e1 in [F(-1, 5), F(3, 11)]:
                for e2 in [F(-1, 9), F(4, 13)]:
                    a, b, e = alpha + e0, alpha + e1, alpha + e2
                    if a - 2 * b + e == 0:
                        continue
                    lhs = (a * e - b * b) / (a - 2 * b + e) - alpha
                    rhs = (e0 * e2 - e1 * e1) / (e0 - 2 * e1 + e2)
                    assert lhs == rhs
                    cases += 1
    # Exactly geometric error is eliminated without knowledge of alpha.
    for alpha in [F(7, 5), F(-5, 11)]:
        for amplitude in [F(1, 3), F(-7, 2)]:
            for ratio in [F(1, 4), F(-1, 2), F(2)]:
                a, b, e = [alpha + amplitude * ratio ** k for k in range(3)]
                assert (a * e - b * b) / (a - 2 * b + e) == alpha
                cases += 1
    return cases


def main():
    candidates = []
    # Generate the complete finite grid before accessing any pi interval.
    for c in CS:
        for method, ns, make in [
            ("trapezoid", TRAP_NS, trapezoid),
            ("simpson", SIMPSON_NS, simpson),
            ("aitken", AITKEN_NS, aitken),
        ]:
            for n in ns:
                r, cost = make(n, c)
                max_grid = 4 * n if method == "aitken" else n
                candidates.append({
                    "method": method, "base_n": n, "c": str(c),
                    "parameter_height": max(abs(c.numerator), c.denominator),
                    "parameter_count": 1,
                    "integer_spinor": {"P": [c.denominator],
                                        "Q": [0, c.numerator + c.denominator, -c.numerator]},
                    "integer_spinor_height": max(c.denominator, abs(c.numerator + c.denominator), abs(c.numerator)),
                    "unique_rational_nodes": max_grid + 1,
                    "max_node_denominator": max_grid,
                    "fraction": str(r), "p": str(r.numerator), "q": str(r.denominator),
                    "p_bits": abs(r.numerator).bit_length(),
                    "q_bits": r.denominator.bit_length(), "cost": cost,
                })
    pi_bounds = pi_interval()
    assert pi_bounds[0] < pi_bounds[1] and pi_bounds[1] - pi_bounds[0] < F(1, 10 ** 95)
    for row in candidates:
        row.update(score(F(row["fraction"]), pi_bounds))
    rows = {(r["method"], r["base_n"], r["c"]): r for r in candidates}
    comparisons = []
    for n in [4, 8, 16]:
        base = rows[("trapezoid", n, "0")]
        new = rows[("simpson", n, "0")]
        old0, old1 = map(F, base["quality_interval"])
        new0, new1 = map(F, new["quality_interval"])
        assert new1 < old0
        comparisons.append({
            "baseline": ["trapezoid", n, "0"],
            "candidate": ["simpson", n, "0"],
            "quality_improvement_factor_interval": [str(old0 / new1), str(old1 / new0)],
            "quality_improvement_factor_display": decimal_display(old0 / new1),
            "same_nodes": True,
        })
    for n in [1, 2, 4]:
        base = rows[("trapezoid", 4 * n, "0")]
        new = rows[("aitken", n, "0")]
        old0, old1 = map(F, base["quality_interval"])
        new0, new1 = map(F, new["quality_interval"])
        assert (new1 < old0) if n < 4 else (new0 > old1)
        comparisons.append({
            "baseline": ["trapezoid", 4 * n, "0"],
            "candidate": ["aitken", n, "0"],
            "quality_improvement_factor_interval": [str(old0 / new1), str(old1 / new0)],
            "quality_improvement_factor_display": decimal_display(old0 / new1),
            "same_nodes": True,
        })
    output = {
        "schema": 1,
        "scope": "Finite exact arithmetic experiment; no new pi record, no asymptotic claim.",
        "candidate_generation_uses_pi": False,
        "classical_methods": ["composite trapezoid", "composite Simpson / Richardson", "Aitken delta-squared"],
        "pi_certificate": {
            "identity": "pi = 16 atan(1/5) - 4 atan(1/239)",
            "alternating_terms": {"atan(1/5)": 70, "atan(1/239)": 20},
            "interval": [str(x) for x in pi_bounds],
            "width": str(pi_bounds[1] - pi_bounds[0]),
        },
        "candidate_count": len(candidates),
        "posthoc_reference_baselines": {
            "22/7": score(F(22, 7), pi_bounds),
            "355/113": score(F(355, 113), pi_bounds),
        },
        "candidates": candidates, "comparisons": comparisons,
        "checks": {
            "all_candidates_fully_reduced": all(gcd(int(r["p"]), int(r["q"])) == 1 for r in candidates),
            "common_multiplier_cases": common_multiplier_checks(),
            "aitken_error_identity_cases": error_identity_checks(),
            "status": "PASS",
        },
    }
    target = Path(__file__).with_name("results.json")
    target.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "candidates": len(candidates),
                      "common_multiplier_cases": output["checks"]["common_multiplier_cases"],
                      "aitken_error_identity_cases": output["checks"]["aitken_error_identity_cases"],
                      "output": str(target)}, indent=2))


if __name__ == "__main__":
    main()
