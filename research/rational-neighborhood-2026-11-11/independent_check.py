#!/usr/bin/env python3
"""Exact, independent finite checks for the rational-neighborhood note.

No floating-point arithmetic and no imports from the note's construction code.
Continued fractions are extracted by exact arithmetic in Q(sqrt(2)); that
extraction is independent of the proposed parent/tail classification.
Finite checks supplement, and do not replace, the general proofs in the note.
"""

from dataclasses import dataclass
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import hashlib
import json


@dataclass(frozen=True)
class Quad:
    """a + b*sqrt(2), with exact rational coefficients."""

    a: F
    b: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "a", F(self.a))
        object.__setattr__(self, "b", F(self.b))

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Quad) else Quad(F(value))

    def __add__(self, other):
        other = self.coerce(other)
        return Quad(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Quad(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        return Quad(self.a * other.a + 2 * self.b * other.b,
                    self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def reciprocal(self):
        norm = self.a * self.a - 2 * self.b * self.b
        if not norm:
            raise ZeroDivisionError
        return Quad(self.a / norm, -self.b / norm)

    def __truediv__(self, other):
        return self * self.coerce(other).reciprocal()

    def __rtruediv__(self, other):
        return self.coerce(other) / self

    def sign(self):
        if not self.b:
            return (self.a > 0) - (self.a < 0)
        if not self.a:
            return (self.b > 0) - (self.b < 0)
        if self.a > 0 and self.b > 0:
            return 1
        if self.a < 0 and self.b < 0:
            return -1
        norm = self.a * self.a - 2 * self.b * self.b
        # A nonzero rational cannot equal a nonzero rational times sqrt(2).
        assert norm
        return ((norm > 0) - (norm < 0)) * ((self.a > 0) - (self.a < 0))

    def __abs__(self):
        return self if self.sign() >= 0 else -self

    def __lt__(self, other):
        return (self - other).sign() < 0

    def floor(self):
        # 1 < sqrt(2) < 2 gives an exact initial integer bracket.
        ends = [self.a + self.b, self.a + 2 * self.b]
        lo = min(ends).numerator // min(ends).denominator
        hi = max(ends).numerator // max(ends).denominator + 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if self < mid:
                hi = mid
            else:
                lo = mid
        assert not self < lo and self < lo + 1
        return lo

    def record(self):
        return {"rational": str(self.a), "sqrt2_coefficient": str(self.b)}


def convergents_through_denominator(x, bound):
    """Extract by the usual floor/reciprocal algorithm, stopping above bound."""
    n0, n1, d0, d1 = 0, 1, 1, 0
    entries = []
    for _ in range(128):
        digit = x.floor()
        numerator = digit * n1 + n0
        denominator = digit * d1 + d0
        entries.append((numerator, denominator, digit))
        if denominator > bound:
            return entries
        x = (x - digit).reciprocal()
        n0, n1, d0, d1 = n1, numerator, d1, denominator
    raise AssertionError("Unexpected continued-fraction iteration limit")


def parents(p, q):
    s = pow(p, -1, q)
    r = (p * s - 1) // q
    return ((r, s, -1), (p - r, q - s, 1))


def run():
    centers = [(p + shift * q, q)
               for q in range(2, 51)
               for p in range(1, q)
               if gcd(p, q) == 1
               for shift in (-2, 0, 3)]
    centers.extend([(22, 7), (355, 113), (-355, 113)])
    centers = sorted(set(centers), key=lambda item: (item[1], item[0]))
    counts = {"centers": len(centers), "parent_determinants": 0,
              "quadratic_chart_cases": 0, "inside_ball_cases": 0,
              "outside_ball_cases": 0, "convergent_classifications": 0,
              "rational_tail_transfers": 0, "cf_entries": 0}
    examples = []
    rational_tails = [F(1, 3), F(1), F(7, 5), F(17, 3)]
    for p, q in centers:
        for r, s, side in parents(p, q):
            assert 0 < s < q
            assert q * r - p * s == side
            assert gcd(r, s) == 1
            counts["parent_determinants"] += 1
            threshold = F(q - s, q)
            tail_cases = [
                ("outside", Quad(threshold, F(-1, 4 * q))),
                ("near_outer_boundary", Quad(threshold, F(1, 4 * q))),
                ("below_cf_threshold", Quad(1, F(-1, 4 * q))),
                ("above_cf_threshold", Quad(1, F(1, 4 * q))),
                ("moderate_tail", Quad(1, 1)),
                ("large_tail", Quad(37, 1)),
            ]
            for label, t in tail_cases:
                assert t.sign() > 0
                x = (p * t + r) / (q * t + s)
                assert x.b and (x - F(p, q)).sign() == side
                assert (r - s * x) / (q * x - p) == t
                quality = q * q * abs(x - F(p, q))
                assert quality == 1 / (t + F(s, q))
                inside = quality < 1
                assert inside == (Quad(threshold) < t)
                counts["inside_ball_cases" if inside else "outside_ball_cases"] += 1
                cf = convergents_through_denominator(x, q)
                found = any((n, d) == (p, q) for n, d, _ in cf)
                assert found == (Quad(1) < t), (p, q, r, s, label, cf)
                counts["quadratic_chart_cases"] += 1
                counts["convergent_classifications"] += 1
                counts["cf_entries"] += len(cf)
                for tail in rational_tails:
                    a, b = tail.numerator, tail.denominator
                    P, Q = p * a + r * b, q * a + s * b
                    assert Q > 0 and gcd(P, Q) == 1
                    lhs = Q * Q * abs(x - F(P, Q))
                    rhs = b * b * abs(t - tail) * (q * tail + s) / (q * t + s)
                    assert lhs == rhs
                    counts["rational_tail_transfers"] += 1
                if (p, q) == (22, 7):
                    examples.append({
                        "center": "22/7", "parent": str(F(r, s)),
                        "side": side, "label": label,
                        "tail": t.record(), "x": x.record(),
                        "inside_ball": inside, "center_is_convergent": found,
                        "convergents_checked": [f"{n}/{d}" for n, d, _ in cf],
                    })
    return {
        "status": "PASS",
        "scope": "Finite exact checks; general statements require the accompanying proofs.",
        "arithmetic": "fractions.Fraction and exact Q(sqrt(2)) arithmetic; no floating point",
        "independence": "CFs are computed by floor/reciprocal, without using the parent/tail criterion.",
        "ranges": {"q_range": [2, 50], "numerators": "coprime residues plus shifts -2q, 0, 3q",
                   "additional_centers": ["22/7", "355/113", "-355/113"],
                   "tail_cases_per_side": 6,
                   "rational_tail_probes": [str(t) for t in rational_tails]},
        "counts": counts,
        "examples": examples,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


if __name__ == "__main__":
    result = run()
    output = Path(__file__).with_name("verification.json")
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "counts": result["counts"]}, sort_keys=True))
