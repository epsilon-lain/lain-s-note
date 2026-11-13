#!/usr/bin/env python3
"""Exact checks for the critical-scale denominator jump example.

Run from any directory:
    python3 check_critical_scale_jump.py

Uses only the Python standard library. It reads the independently certified
runB finite-witness record already stored next to this research package.
It does not recompute pi or certify an infinite statement.
"""
from decimal import Decimal, localcontext
from fractions import Fraction
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = (
    HERE.parent
    / "yohaku-pi-cf-audit-2026-11-12"
    / "runB-independent.json"
)

data = json.loads(SOURCE.read_text(encoding="utf-8"))
assert data["certified_cf_terms"] == 28422
assert data["a28421"] == 78629
assert data["p28420_is_nearest_integer_to_q28420_pi"] is True

a = int(data["a28421"])
lower = Fraction(data["strict_witness_lower_40"])
upper = Fraction(data["strict_witness_upper_40"])
assert 0 < lower < upper

# The convergent identity gives 1/(a+2) < eps_n < 1/a.
coarse_lower = Fraction(1, a + 2)
coarse_upper = Fraction(1, a)
assert coarse_lower < lower < upper < coarse_upper

# 1/eps = R + theta with 0 < theta < 1, hence
# a < R < 1/eps < 1/lower.
reciprocal_lower = 1 / upper
reciprocal_upper = 1 / lower
ratio_lower = Fraction(a, 1)
ratio_upper = reciprocal_upper
assert ratio_lower < ratio_upper
assert ratio_upper < Fraction(a, 1) + Fraction(598, 1000)

with localcontext() as ctx:
    ctx.prec = 60
    ratio_upper_decimal = str(
        Decimal(ratio_upper.numerator) / Decimal(ratio_upper.denominator)
    )
    reciprocal_lower_decimal = str(
        Decimal(reciprocal_lower.numerator) / Decimal(reciprocal_lower.denominator)
    )
    reciprocal_upper_decimal = str(
        Decimal(reciprocal_upper.numerator) / Decimal(reciprocal_upper.denominator)
    )

result = {
    "format": "critical-scale-denominator-jump/1",
    "source": "yohaku-pi-cf-audit-2026-11-12/runB-independent.json",
    "convergent_index": 28420,
    "next_partial_quotient": a,
    "epsilon_interval_open": [
        data["strict_witness_lower_40"],
        data["strict_witness_upper_40"],
    ],
    "coarse_cf_epsilon_interval": [
        f"1/{a + 2}",
        f"1/{a}",
    ],
    "reciprocal_epsilon_interval_open_decimal": [
        reciprocal_lower_decimal,
        reciprocal_upper_decimal,
    ],
    "denominator_ratio_interval_open": [
        str(ratio_lower),
        f"{a + 598}/1000",
    ],
    "denominator_ratio_upper_exact": (
        f"{ratio_upper.numerator}/{ratio_upper.denominator}"
    ),
    "denominator_ratio_upper_decimal": ratio_upper_decimal,
    "checks": {
        "epsilon_inside_cf_interval": True,
        "reciprocal_upper_below_78629_598": True,
        "ratio_lower_from_next_partial_quotient": True,
        "ratio_upper_from_1_over_epsilon": True,
    },
    "scope": (
        "Exact rational checks for one finite witness; no new pi computation, "
        "no tail liminf result, and no proof of P2848."
    ),
}

out = HERE / "results.json"
out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
