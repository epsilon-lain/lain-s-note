#!/usr/bin/env python3
"""Independent verification of the stored rational-neighborhood examples.

Parents are rebuilt with the Euclidean continued-fraction algorithm, without
modular inverses. Pi uses a dyadic enclosure for
    pi = 4 * (arctan(1/2) + arctan(1/3)),
instead of the generating script's Machin formula. All arithmetic is exact.
The generator is neither imported nor executed.
"""

from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

from independent_check import Quad, convergents_through_denominator


def finite_cf(value):
    digits = []
    while True:
        digit = value.numerator // value.denominator
        digits.append(digit)
        remainder = value - digit
        if not remainder:
            return digits
        value = 1 / remainder


def rational_convergents(digits):
    n0, n1, d0, d1 = 0, 1, 1, 0
    values = []
    for digit in digits:
        n0, n1 = n1, digit * n1 + n0
        d0, d1 = d1, digit * d1 + d0
        values.append(F(n1, d1))
    return values


def euclidean_parents(center):
    short = finite_cf(center)
    assert len(short) >= 2 and short[-1] >= 2
    long = short[:-1] + [short[-1] - 1, 1]
    cv_short, cv_long = rational_convergents(short), rational_convergents(long)
    assert cv_short[-1] == cv_long[-1] == center
    return sorted((cv_short[-2], cv_long[-2]))


def dyadic_arctan(inv, bits=512, terms=300):
    """Outward integer rounding termwise, plus the alternating tail bound."""
    scale = 1 << bits
    lo, hi = 0, 0
    for k in range(terms):
        denominator = (2 * k + 1) * inv ** (2 * k + 1)
        numerator = scale if k % 2 == 0 else -scale
        lo += numerator // denominator
        hi += -((-numerator) // denominator)
    denominator = (2 * terms + 1) * inv ** (2 * terms + 1)
    numerator = scale if terms % 2 == 0 else -scale
    lo += min(0, numerator // denominator)
    hi += max(0, -((-numerator) // denominator))
    return F(lo, scale), F(hi, scale)


def independent_pi():
    lo2, hi2 = dyadic_arctan(2)
    lo3, hi3 = dyadic_arctan(3)
    return 4 * (lo2 + lo3), 4 * (hi2 + hi3)


def check_decimal(printed, interval):
    low, high = map(F, printed)
    lo, hi = interval
    assert low <= lo < hi <= high
    assert high - low == F(1, 10 ** 24)


def run():
    folder = Path(__file__).resolve().parent
    path = folder / "examples.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["status"] == "PASS"
    counts = {"stored_centers": 0, "euclidean_parent_pairs": 0,
              "stored_integer_centers": 0, "integer_cf_probes": 0,
              "decimal_enclosures": 0, "tail_transfer_endpoints": 0}
    rebuilt = []
    for center_data in data["centers"]:
        center = F(center_data["center"])
        p, q = center.numerator, center.denominator
        assert (p, q) == (center_data["p"], center_data["q"])
        assert list(map(F, center_data["ball"])) == [center-F(1, q*q), center+F(1, q*q)]
        counts["stored_centers"] += 1
        if q == 1:
            assert list(map(F, center_data["integer_case"]["convergent_interval"])) == [center-F(1, 2), center+1]
            counts["stored_integer_centers"] += 1
            continue
        parent_values = euclidean_parents(center)
        endpoints = []
        radii = []
        for branch, parent, side in zip(center_data["branches"], parent_values, (-1, 1)):
            r, s = parent.numerator, parent.denominator
            assert (r, s) == (branch["r"], branch["s"])
            assert F(branch["parent"]) == parent
            assert p*s-q*r == branch["determinant_p_s_minus_q_r"] == -side
            assert F(branch["tail_lower_bound_strict"]) == F(q-s, q)
            endpoint = F(p+r, q+s)
            radius = abs(endpoint-center)
            assert F(branch["convergent_outer_endpoint"]) == endpoint
            assert F(branch["convergent_radius"]) == radius
            endpoints.append(endpoint)
            radii.append(radius)
        assert list(map(F, center_data["convergent_interval"])) == endpoints
        assert F(center_data["sharp_symmetric_convergent_radius"]) == min(radii)
        assert F(center_data["convergent_length_fraction"]) == (endpoints[1]-endpoints[0])*q*q/2
        rebuilt.append({"center": str(center), "euclidean_parents": list(map(str, parent_values)),
                        "cylinder_endpoints": list(map(str, endpoints))})
        counts["euclidean_parent_pairs"] += 1

    offsets = [Quad(-1, F(1, 100)), Quad(F(-3, 4), F(1, 100)),
               Quad(F(-1, 2), F(-1, 100)), Quad(F(-1, 2), F(1, 100)),
               Quad(0, F(-1, 4)), Quad(0, F(1, 4)), Quad(1, F(-1, 100))]
    for p in range(-5, 6):
        for offset in offsets:
            x = p + offset
            assert Quad(p-1) < x < p+1
            cv = convergents_through_denominator(x, 1)
            found = any((n, d) == (p, 1) for n, d, _ in cv)
            assert found == (Quad(F(2*p-1, 2)) < x)
            counts["integer_cf_probes"] += 1

    pi_data = data["pi_example"]
    independent_lo, independent_hi = independent_pi()
    root_lo, root_hi = map(F, pi_data["pi_interval_exact"])
    assert root_lo < independent_lo < independent_hi < root_hi
    assert F(3) < independent_lo < independent_hi < F(22, 7)
    assert independent_hi < F(355, 113)

    def inverse(x):
        return (3-x)/(7*x-22)

    root_tlo, root_thi = map(F, pi_data["lower_22_over_7_tail_exact"])
    assert root_tlo == inverse(root_lo) and root_thi == inverse(root_hi)
    tlo, thi = inverse(independent_lo), inverse(independent_hi)
    assert root_tlo < tlo < thi < root_thi < 16
    assert F(pi_data["tail_16_image"]) == F(22*16+3, 7*16+1) == F(355, 113)

    exact_intervals = {
        "pi_interval_decimal": (root_lo, root_hi),
        "tail_decimal": (root_tlo, root_thi),
        "tail_minus_16_abs_decimal": (16-root_thi, 16-root_tlo),
        "quality_22_over_7_decimal": (49*(F(22, 7)-root_hi), 49*(F(22, 7)-root_lo)),
        "quality_355_over_113_decimal": (113**2*(F(355, 113)-root_hi), 113**2*(F(355, 113)-root_lo)),
    }
    for key, interval in exact_intervals.items():
        check_decimal(pi_data[key], interval)
        counts["decimal_enclosures"] += 1
    for pi_bound, tail_bound in ((independent_lo, tlo), (independent_hi, thi)):
        transformed = (16-tail_bound)*113/(7*tail_bound+1)
        assert transformed == 113**2*(F(355, 113)-pi_bound)
        counts["tail_transfer_endpoints"] += 1

    return {
        "status": "PASS", "counts": counts,
        "arithmetic": "exact integers, Fraction, Q(sqrt(2)); no floating point",
        "scope": "Stored examples only, plus 77 exact integer-center CF probes; not a formal proof or new pi result.",
        "independence": {"parents": "Euclidean CF predecessors from both finite expansions; no modular inverse",
                         "pi": "4*(arctan(1/2)+arctan(1/3)), each 300 terms with outward rounding to 512-bit dyadics",
                         "pi_containment": "The independent interval is strictly inside the generator's exact Machin interval",
                         "generator": "neighborhood.py was neither imported nor executed",
                         "integer_cases": "Uses only the independent_check.py quadratic arithmetic and direct CF extraction"},
        "independent_pi_interval": [str(independent_lo), str(independent_hi)],
        "rebuilt_centers": rebuilt,
        "examples_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


if __name__ == "__main__":
    result = run()
    Path(__file__).with_name("examples_verification.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "counts": result["counts"]}, sort_keys=True))
