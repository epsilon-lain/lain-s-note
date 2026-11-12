"""Check directed decimal displays against a frozen exact-rational artifact.

This verifies integer arithmetic and file provenance only, not a certificate for
pi or the continued-fraction computation that produced the artifact.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


KNOWN_ARTIFACTS = {
    "324a5a829cd5b95bd972e80cf904a69ee1b7e3bed3ea44e5653c8a663200f4d4",
    "c05395aad616260699357900762832d0d63cf4d2d71b005c5951f5bdf5deae3a",
}
EXPECTED = {
    "runA": "0.0000481305612462894046281170790752902943",
    "runB": "0.0000127178573727574566035844693935801275",
}
EXPECTED_A_18 = "0.000048130561246290"


def check_grid_upper(numerator: int, denominator: int, places: int,
                     expected: str) -> dict:
    """Verify the exact fraction lies between adjacent decimal grid points."""
    scale = 10**places
    lower, remainder = divmod(numerator * scale, denominator)
    upper = lower + bool(remainder)
    display = f"{upper // scale}.{upper % scale:0{places}d}"
    floor_display = f"{lower // scale}.{lower % scale:0{places}d}"
    if display != expected or remainder == 0:
        raise SystemExit(f"Unexpected nonterminating {places}-place display")
    if not (upper - 1) * denominator < numerator * scale < upper * denominator:
        raise SystemExit(f"{places}-place strict bound/minimality comparison failed")
    return {
        "places": places,
        "lower_truncation": floor_display,
        "safe_upper": display,
        "strict_upper_of_saved_fraction": True,
        "previous_grid_point_strictly_below_saved_fraction": True,
        "minimal_grid_upper": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path)
    args = parser.parse_args()
    raw = args.artifact.read_bytes()
    fingerprint = hashlib.sha256(raw).hexdigest()
    if fingerprint not in KNOWN_ARTIFACTS:
        raise SystemExit("Unrecognized artifact SHA-256; no known-result check performed")
    # Known artifacts contain decimal integers longer than Python's default cap.
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    report = json.loads(raw)
    checks = {}
    for key, expected in EXPECTED.items():
        bound = report[key]["min_markov_upper_bound"]
        numerator, denominator = int(bound["upper_num"]), int(bound["upper_den"])
        if numerator < 0 or denominator <= 0:
            raise SystemExit("Invalid nonnegative rational bound")
        check_40 = check_grid_upper(numerator, denominator, 40, expected)
        if bound["decimal_40"] != check_40["lower_truncation"]:
            raise SystemExit(f"{key}: expected historical strict downward truncation")
        checks[key] = {
            "safe_upper_40": check_40["safe_upper"],
            "historical_decimal_is_strict_lower_truncation": True,
            "integer_arithmetic_verified": True,
            "decimal_grid_40": check_40,
        }
        if key == "runA":
            checks[key]["v04c_decimal_grid_18"] = check_grid_upper(
                numerator, denominator, 18, EXPECTED_A_18)
    print(json.dumps({
        "artifact_sha256": fingerprint,
        "checks": checks,
        "scope": "Exact rational display checks, not independent pi certification",
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
