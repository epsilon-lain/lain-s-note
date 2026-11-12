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
    scale = 10**40
    checks = {}
    for key, expected in EXPECTED.items():
        bound = report[key]["min_markov_upper_bound"]
        numerator, denominator = int(bound["upper_num"]), int(bound["upper_den"])
        if numerator < 0 or denominator <= 0:
            raise SystemExit("Invalid nonnegative rational bound")
        lower, remainder = divmod(numerator * scale, denominator)
        upper = lower + bool(remainder)
        display = f"{upper // scale}.{upper % scale:040d}"
        floor_display = f"{lower // scale}.{lower % scale:040d}"
        if display != expected:
            raise SystemExit(f"{key}: unexpected upper display")
        if bound["decimal_40"] != floor_display or remainder == 0:
            raise SystemExit(f"{key}: expected historical strict downward truncation")
        if not lower * denominator < numerator * scale <= upper * denominator:
            raise SystemExit(f"{key}: directed-rounding comparison failed")
        checks[key] = {
            "safe_upper_40": display,
            "historical_decimal_is_strict_lower_truncation": True,
            "integer_arithmetic_verified": True,
        }
    print(json.dumps({
        "artifact_sha256": fingerprint,
        "checks": checks,
        "scope": "Exact rational display checks, not independent pi certification",
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
