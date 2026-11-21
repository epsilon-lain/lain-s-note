"""Replay the finite window and scope audit; the full archive scan is an input."""
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent


def tail_interval(digits):
    """Range of [0; digits..., remaining] using remainder in [0,1]."""
    low, high = Fraction(0), Fraction(1)
    for digit in reversed(digits):
        low, high = 1 / (digit + high), 1 / (digit + low)
    return low, high


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("maxima", type=Path, help="frank_iteration_maxima.json")
    args = parser.parse_args()
    raw = (HERE / "a1_raw.txt").read_bytes()
    window_bytes = (HERE / "a2_window.json").read_bytes()
    witness = json.loads(window_bytes)
    certificate = json.loads((HERE / "a2_exclusion.json").read_bytes())
    batches = json.loads(args.maxima.read_bytes())["batches"]
    assert hashlib.sha256(raw).hexdigest() == "96e50dd3cfff35b9dd511a28e0bfb0aeba013e557d8cb6b2851ee144f0efccbd"
    assert hashlib.sha256(window_bytes).hexdigest() == certificate["authenticated_sources"]["windows_sha256"]
    left = witness["left_nearest_first"]
    right = witness["right_nearest_first"]
    seed_digit, seed_index = witness["digit"], witness["index"]
    assert len(left) == len(right) == 64
    assert list(map(int, raw.splitlines())) == left[::-1] + [seed_digit] + right
    assert certificate["threshold_scan"]["hits"] == [witness]
    assert certificate["threshold_scan"]["all_volume_hit_count"] == 1
    left_low, left_high = tail_interval(left)
    right_low, right_high = tail_interval(right)
    low = seed_digit + left_low + right_low
    high = seed_digit + left_high + right_high
    endpoints = list(map(Fraction, certificate["threshold_scan"]["seed_reciprocal_error_interval_exact"]))
    assert endpoints == [low, high]
    assert Fraction("1.5576774226367261722") < low - seed_digit < high - seed_digit < Fraction("1.5576774226367261723")
    assert seed_digit == 845188380288 and seed_index == 383059989913
    assert left[:3] == [1, 1, 1] and right[:2] == [1, 11]
    assert low > seed_digit + Fraction(91, 60)
    public_count = 653520000000
    scan_end = 0
    for i, batch in enumerate(batches, 1):
        assert batch["iteration"] == i and batch["one_based_start"] == scan_end + 1
        first = scan_end
        scan_end += batch["term_count"]
        if i == 32:
            assert first == 377333339473
            assert seed_index - first + 1 == 5726650441
            assert scan_end - 1 - seed_index == 6445380663
            assert batch["maximum"] == seed_digit
        elif first < public_count:
            assert batch["maximum"] <= seed_digit - 2
    assert len(batches) == 54 and scan_end >= public_count
    count = public_count - 1 - seed_index
    assert count == certificate["exact_exclusion"]["later_convergent_count"]
    assert [seed_index, public_count - 2] == certificate["exact_exclusion"]["later_convergent_n_inclusive"]
    print(f"PASS. Exact seed endpoints matched; {count} post-seed convergents excluded given the source scan and maxima.")
    print(f"First unresolved index: {public_count - 1}; first unknown partial quotient: a_{public_count}.")
    print("A-1 excluded. For a=A, both adjoining coefficients must be 1.")


if __name__ == "__main__":
    main()
