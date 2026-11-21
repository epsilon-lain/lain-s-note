"""Compare q^2*abs(pi-p/q) via exact local continued-fraction windows.

An input JSON object has "seed" and "candidates". The latter is a list.
Each window has "index" (zero-based position of the large partial quotient,
with pi's initial 3 at index 0), "digit", and the two lists
"left_nearest_first" and "right_nearest_first". Never include a_0 on the left.
The fraction being tested is the convergent at index-1. Supply only positive
partial quotients. Each comparison is conditional on the supplied digits.

Example invocation: python3 verify_local_cf_records.py windows.json
Scanner output: python3 verify_local_cf_records.py windows.jsonl --seed-index 383059989913
Built-in synthetic tests: python3 verify_local_cf_records.py --self-test
"""
import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path


def tail_interval(digits, complete=False):
    """Enclose [0; digits..., unknown positive partial quotients]."""
    lower, upper = Fraction(0), Fraction(0 if complete else 1)
    for digit in reversed(digits):
        lower, upper = 1 / (digit + upper), 1 / (digit + lower)
    return lower, upper


@dataclass(frozen=True)
class Window:
    index: int
    digit: int
    left_nearest_first: tuple
    right_nearest_first: tuple

    def __post_init__(self):
        if type(self.index) is not int or self.index < 1:
            raise ValueError("index must be a positive integer")
        if len(self.left_nearest_first) >= self.index:
            raise ValueError("the left window must exclude a_0")
        digits = (self.digit, *self.left_nearest_first, *self.right_nearest_first)
        if any(type(digit) is not int or digit < 1 for digit in digits):
            raise ValueError("partial quotients must be positive integers")

    def score_interval(self):
        left_complete = len(self.left_nearest_first) == self.index - 1
        left_low, left_high = tail_interval(self.left_nearest_first, left_complete)
        right_low, right_high = tail_interval(self.right_nearest_first)
        return self.digit + left_low + right_low, self.digit + left_high + right_high

    def indexed_digits(self):
        left = ((self.index - offset, digit)
                for offset, digit in enumerate(self.left_nearest_first, 1))
        right = ((self.index + offset, digit)
                 for offset, digit in enumerate(self.right_nearest_first, 1))
        return dict([(self.index, self.digit), *left, *right])


def compare(seed, candidate):
    """Return the result for the candidate under the normalized error criterion."""
    if seed.index == candidate.index:
        raise ValueError("seed and candidate must represent different convergents")
    seed_digits = seed.indexed_digits()
    for index, digit in candidate.indexed_digits().items():
        if index in seed_digits and seed_digits[index] != digit:
            raise ValueError(f"inconsistent partial quotients at index {index}")
    seed_low, seed_high = seed.score_interval()
    candidate_low, candidate_high = candidate.score_interval()
    # The unknown right tail is irrational, so both bounds are strict.
    if candidate_low >= seed_high:
        return "breaks_record"
    if candidate_high <= seed_low:
        return "below_record"
    return "need_more_neighbor_digits"


def certificate(seed, candidates):
    def interval(window):
        return [str(endpoint) for endpoint in window.score_interval()]

    reports = [
        {
            "index": candidate.index,
            "digit": candidate.digit,
            "reciprocal_normalized_error": interval(candidate),
            "comparison": compare(seed, candidate),
        }
        for candidate in sorted(candidates, key=lambda window: window.index)
    ]
    first_open = next((report for report in reports
                       if report["index"] > seed.index
                       and report["comparison"] != "below_record"), None)
    return {
        "seed_index": seed.index,
        "seed_reciprocal_normalized_error": interval(seed),
        "first_unexcluded_later_candidate": first_open,
        "candidates": reports,
    }


def self_test():
    giant = 845188380288
    # The three windows from one finite prefix. The latter two digits improve
    # the normalized error for every irrational continuation of that prefix.
    prefix = [3, 100, giant, 100, 7, 100, 1, giant - 1, 1, 100,
              5, 100, 1, giant, 1, 100]

    def window(index, radius):
        return Window(index, prefix[index],
                      tuple(reversed(prefix[max(1, index - radius):index])),
                      tuple(prefix[index + 1:index + radius + 1]))

    seed = window(2, 2)
    lower_digit = window(7, 2)
    same_digit = window(13, 2)
    assert compare(seed, lower_digit) == "breaks_record"
    assert compare(seed, same_digit) == "breaks_record"

    high_score_seed = Window(100, giant, (1, 100), (1, 100))
    plus_one_candidate = Window(200, giant + 1, (100,), (100,))
    assert compare(high_score_seed, plus_one_candidate) == "below_record"
    assert compare(seed, Window(200, giant - 2, (), ())) == "below_record"
    assert compare(seed, Window(201, giant + 2, (), ())) == "breaks_record"
    assert compare(seed, Window(202, giant, (), ())) == "need_more_neighbor_digits"

    # The historical a_4=292 window. This tests the sign and inversion of
    # the scoring formula independently with exact rational pi endpoints.
    actual_seed = Window(4, 292, (1, 15, 7), (1, 1, 1, 2, 1, 3, 1, 14))
    pi_low = Fraction(314159265358979323846, 10**20)
    pi_high = pi_low + Fraction(1, 10**20)
    rational = Fraction(355, 113)
    error_low = 113**2 * (rational - pi_high)
    error_high = 113**2 * (rational - pi_low)
    score_low, score_high = actual_seed.score_interval()
    assert score_low < 1 / error_high < 1 / error_low < score_high

    print(json.dumps(certificate(seed, [lower_digit, same_digit]), indent=2))
    print("Seven interval and comparison fixtures passed.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, nargs="?")
    parser.add_argument("--seed-index", type=int, help="read scanner JSONL with this seed index")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    elif args.input is None:
        parser.error("specify an input JSON file or --self-test")
    elif args.seed_index is not None:
        windows = [Window(**json.loads(line))
                   for line in args.input.read_text(encoding="utf-8").splitlines()]
        by_index = {window.index: window for window in windows}
        if len(by_index) != len(windows) or args.seed_index not in by_index:
            parser.error("expected unique indices and the specified seed in the scanner output")
        seed = by_index[args.seed_index]
        candidates = [window for window in windows if window.index > seed.index]
        print(json.dumps(certificate(seed, candidates), indent=2))
    else:
        data = json.loads(args.input.read_text(encoding="utf-8"))
        seed = Window(**data["seed"])
        candidates = [Window(**candidate) for candidate in data["candidates"]]
        print(json.dumps(certificate(seed, candidates), indent=2))


if __name__ == "__main__":
    main()
