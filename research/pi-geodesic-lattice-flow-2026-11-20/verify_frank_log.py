"""Identify batches requiring a scan for a new normalized-error record.

The input is the original public Zuuv log from Max Frank:
https://drive.google.com/file/d/1rf5VjABLlb4rRd3aI6LbZQ23vhU7TnhB/view

Usage: python3 verify_frank_log.py zuuv_log.txt
Batch exclusions are conditional on the reported maxima.
"""
import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re


LOG_SHA256 = "af850feac7f7a433e9e59255d6bfdaefe53c38f96717bc3f91ea20ebef398a80"
SEED_DIGIT = 845_188_380_288
SEED_INDEX = 383_059_989_914  # One-based from Frank's README.
PUBLIC_LAST_INDEX = 653_520_000_000
ITERATION = re.compile(
    r"^Iteration ([0-9]+)\.\.\.\n---\n"
    r"Computing terms\.\.\.\s+([0-9]+) terms \| [^\n]+ largest: ([0-9]+)$",
    re.MULTILINE,
)


@dataclass(frozen=True)
class Batch:
    iteration: int
    first_index: int
    last_index: int
    maximum: int


def read_batches(path):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != LOG_SHA256:
        raise ValueError("input does not match the public log byte for byte")
    log = data.decode("utf-8")
    batches = []
    position = 0
    for iteration, count, maximum in ITERATION.findall(log):
        iteration, count, maximum = map(int, (iteration, count, maximum))
        if iteration != len(batches) + 1 or count <= 0 or maximum <= 0:
            raise ValueError("invalid iteration number, size or maximum")
        batches.append(Batch(iteration, position + 1, position + count, maximum))
        position += count
    if len(batches) != 54 or position != 657_290_105_207:
        raise ValueError("unexpected number of iterations or total term count")
    return batches


def analyze(batches):
    relevant = []
    skipped = []
    for batch in batches:
        first = max(batch.first_index, SEED_INDEX + 1)
        last = min(batch.last_index, PUBLIC_LAST_INDEX)
        if first > last:
            continue
        if batch.maximum < SEED_DIGIT - 1:
            skipped.append(batch)
        else:
            relevant.append({
                "iteration": batch.iteration,
                "logged_first": batch.first_index,
                "logged_last": batch.last_index,
                "reported_maximum": batch.maximum,
                "first_to_scan": first,
                "last_to_scan": last,
                "terms_to_scan": last - first + 1,
            })
    return {
        "index_basis": "one-based partial-quotient index, counting a_0 as 1",
        "seed_index": SEED_INDEX,
        "public_last": PUBLIC_LAST_INDEX,
        "possible_batches": relevant,
        "skipped_post_seed_batches": [batch.iteration for batch in skipped],
        "largest_skipped_maximum": max(batch.maximum for batch in skipped),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    result = analyze(read_batches(parser.parse_args().log))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
