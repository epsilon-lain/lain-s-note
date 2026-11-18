#!/usr/bin/env python3
"""Verify the rank-two family using the pinned official operation semantics."""

import argparse
import hashlib
import json
from pathlib import Path
import sys
import types


OFFICIAL_COMMIT = "a0fd6e6f52d82c93ccc06ab91d81e8fa3678256e"
OFFICIAL_BLOBS = {
    "canon": "93175152acd876566ff60bbfaf7cd90cd8cd2fa9",
    "core": "0a5da5fe6e3e42c96877734dad46b155a69c02b2",
    "stable_core": "af48e32f766dbcd0cad0becd8c1c4edb4358efac",
}
LIMITS = {"max_path_length": 100000, "max_total_relator_length": 10000, "max_work": 5000000}
PARAMETERS = [*range(1, 65), 100, 1000]
REUSE_MOVES = [4, 12, 5, 16]


def load_official(official_root):
    verifier = official_root / "competition/tools/verifier"
    package = "_operation_family_verifier"
    sys.modules[package] = types.ModuleType(package)
    sys.modules[package].__path__ = [str(verifier)]
    for name, expected in OFFICIAL_BLOBS.items():
        path = verifier / f"{name}.py"
        source = path.read_bytes()
        blob = hashlib.sha1(f"blob {len(source)}\0".encode() + source).hexdigest()
        assert blob == expected, f"Unexpected source blob for {path}"
        qualified = f"{package}.{name}"
        module = types.ModuleType(qualified)
        module.__package__ = package
        sys.modules[qualified] = module
        exec(compile(source, str(path), "exec"), module.__dict__)
    return sys.modules[f"{package}.stable_core"]


def check_family(official):
    assert official.core.NUM_MOVES == 14
    official_apply = official.apply_move
    count = 0

    def counted_apply(state, move):
        nonlocal count
        count += 1
        return official_apply(state, move)

    official.apply_move = counted_apply
    results = []
    for length in PARAMETERS:
        donor = (1,) * length
        initial = (donor, donor + (2,) + (-1,) * length)
        before_delete = (donor, (2,))
        assert initial != before_delete
        assert counted_apply(before_delete, 16) == ((donor,), None)
        assert sum(map(len, initial)) <= LIMITS["max_total_relator_length"]

        # Every prefix of lengths 0, 1 and 2. No merging or length pruning.
        hits = [0, 0, 0]
        for first in range(official.core.NUM_MOVES):
            middle, error = counted_apply(initial, first)
            assert error is None
            hits[1] += middle == before_delete
            for second in range(official.core.NUM_MOVES):
                end, error = counted_apply(middle, second)
                assert error is None
                hits[2] += end == before_delete
        assert (hits[1] > 0) == (length == 1)
        assert (hits[2] > 0) == (length <= 2)

        challenge = {
            "challenge_id": f"synthetic-ci-m-family-l{length}",
            "move_spec_version": official.MOVE_SPEC_VERSION,
            "initial_relators": initial,
            "target_relators": (donor,),
        }
        baseline_moves = [11] * length + [16]
        baseline = official.verify(challenge, baseline_moves, official.MOVE_SPEC_VERSION, LIMITS)
        reuse = official.verify(challenge, REUSE_MOVES, official.MOVE_SPEC_VERSION, LIMITS)
        for result, steps, work in [(baseline, length + 1, 2 * length**2 + 4 * length + 1),
                                    (reuse, 4, 9 * length + 4)]:
            assert result["ok"] and result["length"] == steps and result["work"] == work
            assert result["peak_total_relator_length"] == 3 * length + 1
        results.append({"L": length, "prefix_hits_by_depth_0_1_2": hits,
                        "CI": baseline, "four_step_CIM": reuse})

    return {"official_commit": OFFICIAL_COMMIT, "official_blobs": OFFICIAL_BLOBS,
            "cases": len(results), "prefix_sequences_per_case": 211,
            "route_verifications": 2 * len(results), "official_apply": count,
            "four_step_moves": REUSE_MOVES, "results": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--official-root", type=Path, required=True)
    args = parser.parse_args()
    result = check_family(load_official(args.official_root))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
