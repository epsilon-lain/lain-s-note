"""Verify the fixed-donor free-word formulas and the rank-three empty-target examples."""

import argparse
import json
from pathlib import Path
import sys


def reduced_words(letters, max_length):
    layer = [()]
    for _ in range(max_length):
        layer = [(*word, letter) for word in layer for letter in letters
                 if not word or letter != -word[-1]]
        yield from layer


def run(official_root):
    sys.path.insert(0, str(official_root / "competition/tools"))
    from verifier import stable_core as official

    manifest = json.loads((official_root / "competition/tools/verifier/data/manifest.json").read_bytes())
    conjugate = {(row["relator"], row["conjugator"]): row["id"]
                 for row in official.MOVE_TABLE if row["category"] == "conjugation"}
    allowed = (1, 4, 5, *(conjugate[1, letter] for letter in (1, -1, 2, -2, 3, -3)))
    short_route = [4, conjugate[1, 2], 5, 16]
    assert short_route == [4, 12, 5, 16]

    def apply(state, move):
        result, reason = official.apply_move(state, move)
        assert result is not None and reason is None
        return result

    def relabel(word):
        return tuple(2 if letter == 3 else -2 if letter == -3 else letter for letter in word)

    def verify(initial, target, route, name):
        challenge = {"challenge_id": name, "move_spec_version": official.MOVE_SPEC_VERSION,
                     "initial_relators": initial, "target_relators": target}
        verdict = official.verify(challenge, route, official.MOVE_SPEC_VERSION, manifest["limits"])
        assert verdict["ok"], (challenge, route, verdict)
        return verdict

    words = [*reduced_words((1, -1, 3, -3), 4),
             (3,) * 6 + (1,) + (-3,) * 6, (1, 3, -1, -3) * 4]
    examples = []
    short_prefix_count = 0
    for index, donor in enumerate(words):
        length = len(donor)
        initial = (donor, (*donor, 2, *(-letter for letter in reversed(donor))), (3,))
        before_delete = (donor, (2,), (3,))
        target = (relabel(donor), (2,))
        assert apply(before_delete, 16) == target
        penalty = sum(map(len, initial[2:]))
        inverse_conjugations = [conjugate[1, -letter] for letter in donor]
        direct = [*inverse_conjugations, 16]
        ci = verify(initial, target, direct, f"family-word-{index}-ci")
        cim = verify(initial, target, short_route, f"family-word-{index}-cim")
        assert (ci["length"], ci["work"]) == (length + 1, 2*length**2 + 4*length + 1 + penalty*(length+2))
        assert (cim["length"], cim["work"]) == (4, 9*length+4 + 5*penalty)
        assert ci["peak_total_relator_length"] == cim["peak_total_relator_length"] == 3*length+1+penalty
        hits = [initial == before_delete, 0, 0]
        for first in allowed:
            middle = apply(initial, first)
            hits[1] += middle == before_delete
            for second in allowed:
                hits[2] += apply(middle, second) == before_delete
        short_prefix_count += 1 + len(allowed) + len(allowed)**2
        assert not hits[0] and (hits[1] > 0) == (length == 1) and (hits[2] > 0) == (length <= 2)
        if donor == (1, 3, -1, -3):
            examples.append({"donor": donor, "zero_exponent_example": True, "ci": ci, "cim": cim})

    empty = []
    for k in (1, 2, 4):
        donor = (3,) * k + (1,) + (-3,) * k
        initial = (donor, (*donor, 2, *(-letter for letter in reversed(donor))), (3,))
        ci_prefix = [*(conjugate[1, -letter] for letter in donor), 16]
        suffix = [conjugate[0, -2]] * k + [16, 15]
        ci = verify(initial, (), [*ci_prefix, *suffix], f"family-empty-{k}-ci")
        cim = verify(initial, (), [*short_route, *suffix], f"family-empty-{k}-cim")
        assert ci["length"] == 3*k+4 and cim["length"] == k+6
        empty.append({"k": k, "ci": ci, "cim": cim})
    return {"word_cases": len(words), "short_prefix_sequences": short_prefix_count,
            "exact_delete_route_checks": 2*len(words), "sample_word": examples,
            "empty_target_examples": empty}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--official-root", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.official_root.resolve()), indent=2))


if __name__ == "__main__":
    main()
