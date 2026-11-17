#!/usr/bin/env python3
"""Exact cyclic convolution in F_257, with an explicit arithmetic ledger.

The primary ledger generates the constant 3 from 1 with two field additions.
Every multiplication, addition and subtraction costs one, including products
with 1 or 0. Reduction modulo 257 is included in the field operation. Memory
access, copying, allocation, integer indexing, loop control, comparisons, and
input/output are excluded. Fixed task parameters 0, 1, n and p are available.

This compares two specified algorithms, not unrestricted algorithm lower bounds
or measured execution times. Random tests are evidence, not a general proof.
The n=16 basis-pair tests check the complete bilinear map at that fixed size:
both algorithms are bilinear in their inputs and agree on every basis pair.

This recursive implementation uses O(n log n) total allocation across the call
tree and O(n) peak live field-element storage (plus O(log n) recursion depth).
An iterative in-place radix-two implementation can use O(n) field storage for
vectors and twiddles; rearrangements would still fall outside this ledger.
"""

from dataclasses import dataclass
import json
from pathlib import Path
from random import Random

P = 257
LENGTHS = (16, 32, 64)
SEED = 20261117


@dataclass
class FieldOps:
    multiplications: int = 0
    additions: int = 0
    subtractions: int = 0

    def mul(self, a, b):
        self.multiplications += 1
        return a * b % P

    def add(self, a, b):
        self.additions += 1
        return (a + b) % P

    def sub(self, a, b):
        self.subtractions += 1
        return (a - b) % P

    def ledger(self):
        return {
            "multiplications": self.multiplications,
            "additions": self.additions,
            "subtractions": self.subtractions,
            "additions_and_subtractions": self.additions + self.subtractions,
            "total": self.multiplications + self.additions + self.subtractions,
        }


def ledger_difference(after, before):
    return {key: after[key] - before[key] for key in after}


def normalize_inputs(a, b):
    if len(a) != len(b) or len(a) not in LENGTHS:
        raise ValueError("Equal input lengths must be 16, 32, or 64")
    return [x % P for x in a], [x % P for x in b]


def naive_convolution(a, b):
    a, b = normalize_inputs(a, b)
    ops = FieldOps()
    n = len(a)
    result = []
    for k in range(n):
        value = ops.mul(a[0], b[k])
        for j in range(1, n):
            value = ops.add(value, ops.mul(a[j], b[(k - j) % n]))
        result.append(value)
    return result, ops.ledger()


def ntt_convolution(a, b):
    a, b = normalize_inputs(a, b)
    ops = FieldOps()
    n = len(a)
    stages = {}

    def finish_stage(name, before):
        stages[name] = ledger_difference(ops.ledger(), before)

    before = ops.ledger()
    three = ops.add(ops.add(1, 1), 1)
    finish_stage("generate_constant_three", before)

    before = ops.ledger()
    ladder = {1: three}
    value, exponent = three, 1
    for _ in range(8):
        value = ops.mul(value, value)
        exponent *= 2
        ladder[exponent] = value
    assert ladder[128] == P - 1 and ladder[256] == 1
    omega = ladder[256 // n]
    finish_stage("root_ladder_and_verification", before)

    before = ops.ledger()
    inverse_n = n
    # 255 = (11111111)_2: starting at n leaves seven bits,
    # each requiring one square and one multiplication by n.
    for bit in bin(P - 2)[3:]:
        inverse_n = ops.mul(inverse_n, inverse_n)
        if bit == "1":
            inverse_n = ops.mul(inverse_n, n)
    finish_stage("inverse_n_exponentiation", before)

    before = ops.ledger()
    twiddles = [1, omega]
    for _ in range(2, n):
        twiddles.append(ops.mul(twiddles[-1], omega))
    finish_stage("full_twiddle_table", before)

    def transform(vector, inverse=False):
        size = len(vector)
        if size == 1:
            return vector[:]
        even = transform(vector[::2], inverse)
        odd = transform(vector[1::2], inverse)
        result = [0] * size
        for j in range(size // 2):
            index = j * n // size
            if inverse:
                index = (-index) % n
            product = ops.mul(twiddles[index], odd[j])
            result[j] = ops.add(even[j], product)
            result[j + size // 2] = ops.sub(even[j], product)
        return result

    before = ops.ledger()
    transformed_a = transform(a)
    finish_stage("forward_transform_a", before)
    before = ops.ledger()
    transformed_b = transform(b)
    finish_stage("forward_transform_b", before)
    before = ops.ledger()
    products = [ops.mul(x, y) for x, y in zip(transformed_a, transformed_b)]
    finish_stage("pointwise_product", before)
    before = ops.ledger()
    unscaled = transform(products, inverse=True)
    finish_stage("inverse_transform_unscaled", before)
    before = ops.ledger()
    result = [ops.mul(value, inverse_n) for value in unscaled]
    finish_stage("inverse_normalization", before)
    return result, ops.ledger(), stages, {
        "omega": omega,
        "inverse_n": inverse_n,
        "ladder": {str(k): v for k, v in ladder.items()},
    }


def expected_ledgers(n):
    m = n.bit_length() - 1
    naive = {
        "multiplications": n * n,
        "additions": n * (n - 1),
        "subtractions": 0,
        "additions_and_subtractions": n * (n - 1),
        "total": 2 * n * n - n,
    }
    ntt = {
        "multiplications": 3 * n * m // 2 + 3 * n + 20,
        "additions": 3 * n * m // 2 + 2,
        "subtractions": 3 * n * m // 2,
        "additions_and_subtractions": 3 * n * m + 2,
        "total": 9 * n * m // 2 + 3 * n + 22,
    }
    return naive, ntt


def run_checks():
    rng = Random(SEED)
    records = []
    for n in LENGTHS:
        expected_naive, expected_ntt = expected_ledgers(n)
        count = 0
        representative = None

        def check_pair(a, b, expected=None):
            nonlocal count, representative
            direct, direct_ledger = naive_convolution(a, b)
            transformed, ntt_ledger, stages, constants = ntt_convolution(a, b)
            assert direct == transformed
            assert expected is None or direct == expected
            assert direct_ledger == expected_naive
            assert ntt_ledger == expected_ntt
            assert stages["generate_constant_three"]["additions"] == 2
            assert stages["root_ladder_and_verification"]["multiplications"] == 8
            assert stages["inverse_n_exponentiation"]["multiplications"] == 14
            assert stages["full_twiddle_table"]["multiplications"] == n - 2
            for key in ntt_ledger:
                assert sum(stage[key] for stage in stages.values()) == ntt_ledger[key]
            count += 1
            representative = (stages, constants)

        zero = [0] * n
        one = [1] * n
        delta = [1] + [0] * (n - 1)
        ramp = list(range(n))
        cornercases = [
            (zero, zero, zero),
            (zero, ramp, zero),
            (ramp, zero, zero),
            (delta, ramp, ramp),
            (ramp, delta, ramp),
            (one, one, [n % P] * n),
            ([P - 1] * n, one, [(-n) % P] * n),
            ([(-1) ** j for j in range(n)], ramp, None),
            ([P + j for j in range(n)], [-j for j in range(n)], None),
        ]
        for a, b, expected in cornercases:
            check_pair(a, b, expected)

        for _ in range(100):
            check_pair([rng.randrange(P) for _ in range(n)],
                       [rng.randrange(P) for _ in range(n)])

        basis_pairs = 0
        if n == 16:
            for i in range(n):
                for j in range(n):
                    a, b, expected = [0] * n, [0] * n, [0] * n
                    a[i] = b[j] = 1
                    expected[(i + j) % n] = 1
                    check_pair(a, b, expected)
                    basis_pairs += 1

        stages, constants = representative
        records.append({
            "n": n,
            "log2_n": n.bit_length() - 1,
            "status": "passed",
            "structured_pairs": len(cornercases),
            "random_pairs": 100,
            "basis_pairs": basis_pairs,
            "total_pairs": count,
            "constants": constants,
            "naive_ledger": expected_naive,
            "ntt_ledger_primary": expected_ntt,
            "ntt_total_if_three_is_supplied": expected_ntt["total"] - 2,
            "net_field_operation_saving": expected_naive["total"] - expected_ntt["total"],
            "ntt_stage_ledgers": stages,
        })

    return {
        "status": "all checks passed",
        "field_modulus": P,
        "random_seed": SEED,
        "dependencies": "Python standard library only",
        "operation_model": {
            "charged": "Every field multiplication, addition, and subtraction, including multiplication by 1 or 0; modular reduction is part of the operation.",
            "available_task_parameters": "0, 1, n, and fixed field modulus 257",
            "constant_three": "Generated from 1 with two charged additions on each invocation.",
            "excluded": "Memory access/copies/allocation, integer indexing, loop control, comparisons, input/output, and input residue normalization.",
            "normalization_note": "Formal inputs are elements of F_257 represented by canonical residues; signed/out-of-range structured tests also exercise an uncharged input adapter.",
        },
        "formulas": {
            "naive_multiplications": "n^2",
            "naive_additions": "n(n-1)",
            "naive_total": "2n^2-n",
            "ntt_multiplications": "(3/2)n log2(n)+3n+20",
            "ntt_additions_and_subtractions": "3n log2(n)+2",
            "ntt_total_primary": "(9/2)n log2(n)+3n+22",
            "setup_multiplications": "8+14+(n-2)=n+20",
            "setup_additions": "2",
        },
        "correctness_scope": {
            "random_checks": "Deterministic checks of sampled inputs; not an independent general proof.",
            "basis_checks": "For n=16, all 256 pairs of standard basis vectors agree with the exact cyclic-convolution answer. Given bilinearity of both implementations, this establishes equality of the fixed-size maps on all F_257 inputs.",
            "general_identity": "For omega of order n, T_omega(a*b)=T_omega(a) pointwise-multiplied by T_omega(b); inverse is n^(-1) T_(omega^(-1)).",
            "comparison_limit": "Counts compare these two algorithms, not an unrestricted lower bound or a wall-clock benchmark.",
        },
        "memory": "Recursive code: O(n) peak live field-element storage, O(n log n) total allocation, O(log n) recursion depth. Iterative in-place radix-two transforms can retain O(n) field storage including vectors and twiddles; memory work is excluded from the arithmetic ledger.",
        "results": records,
    }


if __name__ == "__main__":
    checks = run_checks()
    output = Path(__file__).with_name("checks.json")
    output.write_text(json.dumps(checks, indent=2) + "\n", encoding="utf-8")
    for row in checks["results"]:
        print(f"n={row['n']}: {row['total_pairs']} pairs passed; "
              f"direct={row['naive_ledger']['total']}, "
              f"NTT={row['ntt_ledger_primary']['total']}, "
              f"saving={row['net_field_operation_saving']}")
