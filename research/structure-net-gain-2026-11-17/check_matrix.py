#!/usr/bin/env python3
"""Compare two exact algorithms for A^8 (u,v), A=[[a,b],[b,a]].

The theorem holds over any field of characteristic not two. This checker
exhausts F_17^4, with one unit for each scalar +,-,*,/ operation, including
division by the explicitly constructed constant two. Input reads (four)
and final output writes (two) are separately charged. Register accesses,
loop control and comparisons are not additional charges in this model.
These are algorithm counts, not unrestricted lower bounds or timings.
"""
from dataclasses import dataclass
from itertools import product
import json
from pathlib import Path

P = 17


@dataclass
class Ops:
    add: int = 0
    sub: int = 0
    mul: int = 0
    div: int = 0

    def plus(self, a, b):
        self.add += 1
        return (a + b) % P

    def minus(self, a, b):
        self.sub += 1
        return (a - b) % P

    def times(self, a, b):
        self.mul += 1
        return a * b % P

    def divide(self, a, b):
        assert b % P != 0
        self.div += 1
        return a * pow(b, P - 2, P) % P

    def ledger(self):
        arithmetic = self.add + self.sub + self.mul + self.div
        return {'addition': self.add, 'subtraction': self.sub,
                'multiplication': self.mul, 'division': self.div,
                'arithmetic': arithmetic, 'input_reads': 4, 'output_writes': 2,
                'total': arithmetic + 6}


def direct(a, b, u, v):
    ops = Ops()
    for _ in range(8):
        u, v = (ops.plus(ops.times(a, u), ops.times(b, v)),
                ops.plus(ops.times(b, u), ops.times(a, v)))
    return (u, v), ops.ledger()


def transformed(a, b, u, v):
    ops = Ops()
    two = ops.plus(1, 1)
    p, q = ops.plus(a, b), ops.minus(a, b)
    s, d = ops.plus(u, v), ops.minus(u, v)
    for _ in range(3):
        p, q = ops.times(p, p), ops.times(q, q)
    r, t = ops.times(p, s), ops.times(q, d)
    answer = (ops.divide(ops.plus(r, t), two),
              ops.divide(ops.minus(r, t), two))
    return answer, ops.ledger()


def main():
    count = 0
    raw_ledger = new_ledger = None
    for a, b, u, v in product(range(P), repeat=4):
        raw, rc = direct(a, b, u, v)
        new, nc = transformed(a, b, u, v)
        assert raw == new, (a, b, u, v, raw, new)
        assert rc['arithmetic'] == 48 and rc['total'] == 54
        assert nc['arithmetic'] == 17 and nc['total'] == 23
        if count:
            assert raw_ledger == rc and new_ledger == nc
        raw_ledger, new_ledger = rc, nc
        count += 1
    result = {
        'status': 'passed', 'field': 'F_17', 'exhaustive_input_quadruples': count,
        'direct_ledger': raw_ledger, 'transformed_ledger': new_ledger,
        'net_saving_per_instance': 31,
        'constant_two_is_constructed_and_charged': True,
        'division_note': 'Scalar field division is one primitive in BOTH algorithm models; Python modular-inverse implementation is not the counted cost model.',
        'conversion': 'p=a+b,q=a-b,s=u+v,d=u-v (four scalar operations)',
        'restore': 'two sums/differences and two divisions by constructed 2',
        'extra_instance_information': 'none: both two-coordinate maps are bijections',
        'excluded': ['control/index/comparison overhead', 'additional register accesses within scalar instructions', 'program-discovery or proof-discovery cost', 'wall-clock comparison'],
        'general_proof': 'T A T^-1 = diag(a+b,a-b), T=[[1,1],[1,-1]]; char !=2',
        'not_claimed': ['54 is optimal among arbitrary scalar arithmetic algorithms', 'competition improvement', 'new theorem priority', 'native Yohaku execution'],
    }
    Path(__file__).with_name('matrix-checks.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
