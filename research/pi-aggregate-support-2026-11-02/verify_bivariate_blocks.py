"""Exact finite checks for the bivariate normal-profile appendix.

Standard library only. The universal formula is proved in the note.
The values Z used below are formal polynomial specializations.
"""

from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, factorial
from pathlib import Path
import json


def determinant(matrix):
    a = [[F(x) for x in row] for row in matrix]
    out = F(1)
    for k in range(len(a)):
        p = next((i for i in range(k, len(a)) if a[i][k]), None)
        if p is None:
            return F(0)
        if p != k:
            a[k], a[p] = a[p], a[k]
            out = -out
        diagonal = a[k][k]
        out *= diagonal
        for i in range(k + 1, len(a)):
            factor = a[i][k] / diagonal
            for j in range(k + 1, len(a)):
                a[i][j] -= factor * a[k][j]
    return out


def rank(matrix):
    if not matrix or not matrix[0]:
        return 0
    a = [[F(x) for x in row] for row in matrix]
    k = 0
    for c in range(len(a[0])):
        p = next((i for i in range(k, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[k], a[p] = a[p], a[k]
        diagonal = a[k][c]
        for i in range(k + 1, len(a)):
            factor = a[i][c] / diagonal
            for j in range(c + 1, len(a[0])):
                a[i][j] -= factor * a[k][j]
        k += 1
        if k == len(a):
            break
    return k


def multiply(a, b):
    if not a:
        return []
    return [[sum(x * y for x, y in zip(row, col))
             for col in zip(*b)] for row in a]


def transpose(a):
    return [list(row) for row in zip(*a)]


@lru_cache(None)
def exponential_power_coefficient(power, degree):
    """[z^degree](e^z-1)^power, exactly."""
    if degree < power or degree < 0:
        return F(0)
    numerator = sum((-1) ** (power - h) * comb(power, h) * h ** degree
                    for h in range(power + 1))
    return F(numerator, factorial(degree))


def packet(k, v0, v, source_height):
    return [(j, s, l)
            for j in range(k)
            for s in range((source_height - 1) // v0 + 1)
            for l in range((source_height - 1) // v + 1)
            if v0 * s + v * l < source_height]


def columns(w0, w, column_height):
    return [(h, a)
            for h in range(column_height // w0 + 1)
            for a in range(column_height // w + 1)
            if w0 * h + w * a <= column_height]


def raw_entry(row, col, zeta):
    j, s, l = row
    h, a = col
    if l > a:
        return F(0)
    return comb(a, l) * sum(
        F(comb(a - l, d) * h ** (s - d), factorial(s - d))
        * (j * zeta) ** (a - l - d)
        for d in range(min(s, a - l) + 1)
    )


def shifted_entry(row, col, zeta):
    j, s, l = row
    h, a = col
    if l > a:
        return F(0)
    return comb(a, l) * sum(
        comb(a - l, d) * (j * zeta) ** (a - l - d)
        * exponential_power_coefficient(h, s - d)
        for d in range(min(s, a - l) + 1)
    )


def direct_entry(row, col, zeta):
    j, s, l = row
    h, a = col
    if s != h or l > a:
        return F(0)
    return comb(a, l) * (j * zeta) ** (a - l)


def row_change(old, direct):
    j, s, l = old
    jj, r, a = direct
    if j != jj or a < l:
        return F(0)
    return comb(a, l) * exponential_power_coefficient(r, s - a + l)


def column_change(raw, shifted):
    r, b = raw
    h, a = shifted
    if a != b or r > h:
        return 0
    return (-1) ** (h - r) * comb(h, r)


def block_formula(k, w0, w, v0, v, bs, bc):
    values = []
    for h in range((bs - 1) // v0 + 1):
        n = max(0, -((v0 * h - bs) // v))
        d = max(0, 1 + (bc - w0 * h) // w)
        values.append((h, n, d, min(k * n, d)))
    return sum(x[3] for x in values), values


def exact_models():
    records = []
    coordinate_entries = 0
    for k, weights, heights in product(
            (1, 2, 3),
            ((1, 1, 1, 1), (2, 3, 3, 2), (3, 1, 4, 2)),
            ((3, 4), (5, 6))):
        w0, w, v0, v = weights
        bs, bc = heights
        assert v0 >= v
        rows = packet(k, v0, v, bs)
        cols = columns(w0, w, bc)
        left = [[row_change(old, new) for new in rows] for old in rows]
        right = [[column_change(raw, shifted) for shifted in cols]
                 for raw in cols]
        assert determinant(left) == 1
        assert determinant(right) == 1
        expected, blocks = block_formula(k, w0, w, v0, v, bs, bc)
        for zeta in (F(2), F(-1, 2)):
            raw = [[raw_entry(row, col, zeta) for col in cols] for row in rows]
            shifted = [[shifted_entry(row, col, zeta) for col in cols] for row in rows]
            direct = [[direct_entry(row, col, zeta) for col in cols] for row in rows]
            assert multiply(raw, right) == shifted
            assert multiply(left, direct) == shifted
            observed = rank(raw)
            assert observed == rank(direct) == expected
            shifted_gram = multiply(shifted, transpose(shifted))
            direct_gram = multiply(direct, transpose(direct))
            assert determinant(shifted_gram) == determinant(direct_gram)
            coordinate_entries += 2 * len(rows) * len(cols)
            records.append({
                "K": k, "weights": weights, "cutoffs": heights,
                "formal_Z": str(zeta), "rows": len(rows), "columns": len(cols),
                "rank": observed, "blocks": blocks,
            })
    return records, coordinate_entries


def counterexample():
    rows = packet(1, 1, 3, 8)
    cols = columns(2, 2, 8)
    a = [[raw_entry(row, col, F(0)) for col in cols] for row in rows]
    kernel = [F((-1) ** h * factorial(4),
                factorial(h) * factorial(t) * factorial(4 - h - t))
              for h, t in cols]
    assert len(rows) == len(cols) == 15
    assert all(sum(x * y for x, y in zip(row, kernel)) == 0 for row in a)
    assert rank(a) == 14
    minor = [[x for j, x in enumerate(row) if cols[j] != (0, 0)]
             for i, row in enumerate(a) if rows[i] != (0, 0, 1)]
    minor_value = determinant(minor)
    assert minor_value == F(1, 640)
    families = []
    for n in range(1, 9):
        r = packet(1, 1, 3, 2 * n)
        c = columns(2, 2, 2 * n)
        vector = [F((-1) ** h * factorial(n),
                    factorial(h) * factorial(t) * factorial(n - h - t))
                  for h, t in c]
        assert any(vector)
        for row in r:
            assert sum(raw_entry(row, col, F(0)) * coefficient
                       for col, coefficient in zip(c, vector)) == 0
        families.append({"n": n, "rows": len(r), "columns": len(c)})
    return {"rows": 15, "columns": 15, "rank": 14,
            "polynomial": "(T-Y+1)^4", "deleted_row": [0, 0, 1],
            "deleted_column": [0, 0], "minor_determinant": str(minor_value),
            "invisible_polynomial_family": families}


if __name__ == "__main__":
    models, entries = exact_models()
    result = {
        "status": "passed", "exact_models": models,
        "coordinate_identity_entries": entries,
        "counterexample": counterexample(),
        "scope": "Exact finite block ranks under v0>=v, invertible coordinate changes, "
                 "shifted-column Gram equality, and an explicit rank defect outside "
                 "that regime. Rational Z values are formal specializations. "
                 "No critical-family rank formula or pi endpoint is asserted.",
    }
    Path(__file__).with_name("bivariate_block_verification.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": "passed", "exact_models": len(models),
                      "coordinate_identity_entries": entries,
                      "counterexample_rank": result["counterexample"]["rank"],
                      "counterexample_minor": result["counterexample"]["minor_determinant"],
                      "invisible_family_cases": len(result["counterexample"]["invisible_polynomial_family"])},
                     indent=2))
