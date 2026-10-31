#!/usr/bin/env python3
"""Exact finite checks of the delta-adic Schur recursion in README section 8.

Only the Python standard library is used. Input matrices have rational polynomial
entries. The recursion uses truncated formal series; an independent permutation
expansion computes each polynomial determinant for comparison. These finite
checks do not establish a theorem about pi or P2848.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import permutations
import json
from math import comb
from pathlib import Path


Q = Fraction


def poly(*coefficients):
    """Ascending powers of delta, with exact rational coefficients."""
    values = [Q(x) for x in coefficients] or [Q(0)]
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return values


def poly_add(a, b):
    values = [Q(0)] * max(len(a), len(b))
    for i, value in enumerate(a):
        values[i] += value
    for i, value in enumerate(b):
        values[i] += value
    return poly(*values)


def poly_mul(a, b):
    values = [Q(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            values[i + j] += x * y
    return poly(*values)


def poly_scale(a, c):
    return poly(*(Q(c) * x for x in a))


def permutation_sign(order):
    inversions = sum(order[i] > order[j]
                     for i in range(len(order))
                     for j in range(i + 1, len(order)))
    return -1 if inversions % 2 else 1


def determinant_by_permutations(matrix):
    """Independent polynomial determinant; no Schur or Gaussian elimination."""
    size = len(matrix)
    result = poly(0)
    for order in permutations(range(size)):
        term = poly(permutation_sign(order))
        for row, column in enumerate(order):
            term = poly_mul(term, matrix[row][column])
        result = poly_add(result, term)
    return result


def zeros(rows, columns):
    return [[Q(0) for _ in range(columns)] for _ in range(rows)]


def identity(size):
    result = zeros(size, size)
    for i in range(size):
        result[i][i] = Q(1)
    return result


def matrix_add(a, b):
    return [[x + y for x, y in zip(row_a, row_b)]
            for row_a, row_b in zip(a, b)]


def matrix_neg(a):
    return [[-x for x in row] for row in a]


def matrix_mul(a, b):
    rows, inner, columns = len(a), len(b), len(b[0])
    assert len(a[0]) == inner
    result = zeros(rows, columns)
    for i in range(rows):
        for k in range(inner):
            if a[i][k]:
                for j in range(columns):
                    result[i][j] += a[i][k] * b[k][j]
    return result


def determinant_constant(matrix):
    a = [row[:] for row in matrix]
    size = len(a)
    result = Q(1)
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if a[row][column]), None)
        if pivot is None:
            return Q(0)
        if pivot != column:
            a[pivot], a[column] = a[column], a[pivot]
            result = -result
        value = a[column][column]
        result *= value
        for row in range(column + 1, size):
            multiplier = a[row][column] / value
            for j in range(column + 1, size):
                a[row][j] -= multiplier * a[column][j]
            a[row][column] = Q(0)
    return result


def inverse_constant(matrix):
    size = len(matrix)
    unit = identity(size)
    a = [row[:] + unit[i] for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next(row for row in range(column, size)
                     if a[row][column])
        a[pivot], a[column] = a[column], a[pivot]
        value = a[column][column]
        a[column] = [x / value for x in a[column]]
        for row in range(size):
            if row != column:
                multiplier = a[row][column]
                a[row] = [x - multiplier * y
                          for x, y in zip(a[row], a[column])]
    return [row[size:] for row in a]


def rank_permutations(matrix):
    """Find permutations whose leading rank-sized block is nonsingular."""
    a = [row[:] for row in matrix]
    size = len(a)
    row_order, column_order = list(range(size)), list(range(size))
    rank = 0
    for index in range(size):
        pivot = next(((i, j)
                      for i in range(index, size)
                      for j in range(index, size) if a[i][j]), None)
        if pivot is None:
            break
        i, j = pivot
        a[index], a[i] = a[i], a[index]
        row_order[index], row_order[i] = row_order[i], row_order[index]
        for row in a:
            row[index], row[j] = row[j], row[index]
        column_order[index], column_order[j] = column_order[j], column_order[index]
        pivot_value = a[index][index]
        for i in range(index + 1, size):
            multiplier = a[i][index] / pivot_value
            for j in range(index + 1, size):
                a[i][j] -= multiplier * a[index][j]
            a[i][index] = Q(0)
        rank += 1
    return rank, row_order, column_order


def polynomial_matrix_to_series(matrix, precision):
    size = len(matrix)
    return [[[matrix[i][j][degree] if degree < len(matrix[i][j]) else Q(0)
              for j in range(size)] for i in range(size)]
            for degree in range(precision + 1)]


def series_product(a, b, precision):
    rows, columns = len(a[0]), len(b[0][0])
    result = [zeros(rows, columns) for _ in range(precision + 1)]
    for degree in range(precision + 1):
        for left_degree in range(degree + 1):
            result[degree] = matrix_add(
                result[degree],
                matrix_mul(a[left_degree], b[degree - left_degree]))
    return result


def inverse_series(a, precision):
    """Y_0=A_0^-1; Y_m=-A_0^-1 sum_{j=1}^m A_j Y_{m-j}."""
    size = len(a[0])
    y = [inverse_constant(a[0])]
    for degree in range(1, precision + 1):
        total = zeros(size, size)
        for j in range(1, degree + 1):
            total = matrix_add(total, matrix_mul(a[j], y[degree - j]))
        y.append(matrix_neg(matrix_mul(y[0], total)))
    expected = [identity(size)] + [zeros(size, size) for _ in range(precision)]
    assert series_product(a, y, precision) == expected
    assert series_product(y, a, precision) == expected
    return y


def constant_is_zero(matrix):
    return all(value == 0 for row in matrix for value in row)


def schur_saturation(matrix):
    """Return the exact determinant valuation and leading coefficient.

    The determinant degree bound is the sum of the maximum entry degree in each
    column. With accumulated valuation e, only coefficients through D-e are
    retained. A recursion that would force e>D proves the determinant is zero.
    No independent determinant calculation is used by this function.
    """
    size = len(matrix)
    if not size or any(len(row) != size for row in matrix):
        raise ValueError("A nonempty square polynomial matrix is required")
    degree_bound = sum(max(len(matrix[row][column]) - 1 for row in range(size))
                       for column in range(size))
    current = polynomial_matrix_to_series(matrix, degree_bound)
    exponent = 0
    leading_factor = Q(1)
    trace = []
    while True:
        size = len(current[0])
        precision = degree_bound - exponent
        assert len(current) == precision + 1
        rank, rows, columns = rank_permutations(current[0])
        step = {
            "size": size, "precision": precision, "rank_at_zero": rank,
            "valuation_before": exponent,
        }
        if rank == size:
            terminal = determinant_constant(current[0])
            leading_factor *= terminal
            step.update(action="full_rank", constant_determinant=str(terminal))
            trace.append(step)
            return {"zero": False, "degree_bound": degree_bound,
                    "valuation": exponent, "leading_coefficient": str(leading_factor),
                    "trace": trace}

        removed = size - rank
        next_exponent = exponent + removed
        step["valuation_added"] = removed
        if next_exponent > degree_bound:
            step.update(action="degree_bound_forces_zero",
                        forced_valuation_at_least=next_exponent)
            trace.append(step)
            return {"zero": True, "degree_bound": degree_bound,
                    "valuation": None, "leading_coefficient": "0", "trace": trace}
        next_precision = degree_bound - next_exponent

        if rank == 0:
            assert constant_is_zero(current[0])
            current = current[1:next_precision + 2]
            step.update(action="divide_whole_matrix", precision_after=next_precision)
        else:
            permuted = [[[coefficient[i][j] for j in columns] for i in rows]
                        for coefficient in current]
            a = [[row[:rank] for row in coefficient[:rank]] for coefficient in permuted]
            b = [[row[rank:] for row in coefficient[:rank]] for coefficient in permuted]
            c = [[row[:rank] for row in coefficient[rank:]] for coefficient in permuted]
            d = [[row[rank:] for row in coefficient[rank:]] for coefficient in permuted]
            pivot_determinant = determinant_constant(a[0])
            sign = permutation_sign(rows) * permutation_sign(columns)
            leading_factor *= sign * pivot_determinant
            inverse = inverse_series(a, precision)
            product = series_product(series_product(c, inverse, precision), b, precision)
            schur = [matrix_add(x, matrix_neg(y)) for x, y in zip(d, product)]
            assert constant_is_zero(schur[0])
            current = schur[1:next_precision + 2]
            step.update(action="schur_then_divide", row_permutation=rows,
                        column_permutation=columns, permutation_sign=sign,
                        pivot_constant_determinant=str(pivot_determinant),
                        precision_after=next_precision,
                        schur_constant_is_zero=True)
        trace.append(step)
        exponent = next_exponent


def columns_to_matrix(columns):
    return [[columns[column][row] for column in range(len(columns))]
            for row in range(len(columns[0]))]


def endpoint_columns(vectors, reduced):
    """U=2, B=1 block, directly from the original or reduced formula."""
    result = []
    for b in range(2):
        column = [poly(0) for _ in range(4)]
        if reduced:
            # R_b=(-1)^b/binom(2b,b) sum_s binom(3+b,b+s+1)binom(s,b)d^(s-b)V_s.
            terms = [(s, s - b, Q((-1) ** b * comb(3 + b, b + s + 1) * comb(s, b),
                                  comb(2 * b, b))) for s in range(b, 3)]
        else:
            terms = [(s, b + 1 + s, Q(comb(3 + b, 2 - s))) for s in range(3)]
        for s, degree, coefficient in terms:
            for row in range(4):
                term = [Q(0)] * degree + [coefficient * vectors[s][row]]
                column[row] = poly_add(column[row], poly(*term))
        result.append(column)
    return result


def cases():
    p = poly
    examples = [
        ("invertible_constant_3", "Nonunit rational leading coefficient and valuation zero.",
         [[p(0), p(Q(1, 2)), p(0)], [p(0), p(0), p(3)], [p(2), p(1), p(1)]]),
        ("nontrivial_schur_inverse_3", "A(delta)=1+delta requires a nonconstant inverse series.",
         [[p(1, 1), p(1), p(0)], [p(1), p(1, 1, 1), p(0, 1)],
          [p(0), p(0, 1), p(0, 0, 1)]]),
        ("permuted_pivot_and_late_zero_ranks_3",
         "A rank-two pivot uses row/column permutations and a nontrivial Schur denominator.",
         [[p(0), p(0), p(0, 0, 1)], [p(0), p(2, 1), p(1)],
          [p(0, 1), p(1), p(3, 1)]]),
        ("rational_nontrivial_schur_3", "Exact rational data with singular constant matrix.",
         [[p(Q(1, 2), 1), p(Q(1, 3)), p(0)],
          [p(2), p(Q(4, 3), 1), p(0, 0, 1)],
          [p(0), p(0, 1), p(0, 0, 0, 1)]]),
        ("mixed_orders_4", "Repeated rank-zero and Schur steps; unequal invariant orders.",
         [[p(0, 1), p(0, 1), p(0), p(0)],
          [p(0), p(0, 0, 1), p(0, 0, 1), p(0)],
          [p(0), p(0), p(0, 0, 0, 1), p(0, 0, 0, 1)],
          [p(0), p(0), p(0), p(0, 0, 0, 0, 1)]]),
        ("identically_singular_3", "Equal columns; degree bound must certify the zero determinant.",
         [[p(1, 1), p(1, 1), p(0, 1)],
          [p(2), p(2), p(1)], [p(0, 1), p(0, 1), p(0, 0, 1)]]),
    ]
    first = [[Q(1), Q(0), Q(0), Q(0)],
             [Q(0), Q(1), Q(0), Q(0)],
             [Q(0), Q(0), Q(1), Q(0)]]
    second = [[Q(1), Q(1), Q(0), Q(0)],
              [Q(0), Q(1), Q(0), Q(0)],
              [Q(0), Q(0), Q(0), Q(1)]]
    for reduced, name in [(False, "two_endpoint_blocks_original_4"),
                          (True, "two_endpoint_blocks_residual_4")]:
        matrix = columns_to_matrix(endpoint_columns(first, reduced)
                                   + endpoint_columns(second, reduced))
        examples.append((name,
                         "Constructed coefficient-vector model of two U=2,B=1 blocks; not a full paper interpolation matrix.",
                         matrix))
    return examples


def verify_case(name, description, matrix):
    result = schur_saturation(matrix)
    determinant = determinant_by_permutations(matrix)
    expected_valuation = next((i for i, value in enumerate(determinant) if value), None)
    expected_leading = determinant[expected_valuation] if expected_valuation is not None else Q(0)
    assert result["valuation"] == expected_valuation, (name, result, determinant)
    assert Q(result["leading_coefficient"]) == expected_leading, (name, result, determinant)
    assert result["zero"] == (expected_valuation is None)
    assert len(determinant) - 1 <= result["degree_bound"]
    return {"name": name, "description": description,
            "dimension": len(matrix),
            "input_polynomial_matrix": [[[str(x) for x in entry] for entry in row] for row in matrix],
            "independent_determinant_coefficients": [str(x) for x in determinant],
            "recursion": result, "matches_independent_determinant": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("saturation_results.json"))
    args = parser.parse_args()
    results = [verify_case(*case) for case in cases()]
    by_name = {result["name"]: result for result in results}
    original = [Q(x) for x in by_name["two_endpoint_blocks_original_4"]["independent_determinant_coefficients"]]
    residual = [Q(x) for x in by_name["two_endpoint_blocks_residual_4"]["independent_determinant_coefficients"]]
    assert original == poly(*([Q(0)] * 8 + residual))
    assert by_name["two_endpoint_blocks_residual_4"]["recursion"]["valuation"] == 2
    report = {
        "status": "all_exact_checks_passed",
        "scope": "Fixed finite rational polynomial matrices only; no theorem about pi or P2848.",
        "arithmetic": "Python standard-library fractions.Fraction; no floating-point operations.",
        "independent_reference": "Permutation expansion with rational polynomial arithmetic.",
        "precision_rule": "Initial D=sum(maximum entry degree in each column); after accumulated order e retain through D-e.",
        "case_count": len(results),
        "endpoint_crosscheck": {"model": "Constructed coefficient vectors satisfying the block expansion; not a full paper interpolation matrix.",
                                "known_first_reduction_order": 8, "additional_order": 2,
                                "total_order": 10, "leading_coefficient": "-9",
                                "exact_identity": "det(original)=delta^8 det(residual)"},
        "cases": results,
    }
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "case_count": len(results),
                      "endpoint_crosscheck": report["endpoint_crosscheck"],
                      "output": str(args.output)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
