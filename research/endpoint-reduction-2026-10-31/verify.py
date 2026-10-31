#!/usr/bin/env python3
"""Finite exact checks for the two-centre endpoint column reduction.

Uses Python's standard library only.  These tests check formulas and their
implementation on finite examples; they do not prove the general identities,
the original interpolation hypotheses, or any assertion about P2848.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import permutations
import json
from math import comb, gcd
from pathlib import Path


def lcm(a: int, b: int) -> int:
    return abs(a * b) // gcd(a, b) if a and b else 0


def choose(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


def coefficient(u: int, a: int, b: int) -> F:
    """Coefficient c[a,b] in the unit upper-triangular column change."""
    if not 0 <= a <= b <= u:
        raise ValueError("require 0 <= a <= b <= u")
    return F((-1) ** (b + a) * comb(b + a, a) * comb(u + b + 1, b - a),
             comb(2 * b, b))


def transformed_closed(u: int, s: int, b: int) -> F:
    """Closed coefficient after transformation (not including delta powers)."""
    if not 0 <= s <= u or not 0 <= b <= u:
        raise ValueError("require 0 <= s,b <= u")
    if s < b:
        return F(0)
    return F((-1) ** b * comb(s, b) * comb(u + b + 1, u - s),
             comb(2 * b, b))


def transformed_direct(u: int, s: int, b: int) -> F:
    return sum(
        (coefficient(u, a, b) * choose(u + 1 + a, u - s)
         for a in range(b + 1)),
        F(0),
    )


def coefficient_matrix(u: int, bmax: int) -> list[list[F]]:
    return [[F(choose(u + 1 + b, u - s)) for b in range(bmax + 1)]
            for s in range(u + 1)]


def canonical_column_elimination(u: int, bmax: int) -> list[list[F]]:
    """Find the unit upper triangular change independently by elimination."""
    matrix = coefficient_matrix(u, bmax)
    change = [[F(a == b) for b in range(bmax + 1)]
              for a in range(bmax + 1)]
    for b in range(bmax + 1):
        for a in range(b):
            pivot = matrix[a][a]
            if not pivot:
                raise AssertionError("unexpected zero canonical pivot")
            multiplier = matrix[a][b] / pivot
            for s in range(u + 1):
                matrix[s][b] -= multiplier * matrix[s][a]
            for j in range(bmax + 1):
                change[j][b] -= multiplier * change[j][a]
    return change


def determinant(matrix):
    """Exact Gaussian elimination, valid for Fraction and Gaussian below."""
    n = len(matrix)
    if not n:
        return F(1)
    a = [row[:] for row in matrix]
    result = a[0][0] * 0 + 1
    for k in range(n):
        pivot_row = next((r for r in range(k, n) if a[r][k]), None)
        if pivot_row is None:
            return result * 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            result = -result
        pivot = a[k][k]
        result *= pivot
        for r in range(k + 1, n):
            multiplier = a[r][k] / pivot
            for c in range(k + 1, n):
                a[r][c] -= multiplier * a[k][c]
            a[r][k] = a[r][k] * 0
    return result


def polynomial_determinant(u: int, bmax: int) -> dict[int, F]:
    """Expand the small S block as an exact polynomial in delta."""
    n = bmax + 1
    terms: dict[int, F] = {}
    for perm in permutations(range(n)):
        inversions = sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
        value = F((-1) ** inversions)
        exponent = 0
        for s, b in enumerate(perm):
            value *= choose(u + 1 + b, u - s)
            exponent += 1 + b + s
        terms[exponent] = terms.get(exponent, F(0)) + value
    return {power: value for power, value in terms.items() if value}


@dataclass(frozen=True)
class Gaussian:
    real: F = F(0)
    imag: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "real", F(self.real))
        object.__setattr__(self, "imag", F(self.imag))

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Gaussian) else Gaussian(F(value))

    def __bool__(self):
        return bool(self.real or self.imag)

    def __add__(self, other):
        other = self.coerce(other)
        return Gaussian(self.real + other.real, self.imag + other.imag)

    __radd__ = __add__

    def __neg__(self):
        return Gaussian(-self.real, -self.imag)

    def __sub__(self, other):
        return self + -self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        return Gaussian(self.real * other.real - self.imag * other.imag,
                        self.real * other.imag + self.imag * other.real)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        denominator = other.real * other.real + other.imag * other.imag
        if not denominator:
            raise ZeroDivisionError
        return Gaussian((self.real * other.real + self.imag * other.imag) / denominator,
                        (self.imag * other.real - self.real * other.imag) / denominator)

    def __pow__(self, exponent: int):
        if exponent < 0:
            return (Gaussian(1) / self) ** (-exponent)
        result = Gaussian(1)
        base = self
        while exponent:
            if exponent & 1:
                result *= base
            base *= base
            exponent //= 2
        return result

    def is_integral(self):
        return self.real.denominator == self.imag.denominator == 1


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, Gaussian):
        return {"real": str(value.real), "imag": str(value.imag)}
    raise TypeError(type(value).__name__)


def tensor_block(u: int, amax: int, bmax: int, r: Gaussian, delta: Gaussian):
    """Rows (i,s), columns (a,b); s indexes derivative degree U-s."""
    return [[
        (Gaussian(choose(a, i)) * r ** (a - i)
         * choose(u + 1 + b, u - s) * delta ** (1 + b + s))
        if a >= i else Gaussian()
        for a in range(amax + 1) for b in range(bmax + 1)
    ] for i in range(amax + 1) for s in range(bmax + 1)]


def reduced_column(u: int, b: int, z, vectors):
    """Return R_b(z)=sum_{s>=b} T[s,b] z^(s-b) V_s.

    The vectors are the full coefficient vectors V_0,...,V_U.  Their
    components may, for example, be indexed by the original (j,k,beta) rows.
    This routine does not assume those components form a small tensor block.
    """
    if len(vectors) != u + 1 or not 0 <= b <= u:
        raise ValueError("need U+1 vectors and 0 <= b <= U")
    size = len(vectors[0])
    if any(len(vector) != size for vector in vectors):
        raise ValueError("coefficient vectors must have equal lengths")
    z = Gaussian.coerce(z)
    result = [Gaussian() for _ in range(size)]
    for s in range(b, u + 1):
        scalar = transformed_closed(u, s, b) * z ** (s - b)
        for index, value in enumerate(vectors[s]):
            result[index] += scalar * value
    return result


# Sparse polynomials in (t,u1,u2), with weights (1,4,8).  All arithmetic is
# exact and terms of weight >=32 are discarded immediately during products.
def polynomial_add(left, right):
    result = left.copy()
    for exponent, value in right.items():
        new_value = result.get(exponent, Gaussian()) + value
        if new_value:
            result[exponent] = new_value
        else:
            result.pop(exponent, None)
    return result


def polynomial_product(left, right, cutoff=32):
    result = {}
    for e1, c1 in left.items():
        for e2, c2 in right.items():
            exponent = tuple(a + b for a, b in zip(e1, e2))
            if exponent[0] + 4 * exponent[1] + 8 * exponent[2] >= cutoff:
                continue
            value = result.get(exponent, Gaussian()) + c1 * c2
            if value:
                result[exponent] = value
            else:
                result.pop(exponent, None)
    return result


def polynomial_powers(poly, maximum, cutoff=32):
    result = [{(0, 0, 0): Gaussian(1)}]
    for _ in range(maximum):
        result.append(polynomial_product(result[-1], poly, cutoff))
    return result


def check_full_row_block():
    """Integrate the reduction into one actual two-centre coefficient pack.

    This checks the original local row coefficients, unequal coordinate
    weights, unequal logarithm truncations, and row-dependent denominator
    clearing. It does not verify the global interpolation separation theorem.
    """
    u, amax, bmax, remaining = 7, 1, 1, 32
    q1, q2, denominator_lcm = 2, 3, 6
    r1, r2 = Gaussian(0, 1), Gaussian(0, F(4, 3))
    delta = r2 - r1
    d = denominator_lcm * delta
    g1 = {(k, 0, 0): Gaussian(F((-1) ** (k + 1), k)) for k in range(1, 10)}
    g2 = {(k, 0, 0): Gaussian(F((-1) ** (k + 1), k)) for k in range(1, 20)}
    g_difference = polynomial_add(g2, {e: -c for e, c in g1.items()})
    base_difference = {(0, 1, 0): Gaussian(-1), (0, 0, 1): Gaussian(1)}
    v = polynomial_add(base_difference, g_difference)
    v_powers = polynomial_powers(v, u + 1 + bmax)
    no_tail_powers = polynomial_powers(base_difference, u)
    x_without_center = polynomial_add(g1, {(0, 1, 0): Gaussian(1)})
    x_polynomials = [x_without_center,
                     polynomial_add(x_without_center, {(0, 0, 0): r1})]
    difference_polynomials = [v, polynomial_add(v, {(0, 0, 0): delta})]
    differences_powers = [polynomial_powers(poly, u + 1 + bmax)
                          for poly in difference_polynomials]
    rows = [(j, k, beta1, beta2)
            for j in range(2)
            for beta2 in range(4)
            for beta1 in range(8)
            for k in range(remaining)
            if k + 4 * beta1 + 8 * beta2 < remaining]
    assert len(rows) == 560

    def pack(poly, supported_j=1):
        return [poly.get((k, beta1, beta2), Gaussian()) if j == supported_j else Gaussian()
                for j, k, beta1, beta2 in rows]

    original_columns = {}
    coefficient_vectors = {}
    tail_effects = []
    counts = {"original_reconstruction_entries": 0,
              "transformed_reconstruction_entries": 0,
              "cleared_V_entries": 0,
              "cleared_original_entries": 0,
              "cleared_transformed_entries": 0,
              "nonzero_V_entries": 0,
              "nonzero_tail_difference_entries": 0}

    for a in range(amax + 1):
        xpower = [{(0, 0, 0): Gaussian(1)} if a == 0 else xp for xp in x_polynomials]
        for b in range(bmax + 1):
            polynomials = [polynomial_product(xpower[j], differences_powers[j][u + 1 + b])
                           for j in range(2)]
            assert not polynomials[0], "special columns must vanish at j=0"
            original_columns[a, b] = [polynomials[j].get((k, beta1, beta2), Gaussian())
                                       for j, k, beta1, beta2 in rows]
        vectors = []
        for s in range(u + 1):
            actual_poly = polynomial_product(xpower[1], v_powers[u - s])
            without_tail = polynomial_product(xpower[1], no_tail_powers[u - s])
            actual_vector, flat_vector = pack(actual_poly), pack(without_tail)
            vectors.append(actual_vector)
            counts["nonzero_V_entries"] += sum(bool(value) for value in actual_vector)
            differences = [i for i, (actual, flat) in enumerate(zip(actual_vector, flat_vector))
                           if actual != flat]
            counts["nonzero_tail_difference_entries"] += len(differences)
            if differences:
                first = differences[0]
                tail_effects.append({"a": a, "s": s, "changed_entries": len(differences),
                                     "first_row_j_k_beta1_beta2": rows[first],
                                     "coefficient_with_tail": actual_vector[first],
                                     "coefficient_without_tail": flat_vector[first]})
        coefficient_vectors[a] = vectors
        for b in range(bmax + 1):
            reconstructed = [Gaussian() for _ in rows]
            for s in range(u + 1):
                scalar = choose(u + 1 + b, u - s) * delta ** (1 + b + s)
                for i, value in enumerate(vectors[s]):
                    reconstructed[i] += scalar * value
            assert reconstructed == original_columns[a, b]
            counts["original_reconstruction_entries"] += len(rows)
            reduced = reduced_column(u, b, delta, vectors)
            transformed = [sum((coefficient(u, j, b) * delta ** (b - j)
                                 * original_columns[a, j][index]
                                 for j in range(b + 1)), Gaussian())
                           for index in range(len(rows))]
            assert transformed == [delta ** (2 * b + 1) * value for value in reduced]
            counts["transformed_reconstruction_entries"] += len(rows)

    lcm9 = lcm19 = 1
    for k in range(1, 20):
        lcm19 = lcm(lcm19, k)
        if k <= 9:
            lcm9 = lcm(lcm9, k)
    truncation_digits = []
    for a in range(amax + 1):
        truncation_denominator = lcm9 ** (a + u + 2) * lcm19 ** (u + 2)
        truncation_digits.append(len(str(truncation_denominator)))
        cleared_vectors = []
        for s in range(u + 1):
            cleared_vector = []
            for value, (_, _, beta1, beta2) in zip(coefficient_vectors[a][s], rows):
                factor = F(truncation_denominator * q1 ** a * denominator_lcm ** (u - s),
                           q1 ** beta1 * q2 ** beta2)
                cleared = factor * value
                assert cleared.is_integral()
                cleared_vector.append(cleared)
                counts["cleared_V_entries"] += 1
            cleared_vectors.append(cleared_vector)
        original_cleared = {}
        for b in range(bmax + 1):
            vector = []
            for value, (_, _, beta1, beta2) in zip(original_columns[a, b], rows):
                factor = F(truncation_denominator * q1 ** a
                           * denominator_lcm ** (u + 1 + b), q1 ** beta1 * q2 ** beta2)
                cleared = factor * value
                assert cleared.is_integral()
                vector.append(cleared)
                counts["cleared_original_entries"] += 1
            original_cleared[b] = vector
        for b in range(bmax + 1):
            central_denominator = comb(2 * b, b)
            reduced = reduced_column(u, b, d, cleared_vectors)
            transformed = [sum((central_denominator * coefficient(u, j, b)
                                 * d ** (b - j) * original_cleared[j][index]
                                 for j in range(b + 1)), Gaussian())
                           for index in range(len(rows))]
            for actual, residual in zip(transformed, reduced):
                quotient = actual / d ** (2 * b + 1)
                assert actual.is_integral() and quotient.is_integral()
                assert quotient == central_denominator * residual
                counts["cleared_transformed_entries"] += 1

    assert counts["nonzero_tail_difference_entries"] > 0
    return {
        "status": "pass",
        "parameters": {"K": 2, "theta": "1/2", "w1": 2, "w2": 4,
                       "v0": 1, "F0": 5, "T1": 10, "T2": 20,
                       "remaining_budget": 32, "U": 7, "A": 1, "B": 1,
                       "compatible_shift_example": {"H": 160, "h": 128},
                       "p1": 1, "q1": 2, "p2": 2, "q2": 3,
                       "r1": r1, "r2": r2, "delta": delta, "L": 6, "d": d},
        "row_count": len(rows),
        "row_condition": "k + 4*beta1 + 8*beta2 < 32, j in {0,1}; k is s_time-h",
        "counts": counts,
        "truncation_denominator": "lcm(1..9)^(a+U+2) * lcm(1..19)^(U+2)",
        "truncation_denominator_digits_for_a0_a1": truncation_digits,
        "nonzero_G2_minus_G1_effects": tail_effects,
        "scope": "actual local coefficient row pack including unequal weights and unequal truncations; global interpolation separation hypotheses and full-row-rank completion not certified",
    }


def run_checks():
    counts = {"closed_coefficient_identities": 0,
              "canonical_elimination_entries": 0,
              "wedge_polynomials": 0,
              "wedge_determinants": 0,
              "central_denominator_checks": 0,
              "cleared_gaussian_entries": 0,
              "cleared_residual_determinants": 0,
              "tensor_determinants": 0}
    largest_actual_denominator = 0
    largest_actual_example = None
    largest_central_denominator = 0

    for u in range(25):
        for b in range(min(u, 8) + 1):
            central_denominator = comb(2 * b, b)
            numerators = [(-1) ** (b + a) * comb(b + a, a)
                          * comb(u + b + 1, b - a) for a in range(b + 1)]
            content = central_denominator
            actual_lcm = 1
            for a, numerator in enumerate(numerators):
                content = gcd(content, abs(numerator))
                actual_lcm = lcm(actual_lcm, coefficient(u, a, b).denominator)
            assert actual_lcm == central_denominator // content
            assert central_denominator % actual_lcm == 0
            assert central_denominator <= 4 ** b
            counts["central_denominator_checks"] += 1
            if actual_lcm > largest_actual_denominator:
                largest_actual_denominator = actual_lcm
                largest_actual_example = {"U": u, "b": b}
            largest_central_denominator = max(largest_central_denominator, central_denominator)
            for s in range(u + 1):
                assert transformed_direct(u, s, b) == transformed_closed(u, s, b)
                if s < b:
                    assert transformed_direct(u, s, b) == 0
                if s == b:
                    assert transformed_direct(u, s, b) != 0
                counts["closed_coefficient_identities"] += 1
        for bmax in range(min(u, 8) + 1):
            change = canonical_column_elimination(u, bmax)
            for a in range(bmax + 1):
                for b in range(bmax + 1):
                    expected = coefficient(u, a, b) if a <= b else F(0)
                    assert change[a][b] == expected
                    counts["canonical_elimination_entries"] += 1
            pivot_product = F(1)
            for b in range(bmax + 1):
                pivot_product *= transformed_closed(u, b, b)
            top_block = [row[:bmax + 1] for row in coefficient_matrix(u, bmax)[:bmax + 1]]
            assert determinant(top_block) == pivot_product
            counts["wedge_determinants"] += 1
            if bmax <= 4:
                wedge = polynomial_determinant(u, bmax)
                assert wedge == {(bmax + 1) ** 2: pivot_product}
                counts["wedge_polynomials"] += 1
            denominator_product = 1
            for b in range(bmax + 1):
                denominator_product *= comb(2 * b, b)
            for d in (Gaussian(1, 1), Gaussian(0, 2)):
                # The already-reduced block has entry 0 for s<b and
                # T[s,b]*d^(s-b) otherwise. Its column-b denominator is D_b.
                cleared_residual = [
                    [(comb(2 * b, b) * transformed_closed(u, s, b) * d ** (s - b))
                     if s >= b else Gaussian()
                     for b in range(bmax + 1)]
                    for s in range(bmax + 1)
                ]
                assert all(value.is_integral() for row in cleared_residual for value in row)
                residual_det = determinant(cleared_residual)
                quotient = residual_det / denominator_product
                assert quotient.is_integral()
                assert quotient == Gaussian(pivot_product)
                counts["cleared_residual_determinants"] += 1

    # The Gaussian integer numerator d = L*delta is the arithmetic analogue
    # of the analytic difference delta. Clearing the rational column-change
    # coefficients makes its additional powers divide the integer entries.
    gaussian_numerators = [Gaussian(1), Gaussian(-3), Gaussian(0, 2), Gaussian(1, 1)]
    for u in (1, 2, 4, 7, 10, 24):
        for b in range(min(u, 8) + 1):
            column_denominator = comb(2 * b, b)
            for d in gaussian_numerators:
                for scale in (2, 3, 5):
                    delta = d / scale
                    for s in range(u + 1):
                        # Cleared original column: L^(U+1+a) / L^(U-s)
                        # times the Taylor coefficient of (delta+v)^(U+1+a).
                        cleared = [
                            scale ** (1 + a + s)
                            * choose(u + 1 + a, u - s)
                            * delta ** (1 + a + s)
                            for a in range(b + 1)
                        ]
                        assert all(value.is_integral() for value in cleared)
                        direct = sum((column_denominator * coefficient(u, a, b)
                                      * d ** (b - a) * cleared[a]
                                      for a in range(b + 1)), Gaussian())
                        expected = (column_denominator * transformed_closed(u, s, b)
                                    * d ** (1 + b + s))
                        assert direct == expected
                        assert direct.is_integral()
                        quotient = direct / d ** (2 * b + 1)
                        assert quotient.is_integral()
                        if s < b:
                            assert not direct
                        counts["cleared_gaussian_entries"] += 1

    for u, amax, bmax in ((1, 0, 1), (3, 1, 1), (5, 1, 2), (7, 1, 1), (8, 2, 2)):
        base_det = determinant([row[:bmax + 1] for row in coefficient_matrix(u, bmax)[:bmax + 1]])
        exponent = (amax + 1) * (bmax + 1) ** 2
        for r, delta in ((Gaussian(F(2, 3)), Gaussian(F(-3, 5))),
                         (Gaussian(0, F(44, 7)), Gaussian(0, F(2, 791))),
                         (Gaussian(1, F(1, 2)), Gaussian(F(2, 3), F(1, 3)))):
            actual = determinant(tensor_block(u, amax, bmax, r, delta))
            expected = base_det ** (amax + 1) * delta ** exponent
            assert actual == expected
            counts["tensor_determinants"] += 1

    u, amax, bmax = 7, 1, 1
    coefficient_det = determinant([row[:2] for row in coefficient_matrix(u, 1)[:2]])
    special = {
        "U": u, "A": amax, "B": bmax,
        "binomial_block": [[8, 36], [28, 84]],
        "column_change": [[F(1), coefficient(u, 0, 1)], [F(0), F(1)]],
        "pivots": [transformed_closed(u, 0, 0), transformed_closed(u, 1, 1)],
        "block_coefficient_determinant": coefficient_det,
        "tensor_coefficient_determinant": coefficient_det ** 2,
        "delta_exponent": 8,
        "identity": "det = 112896 * delta^8",
    }
    assert coefficient_det == -336
    assert coefficient_det ** 2 == 112896
    return {
        "status": "pass",
        "scope": "finite exact algebra and arithmetic checks; not a proof of P2848 or the general interpolation hypotheses",
        "dependencies": "Python standard library only",
        "counts": counts,
        "parameter_range": {"U": [0, 24], "b": "0 <= b <= min(U,8)"},
        "denominators": {
            "guaranteed": "D_b = binomial(2b,b) <= 4^b",
            "exact_minimal": "D_b / gcd(D_b, all integer numerators in column b)",
            "largest_guaranteed_in_checked_range": largest_central_denominator,
            "largest_actual_in_checked_range": largest_actual_denominator,
            "first_largest_actual_example": largest_actual_example,
        },
        "special_example": special,
        "full_row_block": check_full_row_block(),
        "conclusions_checked": [
            "independent rational elimination matches the displayed closed coefficient formula",
            "transformed column b has first nonzero row s=b and factor delta^(2b+1)",
            "the small block determinant has exact difference order (B+1)^2",
            "the tensor block has exact difference order (A+1)*(B+1)^2",
            "after clearing rational change-of-column coefficients, the corresponding Gaussian integer power d^(2b+1) also divides every transformed column entry",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("results.json"))
    args = parser.parse_args()
    result = run_checks()
    data = json.dumps(result, default=encode, ensure_ascii=False, indent=2) + "\n"
    args.output.write_text(data, encoding="utf-8")
    print(data, end="")


if __name__ == "__main__":
    main()
