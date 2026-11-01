"""Exact finite laboratory for the pinned pi017 determinant construction.

No numerical rank decisions and no replacement of the chosen minor occur.
Z denotes the exact period 2*pi*i in polynomial formulas after using
exp(h*j*Z)=1 for integer h,j. The general rule that every nonzero Q[Z]
polynomial remains nonzero at that period uses its transcendence; this
script does not prove transcendence. Its displayed monomial determinants
need only Z != 0.

Run: python research/pi-construction-mechanisms-2026-11-01/period_lab.py
"""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import permutations, product
import hashlib
import json
import math
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
Z, r1, r2, e1, e2, t = sp.symbols("Z r1 r2 e1 e2 t")


def stringify(value):
    if isinstance(value, dict):
        return {str(k): stringify(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [stringify(v) for v in value]
    if isinstance(value, sp.Basic):
        return str(value)
    return value


def truncated_log(order):
    return sum((sp.Rational((-1) ** (k + 1), k) * t**k
                for k in range(1, order + 1)), sp.Integer(0))


def coefficient(expr, order):
    return sp.Poly(sp.expand(expr), t).nth(order)


@lru_cache(None)
def exact_period_entry(h, alpha, j, ell, a):
    if any(ai > alphai for ai, alphai in zip(a, alpha)):
        return sp.Integer(0)
    factor = sp.prod(sp.binomial(alphai, ai)
                     for ai, alphai in zip(a, alpha))
    degree = sum(alpha) - sum(a)
    return sp.expand(factor * coefficient(
        (1+t)**h * (j*Z + truncated_log(ell))**degree, ell))


def period_matrix(columns, rows):
    return sp.Matrix([[exact_period_entry(h, tuple(alpha), j, ell, tuple(a))
                       for h, alpha in columns] for j, ell, a in rows])


def original_matrix(columns, rows, centers=(r1, r2), truncations=(3, 3)):
    output = []
    for j, s, beta in rows:
        row = []
        for h, alpha in columns:
            if any(bi > ai for bi, ai in zip(beta, alpha)):
                row.append(sp.Integer(0))
                continue
            expression = (1+t)**h
            for center, cutoff, ai, bi in zip(centers, truncations, alpha, beta):
                expression *= sp.binomial(ai, bi) * (
                    j*center + truncated_log(cutoff-1))**(ai-bi)
            row.append(coefficient(expression, s))
        output.append(row)
    return sp.Matrix(output)


def leibniz_determinant(matrix):
    """Independent determinant, without elimination/Bareiss or Matrix.det."""
    n = matrix.rows
    total = sp.Integer(0)
    for p in permutations(range(n)):
        parity = sum(p[i] > p[j] for i in range(n) for j in range(i+1, n))
        total += (-1)**parity * sp.prod(matrix[i, p[i]] for i in range(n))
    return sp.expand(total)


def check_equal(left, right):
    assert sp.expand(left-right) == 0, (left, right)


def source_integrity():
    from source_support import audit_source_files
    return audit_source_files()


def six_by_six():
    columns = [(0, (0, 0)), (1, (0, 0)), (0, (0, 1)),
               (0, (0, 2)), (0, (0, 3)), (0, (1, 0))]
    original_rows = [(j, s, (0, 0)) for j in range(2) for s in range(3)]
    actual = original_matrix(columns, original_rows)
    determinant = sp.factor(actual.det(method="domain-ge"))
    check_equal(determinant, -3*r2**2*(r1-r2))
    check_equal(leibniz_determinant(actual), determinant)
    # Derive the actual jets independently by differentiation of the true log.
    actual_checks = 0
    for row_index, (j, s, _) in enumerate(original_rows):
        for col_index, (h, alpha) in enumerate(columns):
            expression = (1+t)**h * (j*r1+sp.log(1+t))**alpha[0] * (
                j*r2+sp.log(1+t))**alpha[1]
            derivative = sp.diff(expression, t, s).subs(t, 0)/sp.factorial(s)
            check_equal(derivative, actual[row_index, col_index])
            actual_checks += 1

    pure_rows = [(j, s, (0, 0)) for j, s, _ in original_rows]
    pure = period_matrix(columns, pure_rows)
    assert pure.rank() == 5 and pure.det() == 0
    assert actual.rank() == 6
    collapsed = actual.copy()
    collapsed[:, 1] = sp.Matrix([int(s == 0) for _, s, _ in original_rows])
    assert collapsed.rank() == 5 and collapsed.det() == 0
    collapsed_pure = collapsed.subs({r1: Z, r2: Z})
    assert collapsed_pure.rank() == 4

    # Keep the SAME six columns and enumerate all row-translation terms.
    # All tails vanish through order 2, and epsilon=0 at j=0.
    choices = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0)]
    nonzero_terms, zero_terms = [], 0
    reconstructed = sp.Integer(0)
    grouped_by_cost = defaultdict(lambda: sp.Integer(0))
    vector_determinants = {}
    for selection in product(choices, repeat=3):
        row_indices = [(0, s, (0, 0)) for s in range(3)] + [
            (1, s, a) for s, a in enumerate(selection)]
        matrix = period_matrix(columns, row_indices)
        det = sp.factor(matrix.det(method="domain-ge"))
        if det == 0:
            zero_terms += 1
            continue
        scalar = sp.prod(e1**a[0]*e2**a[1] for a in selection)
        cost = sum(sum(a) for a in selection)
        contribution = sp.expand(scalar*det)
        reconstructed += contribution
        grouped_by_cost[cost] += contribution
        counts = Counter([(0, 0)]*3+list(selection))
        nonzero_terms.append(dict(selection=selection, vector_determinant=det,
                                  scalar=scalar, total_transverse_degree=cost,
                                  group_counts=dict(counts),
                                  sum_squared_group_counts=sum(n*n for n in counts.values())))
        vector_determinants[selection] = det
    reconstructed = sp.factor(reconstructed)
    check_equal(reconstructed, determinant.subs({r1: Z+e1, r2: Z+e2}))
    assert min(term["total_transverse_degree"] for term in nonzero_terms) == 1
    # A nonzero individual translated determinant does not prevent cancellation.
    assert any(sp.expand(term["scalar"].subs(e1, e2)*term["vector_determinant"]) != 0
               for term in nonzero_terms)
    check_equal(reconstructed.subs(e1, e2), 0)
    chosen = ((1, 0), (0, 0), (0, 0))
    assert chosen in vector_determinants
    check_equal(vector_determinants[chosen], -3*Z**2)

    # This packet really is the source's entire row simplex for these finite data.
    H, w0, v0, weights, theta, F0 = 15, 5, 5, (2, 5), sp.Rational(1, 10), 21
    enumerated_rows = [(j, s, (b1, b2)) for j in range(2)
                       for s in range(4) for b1 in range(2) for b2 in range(2)
                       if v0*s+(weights[0]*b1+weights[1]*b2)/theta < H]
    assert enumerated_rows == original_rows
    assert all(w0*h+sum(w*a for w, a in zip(weights, alpha)) <= H
               for h, alpha in columns)
    T = [math.ceil(sp.Rational(F0*w, v0)) for w in weights]
    assert T == [9, 21]
    arithmetic_actual = original_matrix(columns, original_rows,
                                         (sp.Rational(44, 7)*sp.I, sp.Rational(710, 113)*sp.I),
                                         tuple(T))
    # Independent use of only the relevant Taylor coefficients gives same matrix.
    assert arithmetic_actual == actual.subs({r1: sp.Rational(44, 7)*sp.I,
                                             r2: sp.Rational(710, 113)*sp.I})
    arithmetic_det = sp.factor(arithmetic_actual.det())
    assert arithmetic_det != 0
    q = (7, 113)
    multipliers = [sp.prod(qi**ai for qi, ai in zip(q, alpha)) for _, alpha in columns]
    scaled = arithmetic_actual * sp.diag(*multipliers)
    cleared = 2*scaled  # A sharper finite clearing than the source's universal D_H.
    assert all(sp.re(x).is_Integer and sp.im(x).is_Integer for x in cleared)
    cleared_det = sp.factor(cleared.det())
    assert sp.re(cleared_det).is_Integer and sp.im(cleared_det).is_Integer
    assert cleared_det != 0
    check_equal(cleared_det, 2**6 * sp.prod(multipliers) * arithmetic_det)
    L = [sp.ilcm(*range(1, cutoff)) for cutoff in T]
    D = sp.prod(Li**(H//wi) for Li, wi in zip(L, weights))
    assert D % 2 == 0
    source_cleared = D*scaled
    assert all(sp.re(x).is_Integer and sp.im(x).is_Integer for x in source_cleared)

    return dict(columns=columns, rows=original_rows, matrix=actual.tolist(),
                determinant=determinant, independent_leibniz_check=True,
                independent_actual_derivative_checks=actual_checks,
                original_rank_over_Q_r1_r2=6,
                pure_period_matrix=pure.tolist(), pure_period_rank_over_Q_Z=5,
                exponential_collapsed_original_rank=5,
                exponential_collapsed_pure_rank=4,
                enumerated_translated_choices=125, zero_vector_determinants=zero_terms,
                nonzero_vector_determinants=len(nonzero_terms),
                minimum_surviving_transverse_degree=1,
                grouped_contributions={cost: sp.factor(v) for cost, v in sorted(grouped_by_cost.items())},
                reconstructed_determinant=reconstructed,
                equal_error_specialization_cancels=True,
                chosen_shifted_determinant=vector_determinants[chosen],
                nonzero_terms=nonzero_terms,
                arithmetic_example=dict(H=H, w0=w0, v0=v0, weights=weights, theta=theta,
                                        F0=F0, T=T, q=q, p=(22, 355),
                                        determinant=arithmetic_det,
                                        sharpened_entry_clearing_factor=2,
                                        column_multipliers=multipliers,
                                        sharpened_cleared_determinant=cleared_det,
                                        source_L=L, source_D_H=D,
                                        entire_row_simplex_used=True,
                                        selected_columns_admissible=True,
                                        interpolation_separation_threshold_claim=False))


def general_translation_checks():
    """Exercise beta binomials, genuine tails, and different truncation orders."""
    checks = 0
    nonzero_tail_summands = 0
    columns = [(0, (1, 2)), (1, (2, 1)), (2, (2, 2)), (0, (0, 3))]
    beta_options = [(0, 0), (1, 0), (0, 1), (1, 1)]
    T = (2, 3)
    log = truncated_log(4)
    tails = [truncated_log(cutoff-1)-log for cutoff in T]
    for h, alpha in columns:
        for beta in beta_options:
            for s in range(5):
                if any(b > aa for b, aa in zip(beta, alpha)):
                    continue
                lhs = original_matrix([(h, alpha)], [(1, s, beta)], truncations=T)[0, 0]
                rhs = sp.Integer(0)
                for a in product(*(range(b, aa+1) for b, aa in zip(beta, alpha))):
                    initial_binomial = sp.prod(sp.binomial(aa, b) for aa, b in zip(a, beta))
                    difference = tuple(aa-b for aa, b in zip(a, beta))
                    for d in product(*(range(n+1) for n in difference)):
                        scalar = initial_binomial * sp.prod(
                            sp.binomial(n, di)*(ri-Z)**(n-di)
                            for n, di, ri in zip(difference, d, (r1, r2)))
                        tail = sp.prod(tau**di for tau, di in zip(tails, d))
                        for k in range(s+1):
                            tail_coefficient = coefficient(tail, k)
                            value = exact_period_entry(h, alpha, 1, s-k, a)
                            summand = scalar*tail_coefficient*value
                            rhs += summand
                            if any(d) and summand != 0:
                                nonzero_tail_summands += 1
                check_equal(lhs, rhs)
                checks += 1
    assert nonzero_tail_summands > 0
    # The smallest omitted log term has the correct positive sign when T=2.
    assert coefficient(tails[0], 2) == sp.Rational(1, 2)
    assert coefficient(tails[1], 3) == -sp.Rational(1, 3)
    derivative_checks = 0
    for h, alpha in columns:
        for a in product(*(range(ai+1) for ai in alpha)):
            for j in range(3):
                for ell in range(4):
                    factor = sp.prod(sp.binomial(ai, aa) for ai, aa in zip(alpha, a))
                    expr = factor*(1+t)**h*(j*Z+sp.log(1+t))**(sum(alpha)-sum(a))
                    expected = sp.diff(expr, t, ell).subs(t, 0)/sp.factorial(ell)
                    check_equal(expected, exact_period_entry(h, alpha, j, ell, a))
                    derivative_checks += 1
    return dict(full_row_identity_checks=checks,
                nonzero_tail_summands_exercised=nonzero_tail_summands,
                independent_period_derivative_checks=derivative_checks,
                truncations=T, tail_first_coefficients=[sp.Rational(1, 2), -sp.Rational(1, 3)])


def collisions():
    columns = [(0, (0,)), (0, (1,)), (0, (2,)), (0, (3,)), (1, (0,))]
    rows = [(0, 0, (0,)), (0, 1, (0,)), (1, 0, (0,)),
            (1, 1, (0,)), (2, 0, (0,))]
    exact = period_matrix(columns, rows)
    check_equal(exact.det(), -6*Z**5)
    check_equal(leibniz_determinant(exact), -6*Z**5)
    assert exact.rank() == 5

    def c(d):
        return [int(d == 0), int(d == 1), int(d == 2), int(d == 3),
                sp.Rational(1, sp.factorial(d))]

    repeated = []
    for degrees in [(0, 1, 2, 3, 3), (0, 0, 1, 2, 4),
                    (0, 1, 2, 4, 4), (1, 2, 3, 4, 1)]:
        matrix = sp.Matrix([c(d) for d in degrees])
        assert matrix.det() == 0
        repeated.append(dict(degrees=degrees, determinant=0))
    minimal = sp.Matrix([c(d) for d in range(5)])
    check_equal(minimal.det(), sp.Rational(1, 24))
    moments = sp.Matrix([[coefficient((j*Z+truncated_log(ell))**d, ell)
                         for d in range(5)] for j, ell, _ in rows])
    check_equal(moments.det(), 4*Z**8)

    # Off-period centers have e^(j*z)=E^j, then E=1 is the exact-period condition.
    x, E = sp.symbols("x E")
    off_period = sp.Matrix([
        [(x*j)**d if ell == 0 else (d*(x*j)**(d-1) if d else 0)
         for d in range(4)] + [E**j]
        for j, ell, _ in rows])
    off_det = sp.factor(off_period.det())
    check_equal(off_det.subs({E: 1, x: Z}), -6*Z**5)
    off_series = sp.series(off_det.subs(E, sp.exp(x)), x, 0, 10)
    assert sp.expand(off_series.removeO()).coeff(x, 8) == sp.Rational(1, 6)
    assert all(sp.expand(off_series.removeO()).coeff(x, degree) == 0 for degree in range(8))

    # Two transverse groups permit repeated degree across groups.
    mixed_columns = [(0, (1, 0)), (0, (0, 1)), (1, (1, 0)), (1, (0, 1))]

    def c_mixed(a, d):
        z = sp.symbols("z")
        output = []
        for h, alpha in mixed_columns:
            if any(ai > alphai for ai, alphai in zip(a, alpha)):
                output.append(sp.Integer(0))
                continue
            factor = sp.prod(sp.binomial(alphai, ai) for alphai, ai in zip(alpha, a))
            f = factor*sp.exp(h*z)*z**(sum(alpha)-sum(a))
            output.append(sp.diff(f, z, d).subs(z, 0)/sp.factorial(d))
        return output

    mixed_labels = [((1, 0), 0), ((0, 1), 0), ((1, 0), 1), ((0, 1), 1)]
    mixed = sp.Matrix([c_mixed(a, d) for a, d in mixed_labels])
    assert mixed.det() == 1
    repeated_same_pair = mixed.copy()
    repeated_same_pair[1, :] = repeated_same_pair[0, :]
    assert repeated_same_pair.det() == 0
    return dict(single_group=dict(columns=columns, rows=rows, exact_period_matrix=exact.tolist(),
                                  exact_period_determinant=-6*Z**5, rank_over_Q_Z=5,
                                  independent_leibniz_check=True,
                                  repeated_degree_checks=repeated,
                                  minimum_distinct_degrees=list(range(5)),
                                  minimum_degree_sum=10, coefficient_determinant=sp.Rational(1, 24),
                                  moment_matrix=moments.tolist(), moment_determinant=4*Z**8,
                                  first_Taylor_selection_contribution=Z**8/6,
                                  off_period_matrix=off_period.tolist(),
                                  off_period_determinant=off_det,
                                  off_period_small_x_series=str(off_series)),
                transverse_groups=dict(columns=mixed_columns, labels=mixed_labels,
                                       coefficient_matrix=mixed.tolist(), determinant=1,
                                       degree_sum=2, group_sizes=[2, 2],
                                       groupwise_minimum_degree_sum=2,
                                       incorrectly_global_distinct_minimum=6,
                                       repeated_same_pair_determinant=0))


def main():
    results = dict(status="passed", sympy_version=sp.__version__, source=source_integrity(),
                   six_by_six=six_by_six(), translation_identity=general_translation_checks(),
                   collisions=collisions(),
                   scope="Exact finite identities and conditional mechanism audit; no asymptotic interpolation or pi endpoint theorem is certified.")
    path = HERE / "results_period_lab.json"
    path.write_text(json.dumps(stringify(results), indent=2)+"\n")
    print(json.dumps(dict(status="passed", output=str(path),
                          six_by_six_determinant=str(results["six_by_six"]["determinant"]),
                          translated_nonzero_terms=results["six_by_six"]["nonzero_vector_determinants"],
                          translation_checks=results["translation_identity"]), default=str, indent=2))


if __name__ == "__main__":
    main()
