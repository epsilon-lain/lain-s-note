"""Symbolic rank oracle for finite pure-period row vectors.

Requires SymPy. Z is a formal variable; no numerical pi value is used.
Finite computations do not establish an asymptotic displacement bound.
"""
from functools import lru_cache
from itertools import product
import json
import sympy as sp
from rank_minor import cheapest_minor

Z, t = sp.symbols('Z t')


@lru_cache(None)
def entry(h, alpha, j, ell, a):
    if len(alpha) != len(a) or any(x > y for x, y in zip(a, alpha)):
        return sp.Integer(0)
    factor = sp.prod(sp.binomial(y, x) for x, y in zip(a, alpha))
    degree = sum(alpha) - sum(a)
    logt = sum(sp.Rational((-1) ** (k + 1), k) * t ** k for k in range(1, ell + 1))
    expr = sp.Poly(sp.expand((1 + t) ** h * (j * Z + logt) ** degree), t)
    return sp.expand(factor * expr.nth(ell))


def row(columns, j, ell, a):
    return [entry(h, tuple(alpha), j, ell, tuple(a)) for h, alpha in columns]


def matrix(columns, rows):
    return sp.Matrix([row(columns, j, ell, a) for j, ell, a in rows])


def main():
    # The existing four-column example, now used as an exact rank-oracle test.
    columns = [(0, (0, 0)), (0, (0, 1)), (0, (0, 2)), (0, (1, 0))]
    original_rows = [(0, 0), (0, 1), (1, 0), (1, 1)]
    no_shift = matrix(columns, [(j, s, (0, 0)) for j, s in original_rows])
    assert no_shift.rank() == 3 and no_shift.det() == 0
    r1, r2 = sp.symbols('r1 r2')
    candidates = [(0, 0), (1, 0), (0, 1), (0, 2)]
    choices = [[(0, 0)], [(0, 0)], candidates, candidates]
    terms, reconstructed = [], sp.Integer(0)
    zero_vector_determinants = 0
    for selection in product(*choices):
        mat = matrix(columns, [(j, s, a) for (j, s), a in zip(original_rows, selection)])
        det = sp.expand(mat.det())
        if not det:
            zero_vector_determinants += 1
            continue
        scalar = sp.prod(((r1 - Z) ** a[0] * (r2 - Z) ** a[1])
                         for (j, _), a in zip(original_rows, selection) if j == 1)
        reconstructed += scalar * det
        terms.append(dict(selection=selection, vector_determinant=str(det),
                          scalar=str(scalar), degree_cost=sum(sum(a) for a in selection)))
    reconstructed = sp.factor(reconstructed)
    assert sp.expand(reconstructed - 2 * r2 * (r2 - r1)) == 0
    minimum_degree = min(term['degree_cost'] for term in terms)
    assert minimum_degree == 1

    # Same row packet, now all original columns of total degree <= 2.
    # Greedy cheap selection can preserve the pure-period term: cheap selection
    # and cancellation must be certified on the SAME columns.
    all_columns = [(h, (a1, a2)) for h in range(3) for a1 in range(3 - h)
                   for a2 in range(3 - h - a1)]
    actual = []
    for j, s in original_rows:
        actual.append([sp.expand((1+t)**h * (2*j+t)**alpha[0] * (3*j+t)**alpha[1]).coeff(t, s)
                       for h, alpha in all_columns])
    costs = [sum(alpha) for h, alpha in all_columns]
    cheapest = cheapest_minor(actual, costs)
    chosen_columns = [all_columns[i] for i in cheapest['chosen']]
    cheap_period = matrix(chosen_columns, [(j, s, (0, 0)) for j, s in original_rows])
    assert cheapest['minimum_cost'] == 2 and cheap_period.rank() == 4

    # Compare symbolic formula against direct z-derivatives at formal periods:
    # e^(h j omega)=1 is used algebraically, not numerically approximated.
    checks = 0
    for h, alpha in [(0, (2, 1)), (1, (1, 2)), (2, (2, 0))]:
        for j, ell, a in product(range(2), range(3), [(0, 0), (1, 0), (0, 1)]):
            if any(x > y for x, y in zip(a, alpha)):
                continue
            # Independent coefficient extraction by differentiation of log(1+t).
            factor = sp.prod(sp.binomial(y, x) for x, y in zip(a, alpha))
            f = factor * (1 + t) ** h * (j * Z + sp.log(1 + t)) ** (sum(alpha) - sum(a))
            direct = sp.diff(f, t, ell).subs(t, 0) / sp.factorial(ell)
            assert sp.expand(direct - entry(h, alpha, j, ell, a)) == 0
            checks += 1
    print(json.dumps(dict(status='passed', sympy_version=sp.__version__,
                         no_shift_rank=3, row_count=4,
                         enumerated_selections=16,
                         zero_vector_determinants=zero_vector_determinants,
                         minimum_surviving_degree=minimum_degree,
                         reconstructed_determinant=str(reconstructed),
                         actual_full_packet_columns=len(all_columns),
                         cheapest_columns=chosen_columns,
                         cheapest_total_horizontal_cost=str(cheapest['minimum_cost']),
                         cheapest_pure_period_rank=cheap_period.rank(),
                         cheapest_pure_period_determinant=str(sp.factor(cheap_period.det())),
                         independent_derivative_checks=checks, nonzero_terms=terms,
                         scope='Exact finite polynomial identities; no pi asymptotic claim.'), indent=2))


if __name__ == '__main__':
    main()
