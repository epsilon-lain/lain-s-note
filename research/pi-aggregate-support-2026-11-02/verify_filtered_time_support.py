"""Symbolic period-jet checks for positive time surplus (requires SymPy)."""
from functools import lru_cache
from itertools import product
from math import comb, factorial
from pathlib import Path
import json
import sympy as sp

Z, zeta, eta = sp.symbols('Z zeta eta')


@lru_cache(None)
def entry(row, column):
    j, t, a = row
    h, alpha = column
    if a > alpha:
        return sp.Integer(0)
    value = 0
    for ell in range(min(t, alpha-a)+1):
        n = t-ell
        numerator = sum((-1)**(h-k)*comb(h, k)*k**n for k in range(h+1))
        time = sp.Rational(numerator, factorial(n))
        value += comb(alpha, a)*comb(alpha-a, ell)*(j*Z)**(alpha-a-ell)*time
    return sp.expand(value)


def matrix(rows, columns):
    return sp.Matrix([[entry(row, column) for column in columns] for row in rows])


def profile(rows, columns):
    top = max([s for _, s, _ in rows]+[h for h, _ in columns])
    return {r: sum(s >= r for _, s, _ in rows)-sum(h >= r for h, _ in columns)
            for r in range(1, top+1)}


def six_row_example():
    rows = [(j, s, 0) for s in range(3) for j in range(2)]
    columns = [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0), (2, 1)]
    B = matrix(rows, columns)
    assert sp.factor(B.det()) == -2*Z**3
    cuts = profile(rows, columns)
    assert cuts == {1: 1, 2: 0}
    counts = {'total_targets': 0, 'aggregate_pass': 0,
              'cut_pass': 0, 'nonzero_period_minors': 0}
    for choices in product(*[range(s+1) for _, s, _ in rows]):
        target = [(j, s-a, a) for (j, s, _), a in zip(rows, choices)]
        gamma = sum(choices)
        cuts_pass = all(0 <= sum(s >= r for _, s, _ in rows)
                       -sum(t >= r for _, t, _ in target) <= capacity
                       for r, capacity in cuts.items())
        determinant = sp.factor(matrix(target, columns).det())
        counts['total_targets'] += 1
        counts['aggregate_pass'] += gamma <= sum(cuts.values())
        counts['cut_pass'] += cuts_pass
        counts['nonzero_period_minors'] += determinant != 0
        assert determinant == 0 or cuts_pass
    assert counts == {'total_targets': 36, 'aggregate_pass': 5,
                      'cut_pass': 3, 'nonzero_period_minors': 3}
    outside = [(j, s-a, a) for j, s, _ in rows for a in range(1, s+1)]
    C = matrix(outside, columns)
    R = sp.Matrix([[(-zeta)**a if j == jj and s == t+a else 0
                    for jj, t, a in outside] for j, s, _ in rows])
    X = C*B.inv()
    N = sp.simplify(R*X)
    marker = sp.factor((sp.eye(len(rows))+eta*N).det())
    assert marker == 1-eta*zeta
    assert N.rank() == 2
    for cut, capacity in cuts.items():
        ri = [i for i, (_, t, _) in enumerate(outside) if t < cut]
        ci = [i for i, (_, s, _) in enumerate(rows) if s >= cut]
        assert X.extract(ri, ci).rank() <= capacity
        nri = [i for i, (_, s, _) in enumerate(rows) if s <= cut]
        assert N.extract(nri, ci).rank() <= capacity
    return {'source_determinant': str(sp.factor(B.det())),
            'cut_capacities': cuts, 'target_counts': counts,
            'relative_boundary_rank': int(N.rank()), 'marker': str(marker)}


def sharp_degree_example():
    K, S = 2, 3
    rows = [(j, s, 0) for s in range(S+1) for j in range(K)]
    columns = [(0, alpha) for alpha in range(len(rows))]
    B = matrix(rows, columns)
    assert sp.factor(B.det()) == Z**16
    outside = [(j, s-a, a) for j, s, _ in rows for a in range(1, s+1)]
    C = matrix(outside, columns)
    # An independently specified exact row-reduction map, checked against C.
    X = sp.Matrix([[comb(t+a, a) if j == jj and s == t+a else 0
                    for jj, s, _ in rows] for j, t, a in outside])
    assert sp.simplify(C-X*B) == sp.zeros(*C.shape)
    R = sp.Matrix([[(-zeta)**a if j == jj and s == t+a else 0
                    for jj, t, a in outside] for j, s, _ in rows])
    N = sp.expand(R*X)
    desired = sp.diag(*[((1-zeta)**s-1) for _, s, _ in rows])
    assert sp.simplify(N-desired) == sp.zeros(len(rows))
    marker = sp.prod((1+eta*((1-zeta)**s-1))**K for s in range(1, S+1))
    cuts = profile(rows, columns)
    rho = sum(min(sum(s == r for _, s, _ in rows), capacity)
              for r, capacity in cuts.items())
    assert sum(cuts.values()) == 12 and rho == 6
    assert sp.Poly(marker, eta).degree() == rho
    assert sp.Poly(marker, zeta).degree() == 12
    return {'source_determinant': str(sp.factor(B.det())),
            'cut_capacities': cuts, 'total_surplus': sum(cuts.values()),
            'marker_degree_bound': rho, 'marker_degree': 6,
            'shear_degree': 12, 'marker': str(sp.factor(marker))}


def raw_filtration_check():
    groups = [[1, 3, 7], [2, 6]]
    determinants = []
    for powers in groups:
        jet = sp.Matrix([[comb(h, i) for h in powers] for i in range(len(powers))])
        assert jet.det() != 0
        determinants.append(int(jet.det()))
        for r in range(len(powers)+1):
            assert jet[:r, :].rank() == min(r, len(powers))
    assert determinants == [24, 4]
    return {'raw_exponent_groups': groups, 'first_jet_determinants': determinants,
            'adapted_time_orders': [[0, 1, 2], [0, 1]]}


def main():
    result = {'status': 'passed', 'six_row_example': six_row_example(),
              'sharp_degree_example': sharp_degree_example(),
              'raw_filtration_check': raw_filtration_check(),
              'scope': 'Symbolic period-jet support and rank checks; no '
                       'quantitative determinant or Diophantine conclusion.'}
    Path(__file__).with_name('filtered_time_verification.json').write_text(
        json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
