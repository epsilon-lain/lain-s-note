#!/usr/bin/env python3
"""Exact same-minor checks of row/column denominator allocation.

Finite algebra only. The examples do not certify the interpolation theorem,
large-H estimates, or a Diophantine property of pi.
"""
from itertools import product
from math import ceil, lcm
from pathlib import Path
import json
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / 'pi-construction-mechanisms-2026-11-01'
sys.path.insert(0, str(PREVIOUS))
from source_support import audit_source_files

t = sp.Symbol('t')


def matrix_entry(h, alpha, row, centers, truncations):
    j, s, beta = row
    if any(a < b for a, b in zip(alpha, beta)):
        return sp.Integer(0)
    expression = (1+t)**h
    for a, b, r, cutoff in zip(alpha, beta, centers, truncations):
        g = sum((sp.Rational((-1)**(k+1), k)*t**k
                 for k in range(1, cutoff)), sp.Integer(0))
        expression *= sp.binomial(a, b)*(j*r+g)**(a-b)
    return sp.Poly(sp.expand(expression), t).nth(s)


def gaussian_integer(z):
    z = sp.expand_complex(z)
    return bool(sp.re(z).is_Integer and sp.im(z).is_Integer)


def matching(support):
    """An independent support matching, not a determinant calculation."""
    occupied = {}
    def augment(row, seen):
        for col, exists in enumerate(support[row]):
            if not exists or col in seen:
                continue
            seen.add(col)
            if col not in occupied or augment(occupied[col], seen):
                occupied[col] = row
                return True
        return False
    assert all(augment(row, set()) for row in range(len(support)))
    return {row: col for col, row in occupied.items()}


def example(H, w0, v0, weights, theta, F0, numerators, denominators, K):
    m = len(weights)
    alpha_options = product(*(range(H//w+1) for w in weights))
    columns = [(h, a) for a in alpha_options
               for h in range(H//w0+1)
               if w0*h+sum(w*x for w, x in zip(weights, a)) <= H]
    columns.sort(key=lambda c: (sum(w*a for w, a in zip(weights, c[1])), c))
    beta_options = list(product(*(range(H//w+1) for w in weights)))
    rows = [(j, s, b) for j in range(K) for b in beta_options
            for s in range(H//v0+1)
            if v0*s+sum(w*x for w, x in zip(weights, b))/theta < H]
    T = [ceil(sp.Rational(F0*w, v0)) for w in weights]
    L = [lcm(*range(1, cutoff)) for cutoff in T]
    centers = [2*sp.I*sp.Rational(p, q) for p, q in zip(numerators, denominators)]
    full = sp.Matrix([[matrix_entry(h, a, row, centers, T)
                       for h, a in columns] for row in rows])
    pivots = full.rref()[1]
    assert len(pivots) == len(rows)
    chosen = [columns[i] for i in pivots]
    A = full[:, list(pivots)]
    determinant = sp.factor(A.det())
    assert determinant != 0
    D = [sum(a[i] for _, a in chosen)-sum(b[i] for _, _, b in rows)
         for i in range(m)]
    witness = matching([[A[r, c] != 0 for c in range(A.cols)]
                       for r in range(A.rows)])
    assert all(all(chosen[c][1][i] >= rows[r][2][i] for i in range(m))
               for r, c in witness.items())
    assert all(d >= 0 for d in D)
    col_factors = [sp.prod((q*l)**a for q, l, a in zip(denominators, L, alpha))
                   for _, alpha in chosen]
    row_factors = [sp.prod(sp.Rational(1, q*l)**b
                           for q, l, b in zip(denominators, L, beta))
                   for _, _, beta in rows]
    cleared = sp.diag(*row_factors)*A*sp.diag(*col_factors)
    assert all(gaussian_integer(x) for x in cleared)
    allocation = sp.prod((q*l)**d for q, l, d in zip(denominators, L, D))
    cleared_det = sp.factor(cleared.det())
    assert sp.expand(cleared_det-determinant*allocation) == 0
    assert gaussian_integer(cleared_det) and cleared_det != 0
    # Independently differentiate the original truncated expression.
    derivative_checks = 0
    for r, (j, s, beta) in enumerate(rows):
        for c, (h, alpha) in enumerate(chosen):
            if any(a < b for a, b in zip(alpha, beta)):
                continue
            expression = (1+t)**h
            for a, b, center, cutoff in zip(alpha, beta, centers, T):
                g = sum((sp.Rational((-1)**(k+1), k)*t**k
                         for k in range(1, cutoff)), sp.Integer(0))
                expression *= sp.binomial(a, b)*(j*center+g)**(a-b)
            actual = sp.diff(expression, t, s).subs(t, 0)/sp.factorial(s)
            assert sp.expand(actual-A[r, c]) == 0
            derivative_checks += 1
    source_DH = sp.prod(l**(H//w) for l, w in zip(L, weights))
    old_factor = source_DH**len(rows)*sp.prod(q**d for q, d in zip(denominators, D))
    assert old_factor % allocation == 0
    norm = len(rows)*H
    chi = sp.Rational(sum(sum(w*a for w, a in zip(weights, alpha))
                          for _, alpha in chosen), norm)
    bbar = sp.Rational(sum(sum(w*b for w, b in zip(weights, beta))
                          for _, _, beta in rows), norm)
    zeta = sp.Rational(w0*sum(h for h, _ in chosen), norm)
    assert chi-bbar == sp.Rational(sum(w*d for w, d in zip(weights, D)), norm)
    assert chi+zeta <= 1
    return dict(H=H, w0=w0, v0=v0, weights=weights, theta=str(theta),
                F0=F0, K=K, p=numerators, q=denominators, T=T, L=L,
                full_shape=full.shape, chosen_columns=chosen, full_rows=rows,
                determinant=str(determinant), nonnegative_degree_differences=D,
                support_matching=witness, all_scaled_entries_gaussian_integers=True,
                allocated_clearing_factor=str(allocation),
                allocated_cleared_determinant=str(cleared_det),
                source_clearing_factor=str(old_factor),
                factor_saving=str(old_factor//allocation),
                derivative_checks=derivative_checks,
                chi=str(chi), bbar=str(bbar), zeta=str(zeta),
                weighted_degree_cost=str(chi-bbar),
                finite_full_row_rank_verified=True,
                asymptotic_separation_or_pi_theorem_claim=False)


def main():
    examples = [example(12, 4, 4, (2,), sp.Rational(1, 3), 7, (22,), (7,), 2),
                example(8, 2, 4, (2,), sp.Rational(1, 3), 7, (22,), (7,), 2)]
    first = examples[0]
    assert first['nonnegative_degree_differences'] == [9]
    assert first['factor_saving'] == str(6**51)
    assert first['chi'] == '13/60' and first['bbar'] == '1/15'
    results = dict(status='passed', scope=__doc__,
                   source_audit=audit_source_files(['build/sections/determinant.tex']),
                   examples=examples)
    path = HERE/'arithmetic_rewrite_results.json'
    path.write_text(json.dumps(results, ensure_ascii=False, indent=2,
                               default=lambda x: int(x) if isinstance(x, sp.Integer) else str(x))+'\n')
    print(json.dumps(dict(status='passed', examples=len(examples), output=path.name)))


if __name__ == '__main__':
    main()
