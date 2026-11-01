#!/usr/bin/env python3
"""Exact nontrivial full-packet tail elimination, with two rational centers."""
from math import comb
from pathlib import Path
import json
import sympy as sp


def truncated_product(a, b, degree):
    return [sp.expand(sum(a[j] * b[k-j] for j in range(k+1)))
            for k in range(degree+1)]


def powers(series, count, degree):
    result = [[sp.Integer(1)] + [sp.Integer(0)] * degree]
    for _ in range(count):
        result.append(truncated_product(result[-1], series, degree))
    return result


def matrix_for(series, rows, columns, centers, degree):
    pp = powers(series, degree, degree)
    def entry(j, s, beta, h, alpha):
        if alpha < beta:
            return sp.Integer(0)
        n = alpha-beta
        return sp.expand(comb(alpha, beta) * sum(
            comb(n, k) * centers[j]**(n-k) * sum(
                comb(h, ell) * pp[k][s-ell]
                for ell in range(min(h,s)+1))
            for k in range(min(n,s)+1)))
    return sp.Matrix([[entry(j,s,beta,h,alpha) for h,alpha in columns]
                      for j,s,beta in rows])


def main():
    H, v0, v1, T = 10, 2, 4, 3
    centers = [sp.Integer(0), 44*sp.I/7]
    columns = [(h,a) for h in range(H+1) for a in range(H-h+1)]
    rows = [(j,s,b) for j in range(2) for s in range(H//v0+1)
            for b in range(H//v1+1) if v0*s+v1*b < H]
    degree = max(s for j,s,b in rows)
    log_series = [sp.Integer(0)] + [sp.Rational((-1)**(k+1),k)
                                                for k in range(1,degree+1)]
    G_series = [x if k<T else sp.Integer(0) for k,x in enumerate(log_series)]
    tau = [sp.expand(g-l) for g,l in zip(G_series,log_series)]
    tau_powers = powers(tau, H//v1, degree)
    log_matrix = matrix_for(log_series, rows, columns, centers, degree)
    G_matrix = matrix_for(G_series, rows, columns, centers, degree)
    def coefficient(output, input_row):
        j,s,b = output
        j2,s2,gamma = input_row
        if j != j2 or s2>s or gamma<b:
            return sp.Integer(0)
        return comb(gamma,b)*tau_powers[gamma-b][s-s2]
    U = sp.Matrix([[coefficient(out,inp) for inp in rows] for out in rows])
    assert U*log_matrix == G_matrix
    order = sorted(range(len(rows)), key=lambda r:(v0*rows[r][1]+v1*rows[r][2],rows[r]))
    ordered_U = U.extract(order,order)
    assert ordered_U.is_lower
    assert all(ordered_U[i,i] == 1 for i in range(len(rows)))
    assert U.det() == 1
    nonidentity = [(i,j,str(U[i,j])) for i in range(len(rows))
                   for j in range(len(rows)) if i!=j and U[i,j]!=0]
    assert nonidentity
    for i,j,_ in nonidentity:
        assert v0*rows[j][1]+v1*rows[j][2] < v0*rows[i][1]+v1*rows[i][2]
    pivots = G_matrix.rref()[1]
    assert len(pivots) == len(rows)
    det_G = sp.expand(G_matrix[:,pivots].det(method="domain-ge"))
    det_log = sp.expand(log_matrix[:,pivots].det(method="domain-ge"))
    assert det_G != 0
    assert det_G == det_log
    single_row = rows.index((0,3,0))
    y3, X = columns.index((3,0)), columns.index((0,1))
    partial_G = G_matrix[single_row,y3] + G_matrix[single_row,X]
    partial_log = log_matrix[single_row,y3] + log_matrix[single_row,X]
    assert partial_G == 1 and partial_log == sp.Rational(4,3)
    result = {
        "parameters": {"H":H,"w0":1,"w1":1,"theta":"1/4",
                       "v0":v0,"v1":v1,"T":T,"F0":5,
                       "centers":[str(x) for x in centers]},
        "rows":len(rows),"columns":len(columns),
        "tail_series_through_maximum_row_t_order":[str(x) for x in tau],
        "nonidentity_U_entries":[{"output":rows[i],"input":rows[j],"value":v}
                                  for i,j,v in nonidentity],
        "U_determinant":"1","U_is_unit_lower_triangular_after_weight_ordering":True,
        "direct_matrix_identity_A_G_equals_U_A_log":True,
        "selected_column_indices":list(pivots),
        "selected_columns":[columns[i] for i in pivots],
        "same_nonzero_minor_G":str(det_G),"same_nonzero_minor_log":str(det_log),
        "full_row_rank":len(pivots),
        "ablation_missing_other_rows":{"P":"Y^3+X","row":[0,3,0],
                                        "G_value":str(partial_G),
                                        "log_value":str(partial_log),
                                        "equal":False},
        "scope":"Exact algebraic packet identity; not a proof of asymptotic interpolation or pi approximation",
    }
    target = Path(__file__).with_name("tail_elimination_results.json")
    target.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"output":str(target),"rows":len(rows),"columns":len(columns),
                      "nonidentity_U_entries":len(nonidentity),"all_checks_passed":True},indent=2))


if __name__ == "__main__":
    main()
