#!/usr/bin/env python3
"""Exact coherent finite-section controls, not a pi approximation theorem."""
from fractions import Fraction as F
from pathlib import Path
import json


def determinant(matrix):
    a=[row[:] for row in matrix]
    total=F(1)
    for j in range(len(a)):
        pivot=next((i for i in range(j,len(a)) if a[i][j]),None)
        if pivot is None:return F(0)
        if pivot!=j:
            a[j],a[pivot]=a[pivot],a[j]
            total=-total
        value=a[j][j]
        total*=value
        for i in range(j+1,len(a)):
            c=a[i][j]/value
            a[i]=[x-c*y for x,y in zip(a[i],a[j])]
    return total


def projection(m):
    return [[F(int(i==j))-F(3,2**(i+j)) for j in range(1,m+1)]
            for i in range(1,m+1)]


def main():
    cases=[]
    previous=None
    for m in range(1,11):
        A=projection(m)
        if previous is not None:
            assert [row[:m-1] for row in A[:m-1]]==previous
        previous=A
        d=determinant(A)
        assert d==F(1,4**m)
        u=[F(1,2**j) for j in range(1,m+1)]
        image=[sum(x*y for x,y in zip(row,u)) for row in A]
        assert image==[F(1,4**m)*x for x in u]
        partial_norm_squared=3*sum(x*x for x in u)
        assert partial_norm_squared==1-F(1,4**m)
        B=[[F(2**(i+j))*A[i-1][j-1] for j in range(1,m+1)]
           for i in range(1,m+1)]
        assert all(x.denominator==1 for row in B for x in row)
        db=determinant(B)
        assert db==4**(m*(m-1)//2)
        assert d.denominator*d==1
        cases.append(dict(m=m,finite_determinant=str(d),
                          smallest_eigenvalue=str(F(1,4**m)),
                          condition_number=1 if m==1 else 4**m,
                          determinant_denominator=d.denominator,
                          cleared_determinant='1',
                          entrywise_integer_matrix_determinant=str(db),
                          exact_eigenvector_verified=True,
                          previous_leading_block_preserved=True,
                          operator_norm_tail_bound=str(F(2,2**m))))
    # A nonzero finite margin would already contradict the two inequalities.
    margin_controls=[]
    for k in range(1,31):
        A,B,G,R=F(k),F(2*k),F(3*k)+F(1,k),F(1,2*k)
        margin=G-A-B-R
        assert margin==F(1,2*k)>0
        margin_controls.append(dict(k=k,strict_finite_margin=str(margin)))
    # Purely algebraic check of the fixed-prefix separation bound.
    prefix=[]
    for theta in (F(1,4),F(1,2),F(3,4)):
        ratio=F(100)
        for m in (1,5,20,100):
            required_lower=2*(m+2)*theta/(1-theta)
            prefix.append(dict(theta=str(theta),m=m,w1_over_v0=str(ratio),
                               strict_lower_required=str(required_lower),
                               this_necessary_condition_satisfied=ratio>required_lower))
    output=Path(__file__).with_name('limit_audit_results.json')
    output.write_text(json.dumps(dict(status='passed',finite_sections=cases,
                       finite_margin_arithmetic_checks=margin_controls,
                       fixed_prefix_bounds=prefix,
                       scope='Exact finite controls for the stated general proofs; no pi tail conclusion.'),indent=2)+'\n')
    print(json.dumps(dict(status='passed',coherent_sections=len(cases),
                          determinant_eliminations=2*len(cases),
                          finite_margin_controls=len(margin_controls))))


if __name__=='__main__':main()
