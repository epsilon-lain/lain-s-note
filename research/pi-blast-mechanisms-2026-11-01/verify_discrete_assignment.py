"""Exact second-moment check of the row-constraint construction.

Uses wi=v0=1, theta=1/2, T=1 and K=2. No random samples: each
row's output distribution and its contribution to the variance are exact.
The infinite statement is proved in ROW_CONSTRAINT_ENDPOINT.md.
"""
from fractions import Fraction as F
from itertools import product
from collections import defaultdict
from pathlib import Path
import json


def check(m, H, repair_zero_center=False):
    K, D = 2, 2**m
    mean_scaled = defaultdict(int)
    linear_sum = F()
    random_rows = total_rows = 0
    for beta in product(range((H-1)//2+1), repeat=m):
        weight = sum(beta)
        for s in range(max(0, H-2*weight)):
            for j in range(K):
                total_rows += 1
                interior = 2*weight+s+2*m+1 <= H
                if interior and not (repair_zero_center and j == 0):
                    random_rows += 1
                    outputs = [tuple(2*b+bit for b,bit in zip(beta,bits))
                               for bits in product((0,1),repeat=m)]
                    assert len(set(outputs)) == D
                    for a in outputs:
                        assert all(x >= y for x,y in zip(a,beta))
                        assert sum(a) <= H
                        mean_scaled[a] += 1
                        linear_sum += F(2*(sum(a)-weight),D)
                else:
                    mean_scaled[beta] += D
    assert sum(mean_scaled.values()) == D*total_rows
    expected_square_sum = (F(sum(v*v for v in mean_scaled.values()),D*D)
                           + F(random_rows*(D-1),D))
    normalized_linear = linear_sum/(total_rows*H)
    normalized_squares = expected_square_sum/(total_rows*H)
    return {
        'm':m,'H':H,'zero_center_repaired':repair_zero_center,
        'rows':total_rows,'random_rows':random_rows,
        'normalized_linear_exact':str(normalized_linear),
        'normalized_occupancy_second_moment_exact':str(normalized_squares),
    }, normalized_linear, normalized_squares


def main():
    cases=[]
    for m in (1,2):
        linear_limit = F(m,m+2)
        square_limit = F(2**(2-m),m+2)
        errors=[]
        for H in (20,40,80):
            row,linear,squares=check(m,H)
            errors.append(abs(linear-linear_limit)+abs(squares-square_limit))
            row['linear_limit']=str(linear_limit)
            row['square_limit']=str(square_limit)
            cases.append(row)
            repaired,repaired_linear,repaired_squares=check(m,H,True)
            assert repaired_linear <= linear
            # A finite, non-asymptotic upper bound on adding the j=0 rows.
            # mu_a <= 2H (a loose bound) and at most H j=0 rows share a type.
            assert repaired_squares-squares <= F(5,2)
            cases.append(repaired)
        assert errors[2] < errors[1] < errors[0]
        assert errors[-1] < F(1,5)
    report={'status':'passed','arithmetic':'exact Fraction expectations; no Monte Carlo',
            'cases':cases,
            'scope':'Finite checks of coordinate constraints, second moments and convergence towards the proved limit; does not assert existence of a nonzero determinant term.'}
    Path(__file__).with_name('discrete_assignment_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','cases':len(cases)}))


if __name__=='__main__':
    main()
