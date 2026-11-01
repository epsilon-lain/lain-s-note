"""Exact checks of boundary coefficient budgets and compatible columns.
No numerical extrapolation is used as a proof of an infinite assertion.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, prod
from pathlib import Path
import json


def determinant(a):
    a = [list(map(F, row)) for row in a]
    result = F(1)
    for k in range(len(a)):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return F(0)
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            result = -result
        p = a[k][k]
        result *= p
        for i in range(k+1, len(a)):
            ratio = a[i][k]/p
            for j in range(k+1, len(a)):
                a[i][j] -= ratio*a[k][j]
    return result


def verify_compatible():
    cases = []
    for n in range(1, 8):
        for r in (F(0), F(6), F(601, 100), F(-7, 3)):
            b = [[F(comb(h,s)) if s <= h else F(0)
                  for h in range(n)] for s in range(n)]
            l = [[sum(F(comb(h,s-k)*(-1)**(k+1), k)
                      for k in range(1,s+1) if s-k <= h)
                  for h in range(n)] for s in range(n)]
            matrix = [b[s]+l[s] for s in range(n)]
            matrix += [b[s]+[l[s][h]+r*b[s][h] for h in range(n)]
                       for s in range(n)]
            result = determinant(matrix)
            assert result == r**n
            cases.append({'n': n, 'r': str(r), 'determinant': str(result)})
    return cases


def verify_coefficients():
    count = 0
    nonzero = 0
    row_sums = 0
    for weights, h, v0 in [((2,), 15, 1), ((2,3), 15, 1), ((3,5), 24, 2)]:
        # theta=1/2, so vi=2wi; all ci are positive integers.
        vi = tuple(2*w for w in weights)
        costs = tuple(v-v0 for v in vi)
        t = F(1,3)
        zeta = tuple(t**c for c in costs)
        types = [a for a in product(*(range(h//w+1) for w in weights))
                 if sum(w*x for w,x in zip(weights,a)) <= h]
        for beta in types:
            for s in range(h//v0+1):
                source_cost = v0*s+sum(v*b for v,b in zip(vi,beta))
                if source_cost >= h:
                    continue
                gap = h-source_cost
                coefficients = []
                for target in types:
                    if any(a<b for a,b in zip(target,beta)):
                        continue
                    gamma = tuple(a-b for a,b in zip(target,beta))
                    if sum(gamma)>s or sum(c*g for c,g in zip(costs,gamma))<gap:
                        continue
                    alt = sum((-1)**sum(k)*prod(comb(g,j) for g,j in zip(gamma,k))
                              for k in product(*(range(g+1) for g in gamma))
                              if sum(c*j for c,j in zip(costs,k))<gap)
                    coefficient = prod(comb(a,b) for a,b in zip(target,beta))
                    coefficient *= prod((-z)**g for z,g in zip(zeta,gamma))*alt
                    # tau=-log(t) exactly, so exp(-tau*gap)=t**gap.
                    bound = 3**(h//min(weights))*t**gap
                    assert abs(coefficient) <= bound
                    coefficients.append(abs(coefficient))
                    count += 1
                    nonzero += bool(coefficient)
                assert sum(coefficients) <= len(types)*3**(h//min(weights))*t**gap
                row_sums += 1
    return {'boundary_coefficients': count, 'nonzero_coefficients': nonzero,
            'row_sum_bounds': row_sums}


def verify_gap_counts():
    cases = []
    for d in range(1,5):
        for h in range(1,10):
            lattice = [a for a in product(range(h),repeat=d) if sum(a)<h]
            count = comb(h+d-1,d)
            total_gap = sum(h-sum(a) for a in lattice)
            assert len(lattice)==count
            assert F(total_gap,count)==F(h+d,d+1)
            cases.append({'coordinates':d,'H':h,'rows':count,
                          'normalized_mean_gap':str(F(total_gap,count*h))})
    return cases


def verify_compatibility_equations():
    cases=0
    for alpha in product(range(5),repeat=3):
        for a in product(*(range(x+1) for x in alpha)):
            ba=prod(comb(x,y) for x,y in zip(alpha,a))
            left=(sum(alpha)-sum(a))*ba
            right=0
            for i in range(3):
                if a[i]<alpha[i]:
                    ap=list(a);ap[i]+=1
                    right+=(a[i]+1)*prod(comb(x,y) for x,y in zip(alpha,ap))
            assert left==right
            cases+=1
    return cases


if __name__=='__main__':
    result={'status':'passed','compatible_column_determinants':verify_compatible(),
            'boundary_budget_checks':verify_coefficients(),
            'simplex_gap_checks':verify_gap_counts(),
            'differential_coefficient_identities':verify_compatibility_equations(),
            'scope':'Exact rational finite identities and coefficient bounds. No pi approximation property, high-dimensional threshold, or uniform mixed-minor estimate is certified.'}
    (Path(__file__).parent/'budget_compatibility_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:(len(v) if isinstance(v,list) else v) for k,v in result.items()},indent=2))
