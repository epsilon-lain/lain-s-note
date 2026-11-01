"""Independent exhaustive minor checks; rational rank-area examples."""
from fractions import Fraction as F
from itertools import combinations, permutations
from random import Random
from math import factorial
import json
from rank_minor import cheapest_minor


def leibniz(matrix):
    n = len(matrix)
    result = F(0)
    for perm in permutations(range(n)):
        sign = (-1) ** sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
        value = F(sign)
        for i in range(n):
            value *= matrix[i][perm[i]]
        result += value
    return result


def main():
    rng = Random(20261102)
    count = full_rank = deficient = 0
    for m in (1, 2, 3):
        for n in (m, m + 1, m + 3):
            for _ in range(20):
                matrix = [[F(rng.randint(-2, 2), rng.randint(1, 3)) for j in range(n)] for i in range(m)]
                costs = [F(rng.randint(0, 5), 6) for j in range(n)]
                expected = []
                for selected in combinations(range(n), m):
                    minor = [[row[j] for j in selected] for row in matrix]
                    if leibniz(minor):
                        expected.append(sum((costs[j] for j in selected), F(0)))
                got = cheapest_minor(matrix, costs)
                if expected:
                    assert got['minimum_cost'] == min(expected)
                    selected = got['chosen']
                    assert got['determinant'] == leibniz([[row[j] for j in selected] for row in matrix])
                    # Layer-cake certificate over [0,1], including the last gap.
                    area, left, rank = F(0), F(0), 0
                    for cost, new_rank in got['profile']:
                        area += (cost - left) * (1 - F(rank, m))
                        left, rank = cost, new_rank
                    area += (1 - left) * (1 - F(rank, m))
                    assert area == got['minimum_cost'] / m
                    full_rank += 1
                else:
                    assert got['minimum_cost'] is None and got['rank'] < m
                    deficient += 1
                count += 1

    # Sharpness: full rank alone permits only the most expensive columns.
    costs = [F(i, 5) for i in range(6)]
    sharp = cheapest_minor([[0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]], costs)
    assert sharp['chosen'] == [4, 5] and sharp['minimum_cost'] == F(9, 5)

    # Exact rational integral of 1-min(1,(t/a)^(m+1)).
    area_checks = 0
    for m in range(1, 9):
        for a in (F(1, 4), F(1, 2), F(3, 4)):
            integral = a - a ** (m + 2) / ((m + 2) * a ** (m + 1))
            assert integral == a * (m + 1) / (m + 2)
            area_checks += 1
    z, rho = F(3, 4), F(5, 32)
    assert 3*z*z - 2*z**3 == 1-rho
    old_bound = (F(1, 2) - 2*z**3 + F(3, 2)*z**4) / rho
    new_bound = z * F(3, 4)
    assert old_bound == F(67, 80) and new_bound == F(9, 16)
    assert old_bound - new_bound == F(11, 40)
    print(json.dumps(dict(status='passed', exhaustive_matrices=count,
                         full_rank=full_rank, rank_deficient=deficient,
                         exact_rank_area_formulas=area_checks,
                         cost_bound_example={'old': str(old_bound), 'new': str(new_bound),
                                             'difference': str(old_bound-new_bound)},
                         worst_case_sharpness=dict(selected=sharp['chosen'], cost=str(sharp['minimum_cost'])),
                         scope='Finite Q-matrix certificates; asymptotic interpolation theorem is not numerically certified.'), indent=2))


if __name__ == '__main__':
    main()
