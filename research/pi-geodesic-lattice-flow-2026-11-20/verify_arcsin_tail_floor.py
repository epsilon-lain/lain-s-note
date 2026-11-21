"""Reconstruct exact primal/dual certificates for the endpoint-tail relaxation.

Place beside verify_common_hypotenuse.py and solve_arcsin_integer_block.cpp.
Needs Python 3, g++ and GMP. Optional --exact prints the exact optima.
"""
import argparse
from fractions import Fraction
from math import comb, gcd, lcm
from pathlib import Path
import subprocess
import sys
import tempfile

from verify_common_hypotenuse import equations, triangle_data

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

# 0-based in the original triangle_data() order, with no reordering.
OMITTED = {
    11: (1,4,5,6,7,8,9,10,11,14,17,18,19,21,22,24,25,26,27,28,29,32,33,35,38,
         40,41,43,44,46,47,48,49,51,52,59,60,63,67,70,71,75,76,78,79,84,85,86,
         87,89,90,92,93,97,100,110,111,118,126,195),
    12: (1,4,5,6,7,8,9,10,11,12,14,17,18,19,21,22,24,25,26,27,28,29,30,32,33,
         35,38,40,41,43,44,46,47,48,49,51,52,57,59,60,63,67,70,71,74,76,78,79,
         85,86,87,89,93,94,97,100,110,111,124,126),
    105: (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
          19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
          37,38,39,40,41,42,43,45,46,47,48,49,50,51,52,53,54,55,
          56,57,58,60,63,67),
    106: (0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
          19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
          37,38,39,40,41,42,43,45,46,47,48,49,50,51,52,53,54,55,
          56,57,58,60,63,67),
}


def clear(values):
    common = lcm(*(v.denominator for v in values))
    return [v.numerator * (common // v.denominator) for v in values], common


def dot(row, vector):
    return sum(a*b for a,b in zip(row,vector))


def solve(binary, matrix, rhs):
    size = len(rhs)
    lines = [str(size)]
    for row, value in zip(matrix, rhs):
        assert len(row) == size
        lines.append(" ".join(map(str, [*row, value])))
    result = subprocess.run([binary], input="\n".join(lines)+"\n", text=True,
                            capture_output=True, check=True)
    solution = [Fraction(s) for s in result.stdout.splitlines()]
    assert len(solution) == size
    return clear(solution)


def verify_order(n, binary, rows, hypotenuse, points, show_exact):
    power = 2*n + 3
    coeff = Fraction(comb(2*n+2,n+1), 4**(n+1)*power)
    def tail(t):
        return t**power/(1-t*t)
    costs = [coeff*(tail(Fraction(x,hypotenuse)) + tail(Fraction(y,hypotenuse)))
             for x,y in points]
    cost, cost_den = clear(costs)
    support = [i for i in range(len(points)) if i not in OMITTED[n]]
    block = [[row[i] for i in support] for row in rows]
    rhs = [0]*(len(rows)-1) + [2]
    assert len(support) == len(rows) == 161 and len(points) == 221

    weight, weight_den = solve(binary, block, rhs)
    assert all(w != 0 for w in weight)
    assert all(dot(row, weight) == b*weight_den for row,b in zip(block,rhs))
    objective = Fraction(sum(abs(w)*cost[i] for w,i in zip(weight,support)),
                         weight_den*cost_den)

    signed_cost = [cost[i] if w>0 else -cost[i] for i,w in zip(support,weight)]
    dual, dual_den = solve(binary, list(zip(*block)), signed_cost)
    bounds = [cost[i]*dual_den - abs(dot(column,dual))
              for i,column in enumerate(zip(*rows))]
    assert min(bounds) >= 0
    assert all(bounds[i] == 0 for i in support)
    assert all(bounds[i] > 0 for i in OMITTED[n])
    assert objective == Fraction(dot(rhs,dual), dual_den*cost_den)
    assert objective.denominator > 1
    floor_value = objective.numerator // objective.denominator
    print(f"N={n}: {floor_value} < minimum U < {floor_value+1}; "
          f"{len(support)} equalities, {len(OMITTED[n])} strict dual inequalities")
    if n in (105,106):
        scaled = (1000*objective.numerator) // objective.denominator
        assert objective != Fraction(scaled,1000)
        print(f"  refined interval: {scaled}/1000 < minimum U < {scaled+1}/1000")
    if show_exact:
        print(objective)


def verify_support_radius(binary, rows, hypotenuse, points):
    order = sorted(range(len(points)), key=lambda i: points[i][1])
    assert len({y for _,y in points}) == len(points)
    rhs = [0]*(len(rows)-1) + [2]
    denoms = []
    for start in (0,1):
        indices = order[start:start+len(rows)]
        block = [[row[i] for i in indices] for row in rows]
        nums, den = solve(binary, block, rhs)
        assert gcd(den,*nums) == 1
        assert all(dot(row,nums) == b*den for row,b in zip(block,rhs))
        denoms.append(den)
    assert denoms[0] > 1 and gcd(*denoms) == 1
    radius = Fraction(points[order[161]][1],hypotenuse)
    assert radius == Fraction(877414027070304,1021090952484265)
    print("least largest endpoint of any feasible integer support:", radius)
    for q_floor in (10**52,10**54):
        # At N=0, c_1=1/6. Advance the one compulsory endpoint's contribution.
        value = q_floor**2 * radius**3 / (6*(1-radius**2))
        n = 0
        while value > 1:
            k = n+1
            value *= Fraction((2*k+1)**2, 2*(k+1)*(2*k+3)) * radius**2
            n += 1
        print(f"Q >= 10^{len(str(q_floor))-1}, Q^2 U <= 1 requires N >= {n}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact", action="store_true")
    args = parser.parse_args()
    primes, hypotenuse, points = triangle_data()
    rows = equations(primes,points)
    source = Path(__file__).resolve().with_name("solve_arcsin_integer_block.cpp")
    with tempfile.TemporaryDirectory() as folder:
        binary = Path(folder)/"solve"
        subprocess.run(["g++","-O2",str(source),"-lgmpxx","-lgmp","-o",str(binary)], check=True)
        verify_support_radius(binary,rows,hypotenuse,points)
        for n in OMITTED:
            verify_order(n, binary, rows, hypotenuse, points, args.exact)


if __name__ == "__main__":
    main()
