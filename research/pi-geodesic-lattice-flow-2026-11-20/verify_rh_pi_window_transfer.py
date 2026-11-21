"""Finite cross-checks for the Gaussian divisor sum and the jet divisor.

The Gaussian calculations are floating-point diagnostics. The jet calculation
uses exact integer arithmetic. Neither loop tests an asymptotic assertion.
"""
from collections import defaultdict
from itertools import combinations
from math import comb, exp, fsum, gcd, log, pi


def variance_by_products(scale, cutoff):
    weights = [exp(-pi * (n / scale)**2) for n in range(cutoff + 1)]
    coefficients = defaultdict(float)
    for a in range(1, cutoff + 1):
        for b in range(1, cutoff + 1):
            coefficients[a*b] += weights[a] * weights[b]
    return 8 * fsum(t*t for t in coefficients.values()) / scale**2


def variance_by_coprime_pairs(scale, cutoff):
    def contribution(a, b):
        rate = pi * (a*a + b*b) / scale**2
        return fsum(exp(-rate * c*c) for c in range(1, cutoff // max(a, b) + 1))**2

    square_sum = fsum(contribution(a, b)
                      for a in range(1, cutoff + 1)
                      for b in range(1, cutoff + 1) if gcd(a, b) == 1)
    return 8 * square_sum / scale**2


def determinant(rows):
    mat = [list(row) for row in rows]
    n = len(mat)
    sign, divisor = 1, 1
    for j in range(n - 1):
        pivot = next((r for r in range(j, n) if mat[r][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            mat[pivot], mat[j] = mat[j], mat[pivot]
            sign = -sign
        a = mat[j][j]
        for i in range(j+1, n):
            for k in range(j+1, n):
                cross = a*mat[i][k] - mat[i][j]*mat[j][k]
                mat[i][k], remainder = divmod(cross, divisor)
                assert remainder == 0
            mat[i][j] = 0
        divisor = a
    return sign * mat[-1][-1]


def truncated_product(a, b, order):
    return [sum(a[i]*b[j-i] for i in range(j+1)) for j in range(order+1)]


def example_jet_minor_gcd():
    # Rows have j=0,1 and s+2*beta<3. The center c=2 and twist t+t^2
    # give the columns [t^s u^beta] (1+t)^h (2*j+t+t^2+u)^a, h+a<=3.
    rows = [(j, s, beta) for j in range(2) for s in range(3)
            for beta in range(2) if s + 2*beta < 3]
    columns = [(h, a) for h in range(4) for a in range(4-h)]
    matrix = []
    for j, s, beta in rows:
        row = []
        for h, a in columns:
            if beta > a:
                row.append(0)
                continue
            coeff = [1, 0, 0]
            for _ in range(a-beta):
                coeff = truncated_product(coeff, [2*j, 1, 1], 2)
            value = sum(coeff[s-i]*comb(h, i) for i in range(min(s, h)+1))
            row.append(comb(a, beta) * value)
        matrix.append(row)
    all_minors = (determinant([[row[j] for j in ids] for row in matrix])
                  for ids in combinations(range(len(columns)), len(rows)))
    return gcd(*all_minors)


def main():
    for scale in (8, 16, 32, 64):
        cutoff = 3 * scale
        products = variance_by_products(scale, cutoff)
        coprime = variance_by_coprime_pairs(scale, cutoff)
        assert abs(products - coprime) < 1e-10
        print(f"X={scale:3d}: variance={products:.10f}, "
              f"variance-(6/pi)log(X)={products - (6/pi)*log(scale):.10f}")
    divisor = example_jet_minor_gcd()
    assert divisor and divisor % (2**6) == 0 and divisor % (2**7) != 0
    print(f"8 x 10 jet matrix: gcd of full-row minors = {divisor}, v2 = 6")


if __name__ == "__main__":
    main()
