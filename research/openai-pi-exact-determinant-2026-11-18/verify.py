"""Exact checks of full-row determinant invariance and denominator budgets."""

from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import comb, factorial, lcm, prod


def mul(a, b, order):
    return [sum(a[j]*b[k-j] for j in range(k+1)) for k in range(order+1)]


def power(a, exponent, order):
    answer = [F(1)] + [F(0)]*order
    for _ in range(exponent):
        answer = mul(answer, a, order)
    return answer


def make_matrix(rows, columns, centers, degree_limits):
    matrix = []
    for j, s, beta in rows:
        line = []
        for h, alpha in columns:
            if any(b > a for a, b in zip(alpha, beta)):
                line.append(F(0))
                continue
            scalar = prod(comb(a,b) for a,b in zip(alpha,beta))
            pol = [F(comb(h,k)) if k <= h else F(0) for k in range(s+1)]
            for i, (a, b) in enumerate(zip(alpha,beta)):
                base = [j*centers[i]] + [
                    F((-1)**(k+1),k) if k <= degree_limits[i] else F(0)
                    for k in range(1,s+1)
                ]
                pol = mul(pol, power(base,a-b,s), s)
            line.append(scalar*pol[s])
        matrix.append(line)
    return matrix


def eliminate(source):
    matrix = [row.copy() for row in source]
    pivots = []
    determinant = F(1)
    for col in range(len(matrix[0])):
        row = len(pivots)
        pivot = next((k for k in range(row,len(matrix)) if matrix[k][col]), None)
        if pivot is None:
            continue
        if pivot != row:
            matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
            determinant = -determinant
        lead = matrix[row][col]
        determinant *= lead
        pivots.append(col)
        for i in range(row+1,len(matrix)):
            factor = matrix[i][col]/lead
            for k in range(col,len(matrix[i])):
                matrix[i][k] -= factor*matrix[row][k]
        if len(pivots) == len(matrix):
            return determinant, pivots
    return F(0), pivots


def valuation(value, prime):
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def primes_up_to(n):
    return [p for p in range(2,n+1) if all(p%d for d in range(2,p) if d*d <= p)]


def prime_budget(limits, q, differences, total_s):
    result = 1
    for prime in primes_up_to(max(limits)):
        costs = [
            max(max(0,valuation(k,prime)-valuation(qi,prime))
                for k in range(1,limit+1))
            for qi,limit in zip(q,limits)
        ]
        rate = max((F(max(0,valuation(k,prime)-valuation(qi,prime)),k)
                    for qi,limit in zip(q,limits) for k in range(1,limit+1)),
                   default=F(0))
        exponent = min(sum(x*y for x,y in zip(costs,differences)), total_s*rate//1)
        result *= prime**exponent
    return result


def check_case(height, widths, theta, v0, count, q, p):
    dim = len(widths)
    indices = list(product(*(range(height//w+1) for w in widths)))
    rows = [(j,s,b) for j in range(count) for s in range(height)
            for b in indices if v0*s + sum(w*k for w,k in zip(widths,b))/theta < height]
    rows.sort(key=lambda row: (row[0],v0*row[1]+sum(w*k for w,k in zip(widths,row[2]))/theta))
    cols = [(h,a) for h in range(height+1) for a in indices
            if h+sum(w*k for w,k in zip(widths,a)) <= height]
    short = [int(w/(v0*theta)) for w in widths]
    max_s = max(s for _,s,_ in rows)
    # Real centers give an exact-Q instance of the same local row operations.
    centers = [F(2*pi,qi) for pi,qi in zip(p,q)]
    arithmetic = make_matrix(rows,cols,centers,short)
    a, pivots = eliminate(arithmetic)
    if not a:
        raise ArithmeticError("Test matrix did not have full row rank")
    selected = [cols[j] for j in pivots]
    b, _ = eliminate(make_matrix(rows,selected,centers,[max_s]*dim))
    if a != b:
        raise ArithmeticError("Filtered substitutions changed a full-row determinant")
    differences = [sum(alpha[i] for _,alpha in selected)-sum(beta[i] for _,_,beta in rows)
                   for i in range(dim)]
    s_total = sum(s for _,s,_ in rows)
    d = prime_budget(short,q,differences,s_total)
    scaled = a*prod(qi**ci for qi,ci in zip(q,differences))
    if d % scaled.denominator:
        raise ArithmeticError("Exact denominator exceeds prime budget")
    fibers = Counter((s,beta[:-1]) for j,s,beta in rows if j == 0)
    exponent = sum(size*size for size in fibers.values())
    vandermonde = prod(factorial(j) for j in range(1,count))
    integer = d*scaled
    if integer.denominator != 1 or integer.numerator % vandermonde**exponent:
        raise ArithmeticError("Confluent Vandermonde divisibility failed")
    f0 = int(2/theta)+1
    old = prod(lcm(*range(1, (f0*w+v0-1)//v0))**((height//w)*len(rows))
               for w in widths)
    print(f"{len(rows)} rows, {len(cols)} columns, C={differences}, sum(s)={s_total}:")
    print(f"  denominator actual {scaled.denominator}; new cap has {len(str(d))} digits; old cap has {len(str(old))} digits")


    if vandermonde > 1:
        print(f"  cleared determinant divisible by {vandermonde}^{exponent}")


def distinct_degree_checks():
    for size in range(1,5):
        ratio = F(1,3)
        exact = factorial(size)*ratio**comb(size,2)/prod(1-ratio**j for j in range(1,size+1))
        bounded = factorial(size)*sum(ratio**sum(degrees)
                                      for degrees in product(range(10), repeat=size)
                                      if all(a < b for a,b in zip(degrees,degrees[1:])))
        if not (F(999,1000)*exact < bounded < exact):
            raise ArithmeticError("Distinct-degree identity failed its finite check")


if __name__ == "__main__":
    distinct_degree_checks()
    check_case(9,(3,4),F(1,2),1,2,(6,11),(19,35))
    check_case(8,(2,4),F(1,3),1,2,(7,10),(22,31))
    check_case(9,(3,4),F(1,2),1,2,(7,11),(22,35))
    check_case(12,(3,4),F(1,4),2,3,(7,10),(22,31))
