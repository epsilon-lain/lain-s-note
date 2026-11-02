"""Exact audits of DIRECT_INTEGRAL_GCD_THRESHOLD.md; standard library only.

Finite checks do not establish the Laplace asymptotic or an infinite gcd
subsequence. Small integrals are independently reduced to rational + pi/4.
"""
from fractions import Fraction as F
from math import comb, factorial, gcd
from pathlib import Path
import json


def solve(matrix, rhs):
    a = [list(map(F, row)) + [F(b)] for row, b in zip(matrix, rhs)]
    size = len(a)
    for col in range(size):
        pivot = next(i for i in range(col, size) if a[i][col])
        a[col], a[pivot] = a[pivot], a[col]
        d = a[col][col]
        a[col] = [x/d for x in a[col]]
        for row in range(size):
            if row != col and a[row][col]:
                d = a[row][col]
                a[row] = [x-d*y for x, y in zip(a[row], a[col])]
    return [row[-1] for row in a]


def hermite_integral(n):
    """Solve N=(1+x²)S'-2nxS+c(1+x²)^n independently of recurrence."""
    size = 2*n+1
    matrix = [[0]*size for _ in range(size)]
    rhs = [0]*size
    for j in range(2*n):
        if j:
            matrix[j-1][j] += j
        matrix[j+1][j] += j-2*n
    for k in range(n+1):
        matrix[2*k][-1] = comb(n, k)
        rhs[n+k] = (-1)**k*comb(n, k)
    answer = solve(matrix, rhs)
    s, c = answer[:-1], answer[-1]
    # Check the polynomial identity again, with the solved coefficients.
    assert all(sum(F(x)*y for x,y in zip(row,answer)) == b
               for row,b in zip(matrix,rhs))
    rational = sum(s)/2**n-s[0]
    return c, rational


def arctan_interval(a, terms):
    total = sum((F((-1)**k, (2*k+1)*a**(2*k+1))
                 for k in range(terms)), F(0))
    next_term = F((-1)**terms, (2*terms+1)*a**(2*terms+1))
    return min(total,total+next_term), max(total,total+next_term)


def pi_interval():
    # Exact tangent identity: (5+i)^4(239-i) has equal positive parts.
    def multiply(z,w):
        return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])
    z = (1,0)
    for _ in range(4):
        z = multiply(z,(5,1))
    z = multiply(z,(239,-1))
    assert z == (114244,114244)
    # The relevant angle is between 0 and 1; pi/2>1 follows from
    # pi=4 integral_0^1 1/(1+x²) dx>2, so this selects the correct branch.
    assert 4*(F(1,5)-F(1,3*5**3))-F(1,239) > 0
    assert F(4,5) < 1
    a,b = arctan_interval(5,380),arctan_interval(239,120)
    lo,hi = 16*a[0]-4*b[1],16*a[1]-4*b[0]
    assert F(3) < lo < hi < F(22,7)
    assert hi-lo < F(1,10**520)
    return lo,hi


def decimal_enclosure(lo,hi,places=10):
    scale = 10**places
    return {"lower_scaled_integer":(lo*scale).__floor__(),
            "upper_scaled_integer":-(-hi*scale).__floor__(),
            "decimal_places":places}


def main():
    max_n = 100
    a,b = [1,1],[0,1]
    for n in range(1,max_n):
        a.append((2*n+1)*a[n]+n*n*a[n-1])
        b.append((2*n+1)*b[n]+n*n*b[n-1])
    ell = [F(a[n],factorial(n)) for n in range(max_n+1)]
    w = [F(b[n],factorial(n)) for n in range(max_n+1)]
    for n in range(1,max_n):
        assert (n+1)*ell[n+1] == (2*n+1)*ell[n]+n*ell[n-1]
        assert (n+1)*w[n+1] == (2*n+1)*w[n]+n*w[n-1]
    integral_checks = 0
    for n in range(1,13):
        c,rational = hermite_integral(n)
        assert (-2)**n*c == ell[n]
        assert (-2)**n*rational == -w[n]
        integral_checks += 1
    plo,phi = pi_interval()
    gs = [gcd(x,4*y) for x,y in zip(a,b)]
    samples = []
    all_quality_bounds = []
    for n in range(max_n+1):
        g = gs[n]
        p,q = 4*b[n]//g,a[n]//g
        assert q>0 and gcd(p,q)==1
        assert F(p,q)==4*w[n]/ell[n]
        assert F(q,1)/ell[n]==F(factorial(n),g)
        r = F(p,q)
        if n%2:
            elo,ehi = r-phi,r-plo
        else:
            elo,ehi = plo-r,phi-r
        assert 0<elo<ehi
        qlo,qhi = q*q*elo,q*q*ehi
        all_quality_bounds.append((qlo,qhi))
        assert qlo/F(factorial(n),g)**2==ell[n]**2*elo
        assert qhi/F(factorial(n),g)**2==ell[n]**2*ehi
        if n in (1,2,4,8,16,32,64,100):
            samples.append({"n":n,"p":str(p),"q":str(q),"gcd":str(g),
                            "Q":decimal_enclosure(qlo,qhi),
                            "analytic_normalization_n_ell_squared_error":
                                decimal_enclosure(n*ell[n]**2*elo,n*ell[n]**2*ehi),
                            "h_squared":str(F(n*g*g,factorial(n)**2))})
    for n in range(1,max_n+1):
        det = a[n]*b[n-1]-a[n-1]*b[n]
        assert det==(-1)**n*factorial(n-1)**2
        assert (4*factorial(n-1)**2)%(gs[n]*gs[n-1])==0
    # Identify the best finite tested value by disjoint rigorous intervals.
    best = min(range(1,max_n+1),key=lambda n:all_quality_bounds[n][1])
    assert all(all_quality_bounds[best][1]<all_quality_bounds[n][0]
               for n in range(1,max_n+1) if n!=best)
    assert best==2
    result = {
        "status":"passed","arithmetic":"Python standard-library exact rationals",
        "maximum_n":max_n,
        "independent_rational_integral_reductions":integral_checks,
        "rational_recurrence_checks":2*(max_n-1),
        "reduced_fraction_and_joint_cost_checks":max_n+1,
        "integer_wronskian_checks":max_n,
        "adjacent_gcd_product_divisibility_checks":max_n,
        "pi_interval_width_less_than":"1e-520",
        "finite_best_n":best,
        "finite_best_fraction":"3/1",
        "finite_best_quality":decimal_enclosure(*all_quality_bounds[best]),
        "samples":samples,
        "scope":"Exact finite algebraic and error certificates. The Laplace asymptotic, gcd threshold equivalence and limiting neighbor bound are proved in the note, not inferred from these samples. No infinite subsequence or pi bad-approximability conclusion."
    }
    Path(__file__).with_name('direct_integral_verification.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='samples'},indent=2))


if __name__=='__main__':
    main()
