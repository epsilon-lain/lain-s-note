#!/usr/bin/env python3
"""Exact finite tests of the pinned pi017 interpolation mechanism.

These are ablations and consistency checks, not a proof of its asymptotic
separated-weight interpolation theorem. All arithmetic is rational.
"""
from fractions import Fraction
from math import comb, factorial
import json
from pathlib import Path
import sympy as sp


def jet_matrix(H, w0, w1, v0, v1, centers):
    """P(Y,X) -> [t^s u^b]P(1+t,c+u+log(1+t))."""
    H, w0, w1, v0, v1 = map(Fraction, (H, w0, w1, v0, v1))
    columns = [(h, a) for h in range(int(H/w0)+1)
               for a in range(int((H-w0*h)/w1)+1)]
    indices = [(s, b) for s in range(int(H/v0)+1)
               for b in range(int(H/v1)+1) if v0*s+v1*b < H]
    max_s = max(s for s, _ in indices)
    # log(1+t)^p, truncated exactly, for every potentially relevant p.
    log = [sp.Rational(0)] + [sp.Rational((-1)**(k+1), k)
                                      for k in range(1, max_s+1)]
    lp = [[sp.Rational(1)]+[sp.Rational(0)]*max_s]
    for _ in range(max_s):
        previous = lp[-1]
        lp.append([sum(previous[i]*log[k-i] for i in range(k+1))
                   for k in range(max_s+1)])

    def coeff(h, a, c, s, b):
        if b > a:
            return sp.Rational(0)
        c = sp.Rational(c)
        n = a-b
        # binom(a,b) (1+t)^h (c+log(1+t))^(a-b).
        return comb(a,b)*sum(
            comb(n,p)*c**(n-p)*sum(
                comb(h,l)*lp[p][s-l] for l in range(min(h,s)+1))
            for p in range(min(n,max_s)+1))
    rows = [(j,s,b) for j in range(len(centers)) for s,b in indices]
    matrix = sp.Matrix([[coeff(h,a,centers[j],s,b) for h,a in columns]
                        for j,s,b in rows])
    return matrix, rows, columns


def check_direct_series():
    Y,X,t,u = sp.symbols('Y X t u')
    matrix, rows, columns = jet_matrix(2,1,Fraction(1,10),Fraction(1,2),Fraction(2,5),[0])
    # Verify selected mixed coefficients directly against symbolic substitution.
    for ci,(h,a) in enumerate(columns):
        if a > 5:
            continue
        direct = sp.expand((1+t)**h*(u+sp.log(1+t))**a).series(t,0,4).removeO().expand()
        for ri,(_,s,b) in enumerate(rows):
            assert matrix[ri,ci] == direct.coeff(t,s).coeff(u,b)
    # At t=0, normalized t-derivatives equal falling factorials of D0.
    P = Y**2*X**3+Y*X+X**2
    D0 = lambda f: sp.expand(Y*sp.diff(f,Y)+sp.diff(f,X))
    falling = P
    direct = sp.expand(P.subs({Y:1+t,X:sp.Rational(2,3)+u+sp.log(1+t)})).series(t,0,5).removeO().expand()
    for s in range(5):
        if s:
            falling = D0(falling)-(s-1)*falling
        expected = (falling/factorial(s)).subs({Y:1,X:sp.Rational(2,3)+u}).expand()
        assert sp.expand(expected-direct.coeff(t,s)) == 0


def modular_rank(matrix, prime):
    """A nonzero minor modulo a prime certifies the same rational rank lower bound."""
    a = [row[:] for row in matrix]
    rank = 0
    for column in range(len(a[0])):
        pivot = next((i for i in range(rank,len(a)) if a[i][column] % prime),None)
        if pivot is None:
            continue
        a[rank],a[pivot] = a[pivot],a[rank]
        inv = pow(a[rank][column],-1,prime)
        a[rank] = [(x*inv)%prime for x in a[rank]]
        for i in range(rank+1,len(a)):
            factor = a[i][column]
            if factor:
                a[i] = [(x-factor*y)%prime for x,y in zip(a[i],a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def separated_t_only_test(H, prime=1000003):
    # Genuine source separation, m=K=1, sigma=1/10, w1=31>30.
    # For H<124, no transverse row occurs. z=log(1+t) is an invertible
    # triangular jet change. Use d^n/dz^n[e^(h*z)z^a] at z=0 modulo p.
    assert H < 124 and 2*H < prime
    columns = [(h,a) for h in range(H+1) for a in range((H-h)//31+1)]
    matrix = [[0 if n<a else factorial(n)//factorial(n-a)*pow(h,n-a,prime)%prime
               for h,a in columns] for n in range(2*H)]
    rank = modular_rank(matrix,prime)
    assert rank == min(2*H,len(columns))
    return {'H':H,'rows':2*H,'columns':len(columns),
            'certified_rank':rank,'certificate_prime':prime}


def main():
    check_direct_series()
    first, rows, columns = jet_matrix(3,1,1,2,2,[0,1,2])
    assert first.shape == (9,10)
    assert first.rank() == 7
    results = {
        'checks': 'direct Taylor coefficient and falling-D0 identities pass',
        'ambient_volume_only': {
            'H':3, 'W':['1','1'], 'V':['2','2'], 'K':3,
            'theta':'1/2', 'ambient_volume_ratio':'3/4',
            'fiber_volume_ratio':'3/2',
            'rows':first.rows,'columns':first.cols,'rank':first.rank(),
            'kernel_dimension':len(first.nullspace()),
        },
        'both_volume_conditions_but_unseparated': [],
        'illustrative_increased_weight': [],
        'genuinely_separated_t_only': [separated_t_only_test(50),separated_t_only_test(100)],
    }
    for H in range(1,6):
        matrix, rows, columns = jet_matrix(H,1,Fraction(1,10),Fraction(1,2),Fraction(2,5),[0])
        rank = matrix.rank()
        assert rank < matrix.rows
        results['both_volume_conditions_but_unseparated'].append(
            {'H':H,'rows':matrix.rows,'columns':matrix.cols,'rank':rank,
             'rank_deficit':matrix.rows-rank})
        repaired, _, _ = jet_matrix(H,1,1,Fraction(1,2),4,[0])
        # This finite check does not certify the source's large threshold.
        assert repaired.rank() == repaired.rows
        results['illustrative_increased_weight'].append(
            {'H':H,'rows':repaired.rows,'columns':repaired.cols,'rank':repaired.rank()})
    output = Path(__file__).with_name('interpolation_lab_results.json')
    output.write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps(results,indent=2))

if __name__ == '__main__':
    main()
