#!/usr/bin/env python3
"""Exact arithmetic checks of the nonzero, unshifted periodic jet blocks."""
from fractions import Fraction
from math import comb, prod


def convolve(left, right, cutoff):
    result = [Fraction(0)] * (cutoff + 1)
    for x, a in enumerate(left):
        for y, b in enumerate(right[:cutoff + 1 - x]):
            result[x+y] += a*b
    return result


def pencil_value(center, power, cutoff):
    log_series = [center] + [
        Fraction((-1)**(j+1), j) for j in range(1, cutoff+1)
    ]
    result = [Fraction(1)] + [Fraction(0)]*cutoff
    for _ in range(power):
        result = convolve(result, log_series, cutoff)
    return result


def determinant(matrix):
    work = [row[:] for row in matrix]
    value = Fraction(1)
    for i in range(len(work)):
        pivot = next((j for j in range(i, len(work)) if work[j][i]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != i:
            work[i], work[pivot] = work[pivot], work[i]
            value = -value
        diagonal = work[i][i]
        value *= diagonal
        for j in range(i+1, len(work)):
            if work[j][i]:
                ratio = work[j][i]/diagonal
                for k in range(i+1, len(work)):
                    work[j][k] -= ratio*work[i][k]
    return value


def block_data(centers, cutoff, transverse_weight=2, time_weight=2):
    # theta = 1/(2*centers), w0=time_weight/2 as in the proof.
    beta_counts = [
        (cutoff - time_weight*s + 2*centers*transverse_weight - 1)
        // (2*centers*transverse_weight)
        for s in range((cutoff+time_weight-1)//time_weight)
    ]
    return beta_counts


def verify(centers, cutoff):
    step = Fraction(2)
    counts = block_data(centers, cutoff)
    rows = [(s, j, beta) for s,b in enumerate(counts)
            for j in range(centers) for beta in range(b)]
    columns = [(h, d) for h,b in enumerate(counts) for d in range(centers*b)]
    assert all(h+2*d <= cutoff for h,d in columns)
    newton, ordinary = [], []
    for s,j,beta in rows:
        newton_row, ordinary_row = [], []
        for h,d in columns:
            if d < beta:
                newton_row.append(Fraction(0))
                ordinary_row.append(Fraction(0))
                continue
            coeffs = pencil_value(j*step, d-beta, s)
            choose = comb(d,beta)
            newton_row.append(choose*coeffs[s-h] if h <= s else Fraction(0))
            ordinary_row.append(choose*sum(
                comb(h,k)*coeffs[s-k] for k in range(min(s,h)+1)
            ))
        newton.append(newton_row)
        ordinary.append(ordinary_row)
    exponent = sum(b*b for b in counts)
    base = prod(step*(j-i) for j in range(centers) for i in range(j))
    predicted = base**exponent
    assert determinant(newton) == determinant(ordinary) == predicted
    print(f"K={centers}, H={cutoff}: size {len(rows)}, "
          f"det={int(base)}^{exponent}, exact")


def verify_two_coordinates():
    K, cutoff, time_weight, light, heavy = 2, 16, 2, 2, 3
    rows = []
    columns = []
    exponent = 0
    step = Fraction(2)
    for s in range((cutoff + time_weight - 1)//time_weight):
        for eta in range(cutoff//(2*K*light)+1):
            budget = cutoff-time_weight*s-2*K*light*eta
            b = (budget+2*K*heavy-1)//(2*K*heavy)
            if b <= 0:
                continue
            rows.extend((s,j,eta,k) for j in range(K) for k in range(b))
            columns.extend((s,eta,d) for d in range(K*b))
            exponent += b*b

    matrix = []
    for s,j,eta,k in rows:
        row = []
        for h,gamma,d in columns:
            assert h+heavy*(gamma+d) <= cutoff
            if h > s or gamma != eta or d < k:
                row.append(Fraction(0))
            else:
                row.append(comb(d,k)*pencil_value(j*step,d-k,s-h)[s-h])
        matrix.append(row)
    assert determinant(matrix) == step**exponent
    print(f"Two normal coordinates: size {len(rows)}, "
          f"det={int(step)}^{exponent} in the transformed basis, exact.")


def relative_collision_mass(centers, cutoff):
    counts = block_data(centers, cutoff)
    populations = [centers*sum(b > beta for b in counts)
                   for beta in range(max(counts))]
    return Fraction(sum(n*n for n in populations), cutoff*sum(populations))


if __name__ == "__main__":
    for centers, cutoff in ((2,12), (3,18), (4,22)):
        verify(centers, cutoff)
    verify_two_coordinates()
    exact_mass = relative_collision_mass(20, 100000)
    limit = Fraction(20, 3)
    assert abs(exact_mass-limit) < Fraction(1, 50)
    print(f"Population test: {float(exact_mass):.8f}, limit {float(limit):.8f}.")
