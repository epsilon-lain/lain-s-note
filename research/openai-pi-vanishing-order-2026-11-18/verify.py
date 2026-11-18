"""Exact verification of a full weighted row packet with ten vanishing orders."""

from fractions import Fraction
from math import comb


def multiply(left, right, order):
    return [
        sum(left[j] * right[k-j] for j in range(k+1))
        for k in range(order+1)
    ]


def row_entry(center, s, beta, exponent, r):
    if any(a < b for a, b in zip(exponent, beta)):
        return Fraction(0)
    result = [Fraction(1)] + [Fraction(0)]*s
    factor = 1
    for a, b, coordinate in zip(exponent, beta, r):
        factor *= comb(a,b)
        base = [center*coordinate] + [Fraction((-1)**(k+1),k) for k in range(1,s+1)]
        for _ in range(a-b):
            result = multiply(result,base,s)
    return factor*result[s]


def rank_and_determinant(rows):
    a = [[Fraction(value) for value in row] for row in rows]
    rank = 0
    determinant = Fraction(1)
    for col in range(len(a[0])):
        pivot = next((k for k in range(rank,len(a)) if a[k][col]), None)
        if pivot is None:
            continue
        if pivot != rank:
            a[pivot], a[rank] = a[rank], a[pivot]
            determinant = -determinant
        leading = a[rank][col]
        determinant *= leading
        for i in range(rank+1,len(a)):
            factor = a[i][col]/leading
            a[i][col:] = [x-factor*y for x,y in zip(a[i][col:],a[rank][col:])]
        rank += 1
    return rank, determinant if rank == len(a) else Fraction(0)


def raw_columns(b):
    times = 3-b
    extra = (4,5,6) if b == 0 else (4,)
    return [(a,c,b) for c in range(times) for a in range(4-c)] + [(a,0,b) for a in extra]


def check():
    theta = Fraction(1,12)
    time_weight = Fraction(24,5)
    height = 12
    weights = (1,1,Fraction(1,2))
    rows = []
    for b in range(2):
        for s in range(3):
            if time_weight*s + (weights[2]*b)/theta < height:
                rows.extend((j,s,(0,0,b)) for j in range(4))
    # This enumeration is complete: beta_1, beta_2, s>=3 or beta_3>=2
    # each individually forces a weight at least H.
    columns = raw_columns(0) + raw_columns(1)
    if len(rows) != len(columns) or any(
        sum(w*a for w,a in zip(weights,column)) > height for column in columns
    ):
        raise ArithmeticError("Incorrect complete row packet or column weights")
    cases = ((1,2,3),(2,-1,5),(1,1,3),(2,2,2),
             (Fraction(22,7),Fraction(355,113),Fraction(19,6)),
             (Fraction(1,2),0,Fraction(17,3)))
    for raw_r in cases:
        r = tuple(map(Fraction,raw_r))
        entries = [[row_entry(j,s,beta,a,r) for a in columns] for j,s,beta in rows]
        rank, determinant = rank_and_determinant(entries)
        expected = -597196800*r[0]**22*(r[0]-r[1])**10
        expected_rank = 12 if r[0] == r[1] else 20
        if determinant != expected or rank != expected_rank:
            raise ArithmeticError("Exact determinant or rank disagreement")
    orders = [a[1] for a in columns]
    shifted_profiles = [sum(n >= k for n in orders) for k in (1,2,3)]
    if shifted_profiles != [8,2,0]:
        raise ArithmeticError("Incorrect column filtration sizes")
    print("Six exact 20-by-20 row packets verified.")
    print("Rank on diagonal 12, order 10; filtration dimensions 8, 2, 0.")
    p1,q1,p2,q2 = 22,7,355,113
    cross = p1*q2-p2*q1
    if abs(cross) != 1:
        raise ArithmeticError("Incorrect cross determinant for the test fractions")
    print(f"22/7 and 355/113 have cross determinant {cross}; the tenth power is 1.")


if __name__ == "__main__":
    check()
