"""Exact checks of degree drops and reduced denominators of center minors."""

import json
import runpy
from fractions import Fraction
from itertools import product
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
ROW_TOOLS = runpy.run_path(str(ROOT / "openai-pi-vanishing-order-2026-11-18/verify.py"))
rank_and_determinant = ROW_TOOLS["rank_and_determinant"]
row_entry = ROW_TOOLS["row_entry"]


def verify_actual_minor(p1, q1, p2, q2):
    r1, r2 = Fraction(p1, q1), Fraction(p2, q2)
    columns = ((0,0), (0,1), (0,2), (1,0))
    rows = [
        [row_entry(j, s, (0,0), column, (r1,r2)) for column in columns]
        for j, s in product(range(2), repeat=2)
    ]
    rank, determinant = rank_and_determinant(rows)
    if rank != 4 or determinant != -2*r2*(r1-r2):
        raise ArithmeticError("Incorrect raw four-row determinant")
    # This minor is homogeneous of degree two. The actual r_i=2i*p_i/q_i
    # therefore multiply the real-center determinant by (2i)^2=-4.
    return -4*determinant


def verify_old_minor(p1, q1, p2, q2):
    centers = (Fraction(p1,q1), Fraction(p2,q2), Fraction(19,6))
    columns = ROW_TOOLS["raw_columns"](0) + ROW_TOOLS["raw_columns"](1)
    rows = [
        [row_entry(j, s, (0,0,b), column, centers) for column in columns]
        for b in range(2)
        for s in range(3-b)
        for j in range(4)
    ]
    rank, determinant = rank_and_determinant(rows)
    r1, r2, _ = centers
    if rank != 20 or determinant != -597196800*r1**22*(r1-r2)**10:
        raise ArithmeticError("Incorrect old 20-row determinant")
    return 2**32*determinant


def check_polynomial_bound():
    count = 0
    for d in range(1,4):
        for coefficient in product(range(-2,4), repeat=d+1):
            if coefficient[-1] <= 0:
                continue
            for q in range(2,16):
                p = 3*q + 1
                numerator = sum(a*p**k*q**(d-k) for k,a in enumerate(coefficient))
                reduced = Fraction(numerator, q**d).denominator
                c = gcd(q,coefficient[-1])
                if reduced % (q//c)**d or gcd(numerator,q) != c:
                    raise ArithmeticError("Failure in the exact gcd lemma")
                count += 1
    print(f"Integer-polynomial denominator checks: {count}.")


def check_consecutive_convergents(data):
    scale = 10**data["pi_grid_digits"]
    lo, hi = (Fraction(data[key],scale)
              for key in ("pi_lower_integer","pi_upper_integer"))
    previous, last = (0,1), (1,0)
    signed = []
    while len(signed) < 6:
        a = lo.numerator//lo.denominator
        if a != hi.numerator//hi.denominator:
            raise ArithmeticError("Partial quotient not interval-certified")
        previous, last = last, tuple(a*x+y for x,y in zip(last,previous))
        signed.append(last)
        lo, hi = 1/(hi-a), 1/(lo-a)
    vectors = signed[1:5]
    if vectors != [(22,7),(333,106),(355,113),(103993,33102)]:
        raise ArithmeticError("Unexpected convergent prefix")
    def bracket(i,j):
        p,q = vectors[i-1]
        r,s = vectors[j-1]
        return p*s-r*q
    values = (bracket(1,2)*bracket(3,4), bracket(1,3)*bracket(2,4),
              bracket(1,4)*bracket(2,3))
    if values != (1,-292,-293):
        raise ArithmeticError("Incorrect consecutive-convergent cross products")
    print(f"Four certified consecutive π convergents: {values}.")


if __name__ == "__main__":
    data = json.loads((ROOT / "pi-first-hit-window-2026-11-18/certificate.json").read_text())
    record = data["records"][-1]
    p, q = record["p"], record["q"]
    new = verify_actual_minor(355,113,p,q)
    old = verify_old_minor(355,113,p,q)
    if gcd(new.denominator,q**3) != q**2 or old.denominator != 113**32 * q**10:
        raise ArithmeticError("Unexpected reduced denominator at the 155-digit witness")
    print(f"155-digit witness: four-row denominator retains exactly q^2 (nominal q^3).")
    print(f"20-row denominator retains all q^10; total {len(str(old.denominator))} digits.")
    check_polynomial_bound()
    check_consecutive_convergents(data)
