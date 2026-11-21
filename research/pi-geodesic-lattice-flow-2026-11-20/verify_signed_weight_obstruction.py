"""Finite rational regression checks for the signed-weight integral obstruction."""
from collections import Counter
from fractions import Fraction as F
from math import comb, gcd


def multiply(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for j, x in enumerate(a):
        for k, y in enumerate(b):
            out[j+k] += x*y
    return out


def integrated_quotient(k, ell, h, u, v):
    polynomial = [0] * (4*k) + [(-1)**j * comb(4*ell, j) for j in range(4*ell+1)]
    weight = [-v, 4*v, -6*v, 4*v, u-v]
    for _ in range(2*h):
        polynomial = multiply(polynomial, weight)
    J = F(0)
    for i in range(len(polynomial)-1, 1, -1):
        coefficient = polynomial[i]
        polynomial[i-2] -= coefficient
        J += F(coefficient, i-1)
    assert not polynomial[1]
    return J, polynomial[0]


def atan_bounds(q, even_terms=120):
    lower = sum((F((-1)**j, (2*j+1)*q**(2*j+1)) for j in range(even_terms)), F(0))
    return lower, lower+F(1, (2*even_terms+1)*q**(2*even_terms+1))


def main():
    low5, high5 = atan_bounds(5)
    low239, high239 = atan_bounds(239)
    pi_low, pi_high = 16*low5-4*high239, 16*high5-4*low239
    assert F(3) < pi_low < pi_high < 4
    counts = Counter()
    for k in (1, 2):
        for ell in (1, 2):
            for h in (0, 1, 2):
                for u in range(-5, 6):
                    for v in range(-5, 6):
                        J, R = integrated_quotient(k, ell, h, u, v)
                        assert R == (-4)**ell * (u+4*v)**(2*h)
                        counts["coefficient checks"] += 1
                        if R == 0:
                            assert h > 0 and u+4*v == 0
                            counts["zero residue"] += 1
                            continue
                        rho = -4*J/R
                        lower = F(1, 24 * 4**(4*(k+ell)+ell) * 27**(2*h))
                        I_over_R = ((pi_low if R > 0 else -pi_high)
                                    + (4*J/abs(R)))/4
                        assert I_over_R >= lower
                        counts["integral bounds"] += 1
                        if not rho:
                            counts["zero rho"] += 1
                            continue
                        P, Q = (1/rho).as_integer_ratio()
                        # Test the intermediate implication on a rational bound for X.
                        X = 12 / lower
                        delta0, delta1 = P - Q/pi_low, P - Q/pi_high
                        assert delta0 * delta1 > 0
                        assert X*min(abs(delta0), abs(delta1)) > Q+2
                        a, r = 113, 355
                        out = (P+F(a, r))/(Q+1)
                        assert out.denominator == r*(Q+1)//gcd(r*P+a, Q+1)
                        counts["residual and reduced denominator checks"] += 1
    assert counts["coefficient checks"] == 1452
    print(dict(sorted(counts.items())))


if __name__ == "__main__":
    main()
