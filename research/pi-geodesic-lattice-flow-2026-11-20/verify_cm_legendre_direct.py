#!/usr/bin/env python3
"""Check direct and total-degree truncations of the square-CM Legendre bridge."""
from fractions import Fraction
from math import comb, log10

from cm_legendre_pi_bridge_experiment import (
    A, P, count_points_mod_p, mul, mul_with_cocycle, xi,
)


def atan_bounds(base, terms=60):
    partial = sum(Fraction((-1)**r, (2*r+1)*base**(2*r+1))
                  for r in range(terms))
    neighbor = partial + Fraction((-1)**terms, (2*terms+1)*base**(2*terms+1))
    return min(partial, neighbor), max(partial, neighbor)


def pi_bounds():
    low5, high5 = atan_bounds(5)
    low239, high239 = atan_bounds(239)
    return 16*low5 - 4*high239, 16*high5 - 4*low239


def local_coefficients(z, cutoff):
    logs = [z]
    zetas = [1/z]
    for r in range(1, cutoff + 1):
        power = (-1)**r * comb(2*r, r) * A**r
        logs.append(Fraction(power, 4*r+1) * z**(4*r+1))
        zetas.append(-Fraction(power, 2*(4*r-1)) * z**(4*r-1))
    return logs, zetas


def approximation_vectors(n1, m1, n2, m2, cutoff):
    d = n2*m1-n1*m2
    assert abs(d) == 1
    values = []
    for n in (n1, n2):
        point, cocycle = mul_with_cocycle(n)
        assert point == mul(n)
        z = -point[0]/point[1]
        values.append((z, cocycle, *local_coefficients(z, cutoff)))
    (z1, r1, l1, h1), (z2, r2, l2, h2) = values
    omega = [(n1*b - n2*a)/d for a, b in zip(l1, l2)]
    quasi = [(n1*b - n2*a)/d for a, b in zip(h1, h2)]
    quasi[0] += (n2*r1-n1*r2)/d
    assert 144*max(abs(z1), abs(z2))**4 < 1
    return omega, quasi


def joint_value(omega, quasi, cutoff):
    return sum(omega[r]*quasi[s]
               for r in range(cutoff+1) for s in range(cutoff+1-r))


def signed_log_error(value, low_pi, high_pi):
    lower, upper = value-low_pi, value-high_pi
    assert lower * upper > 0
    low, high = sorted((abs(lower), abs(upper)))

    def ratio_log(rational):
        return log10(rational.numerator) - log10(rational.denominator)

    at_low = ratio_log(low)
    at_high = ratio_log(high)
    assert f"{at_low:.8f}" == f"{at_high:.8f}"
    return ("+" if lower > 0 else "-"), f"{at_low:.8f}"


def valuation(rational, prime):
    if rational == 0:
        return float("inf")
    num, den = abs(rational.numerator), rational.denominator
    result = 0
    while num % prime == 0:
        num //= prime
        result += 1
    while den % prime == 0:
        den //= prime
        result -= 1
    return result


def main():
    assert mul(1) == P
    assert mul(2) == (Fraction(25, 4), Fraction(35, 8))
    assert [count_points_mod_p(p) for p in (5, 7, 11, 13)] == [8, 8, 12, 20]
    low, high = pi_bounds()
    logs, zetas = approximation_vectors(7, 2, 11, 3, 12)
    assert Fraction(6, 100) < abs(xi(7)) < Fraction(64, 1000)
    assert 0 < xi(11) < Fraction(54, 1000) < abs(xi(7))
    assert valuation(xi(7), 83) == -1 and valuation(xi(11), 83) == 0
    assert all(valuation(mul_with_cocycle(n)[1], 83) >= 0 for n in (7, 11))
    for cutoff in (0, 1, 2, 3, 4, 5, 6, 9, 12):
        direct = sum(logs[:cutoff+1])*sum(zetas[:cutoff+1])
        joint = joint_value(logs, zetas, cutoff)
        for label, value in (("direct", direct), ("joint", joint)):
            sign, log_error = signed_log_error(value, low, high)
            print(f"K={cutoff:2d} {label:6s}: denominator bits="
                  f"{value.denominator.bit_length():5d}, "
                  f"log10|error|={log_error}, sign={sign}")

    # Bound pi directly from exact interval evaluations of both local series.
    # Each of L and Z has alternating terms with decreasing absolute magnitudes.
    def local_interval(z):
        lterms, hterms = local_coefficients(z, 13)
        return [sorted((sum(terms[:13]), sum(terms))) for terms in (lterms, hterms)]

    l7, h7 = local_interval(xi(7))
    l11, h11 = local_interval(xi(11))
    cocycle = 11*mul_with_cocycle(7)[1] - 7*mul_with_cocycle(11)[1]
    omega_low, omega_high = 7*l11[0]-11*l7[1], 7*l11[1]-11*l7[0]
    quasi_low, quasi_high = cocycle+7*h11[0]-11*h7[1], cocycle+7*h11[1]-11*h7[0]
    assert omega_low > 1 and 0 < quasi_low < quasi_high < 4
    assert omega_low*quasi_low < low < high < omega_high*quasi_high
    print("Legendre and independent Machin intervals agree.")

    weights = {2: 1478890, 3: -241928, 4: 15301}
    norm = 2293200
    for moment, target in ((1, Fraction(1)), (5, Fraction(11, 6)),
                           (9, Fraction(39, 14))):
        assert sum(Fraction(weight*j**moment, norm)
                   for j, weight in weights.items()) == target

    def three_point_filter(n):
        return sum(weight*xi(j*n) for j, weight in weights.items()) / (norm*xi(n)**2)

    filtered_quasi = cocycle+7*three_point_filter(11)-11*three_point_filter(7)
    filtered_value = sum(logs[:2])*filtered_quasi
    _, filtered_error = signed_log_error(filtered_value, low, high)
    assert filtered_value.denominator.bit_length() == 7270
    print(f"Three-point zeta filter: 7270 denominator bits, log10|error|={filtered_error}")

    # Exercise the exact local prediction at the first p-power.
    logs, zetas = approximation_vectors(7, 2, 11, 3, 83)
    assert valuation(sum(logs)*sum(zetas), 83) == -8*83
    assert valuation(joint_value(logs, zetas, 83), 83) == -4*83-1
    print("At K=83, prime-83 denominator valuations are 664 and 333.")


if __name__ == "__main__":
    main()
