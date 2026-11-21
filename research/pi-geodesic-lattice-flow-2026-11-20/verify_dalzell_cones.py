"""Exact certificates for the two Dalzell prime-band estimates. Python 3."""
from fractions import Fraction as F
from math import comb, isqrt


def log_interval(x):
    """Enclose log(x) for a positive rational by the atanh series."""
    x = F(x)
    z = (x - 1) / (x + 1)
    terms = 80
    value = 2 * sum(z**(2*i + 1) / (2*i + 1) for i in range(terms))
    tail = 2 * abs(z)**(2*terms + 1) / ((2*terms + 1) * (1-z*z))
    return value - tail, value + tail


def coefficient(n, reflected):
    """Coefficient of x**(n-1) in the corresponding reciprocal quadratic."""
    if reflected:
        sign = (0, 1, 1, 1, 0, -1, -1, -1)[n % 8]
        return F(sign, 2**((n+1)//2))
    return F((0, 1, 0, -1)[n % 4])


def harmonic(j, reflected):
    return sum(coefficient(b, reflected) / b for b in range(1, j+1))


def verify_budget(h, orders, reflected, expected_width, margin, max_sum):
    bounds = [(4*(1+h)/(j+1), F(4, j)) for j in range(1, orders+1)]
    assert all(lo < hi and lo >= F(44, 100) for lo, hi in bounds)
    assert all(bounds[j][1] < bounds[j-1][0] for j in range(1, orders))
    width = sum(hi-lo for lo, hi in bounds)
    endpoint_sum = sum(hi+lo for lo, hi in bounds)
    assert width == expected_width and endpoint_sum < max_sum
    entropy_upper = 4 * ((1+h)*log_interval(1+h)[1] - h*log_interval(h)[0])
    cost_upper = entropy_upper + (1 if reflected else h)*log_interval(4)[1]
    assert width-cost_upper > margin
    assert margin - F(max_sum, 1000) - F(1, 1000) > F(1, 200)
    for j in range(1, orders+1):
        val = harmonic(j, reflected)
        assert 0 < val.numerator <= 10080 and val.denominator <= 10080


def exact_integral(k, ell):
    row = [0] * (4*(k+ell)+1)
    for j in range(4*ell+1):
        row[4*k+j] = (-1)**j * comb(4*ell, j)
    result = F(0)
    for degree in range(len(row)-1, 1, -1):
        value = row[degree]
        row[degree-2] -= value
        result += F(value, degree-1)
    assert row[:2] == [(-4)**ell, 0]
    return result


def verify_prime_bands():
    pairs = [(30, 2), (75, 11), (90, 13), (50, 9), (37, 3), (2, 30),
             (5, 51), (9, 90), (6, 72), (3, 37), (1, 121), (122, 1)]
    counts = [0, 0]
    for k, ell in pairs:
        result = exact_integral(k, ell)
        residue = (-4)**ell
        reduced_denominator = (-4*result/residue).denominator
        for reflected, scale in [(False, k), (True, ell)]:
            for j in range(1, 12):
                value = harmonic(j, reflected)
                for p in range(3, 4*scale//j+1, 2):
                    if any(p % i == 0 for i in range(2, isqrt(p)+1)):
                        continue
                    if p <= j or value.numerator % p == 0:
                        continue
                    if (j+1)*p <= 4*(k+ell) or j*p >= 4*scale:
                        continue
                    character = 1 if p % 4 == 1 else -1
                    diff = p*result + residue*character*value
                    assert diff.denominator % p and diff.numerator % p == 0
                    assert reduced_denominator % p == 0
                    assert reduced_denominator % (p*p) != 0
                    counts[reflected] += 1
    assert counts == [217, 235]
    return counts


if __name__ == "__main__":
    verify_budget(F(9, 50), 5, False, F(1717, 750), F(23, 1000), 16)
    verify_budget(F(1, 10), 9, True, F(17819, 6300), F(101, 1000), 20)
    assert 4400000 > 2**22 and (22*log_interval(2)[0])**2 > 200
    assert log_interval(10)[1] < 3
    print("Prime congruences (direct, reflected):", *verify_prime_bands())
    print("Both entropy margins and explicit thresholds verified.")
