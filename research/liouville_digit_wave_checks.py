"""Exact finite checks of the Liouville base-wave certificate."""
from fractions import Fraction as F
from math import factorial


def main():
    for k in range(2, 6):
        q = 10**factorial(k)
        p = q+1
        x0 = sum((F(1, 10**factorial(j)) for j in range(1, k+1)), F())
        n = int(q*x0)
        # Bounds for L - x0, obtained from the factorial gaps.
        low = x0+F(1, q**(k+1))
        high = x0+F(11, 10*q**(k+1))
        digits = []
        for _ in range(k+1):
            first, last = p*low, p*high
            dl = first.numerator//first.denominator
            du = last.numerator//last.denominator
            assert dl == du
            digits.append(dl)
            low, high = first-dl, last-du
        assert digits == [n]*k+[n+1]
        assert 2*q*q < p**(k+1)
        print(f'k={k}: k equal digits then +1; certified E < 1/q²')

    # For any q <= 1000, a/q != 11/100 is >= 1/(100q) away.
    # L - 11/100 lies strictly between 10^-6 and 2*10^-6.
    assert F(1, 100*1000)-F(2, 10**6) > F(2, 10**6)
    print('q <= 1000: certified distance > 1/1000²')


if __name__ == '__main__':
    main()
