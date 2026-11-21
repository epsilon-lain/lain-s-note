"""Exact finite endpoints for the fixed 1/2, 1/3 Machin truncation theorem.

All divisions are rational, including the independent 333/106 lower bound
for pi. The infinite prime-product estimate is proved in the note.
"""
from fractions import Fraction
from math import isqrt


def atan_interval(base, terms):
    partial = sum((Fraction((-1) ** k, (2*k+1)*base**(2*k+1))
                   for k in range(terms)), Fraction())
    next_term = Fraction((-1)**terms, (2*terms+1)*base**(2*terms+1))
    return min(partial, partial+next_term), max(partial, partial+next_term)


def main():
    a, b = atan_interval(5, 8)
    c, f = atan_interval(239, 8)
    assert 16*a-4*f > Fraction(333, 106) > Fraction(505, 161)
    value = Fraction(0)
    prime_checks = 0
    checked_small = 0
    for n in range(1, 129):
        d = 2*n-1
        value += Fraction(4*(-1)**(n-1), d) * (
            Fraction(1, 2**d) + Fraction(1, 3**d))
        q, p = value.numerator, value.denominator
        if n == 2:
            assert (p, q) == (162, 505)
        if 3 <= n <= 48:
            assert 3*p > 4*(d+2)*2**(d+2)
            checked_small += 1
        for prime in range(max(7, d//2+1), d+1):
            if any(prime % v == 0 for v in range(2, isqrt(prime)+1)):
                continue
            scaled = prime*value
            expected = Fraction(10*(-1)**((prime-1)//2), 3)
            assert scaled.denominator % prime
            assert (scaled-expected).numerator % prime == 0
            assert p % prime == 0 and p % (prime*prime) != 0
            prime_checks += 1
    assert 3*2**21 > 64*99
    print(f'N=2 certified; {checked_small} further small endpoints; '
          f'{prime_checks} independent local prime checks through N=128.')


if __name__ == '__main__':
    main()
