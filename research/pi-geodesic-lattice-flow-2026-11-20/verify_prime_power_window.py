#!/usr/bin/env python3
"""Check the enlarged odd-prime-power window with integer arithmetic."""

from math import gcd


def valuation(value, prime):
    if not value:
        return float('inf')
    exponent = 0
    while value % prime == 0:
        exponent += 1
        value //= prime
    return exponent


def factorial_valuation(bound, prime):
    exponent = 0
    while bound:
        bound //= prime
        exponent += bound
    return exponent


def check_even_numerators():
    cases = 0
    for p in (2, 4, 6, 8, 10, 12, 18, 20, 42, 64):
        for q in (1, 3, 5, 7, 9, 11, 13):
            if gcd(p, q) != 1:
                continue
            prev_real, prev_imag = 1, 0
            real, imag = 2 * q, p
            for n in range(1, 61):
                assert valuation(gcd(real, imag), 2) == n, (p, q, n)
                cases += 1
                factor = (4 * n + 2) * q
                prev_real, real = real, factor * real - p*p * prev_real
                prev_imag, imag = imag, factor * imag - p*p * prev_imag
    return cases


def main():
    cases = beyond_original = 0
    for ell in (3, 5, 7, 11, 13, 17, 19, 23, 29):
        multipliers = (1, 2, 3 if ell == 7 else 7)
        for power in (1, 2):
            for multiplier in multipliers:
                p = multiplier * ell**power
                for q in (1, 2, 5, 11):
                    if gcd(p, q) != 1:
                        continue
                    final_order = min((3 * ell**(2 * power) - 1) // 2, 600)
                    prev_real, prev_imag = 1, 0
                    real, imag = 2 * q, p
                    for n in range(1, final_order + 1):
                        if n % 2 == 0 and 2 * n - 1 < 3 * ell**(2 * power):
                            baseline = factorial_valuation(2*n, ell) - factorial_valuation(n, ell)
                            assert valuation(real, ell) == baseline, (p, q, ell, n)
                            assert valuation(imag, ell) == baseline + power, (p, q, ell, n)
                            assert valuation(gcd(real, imag), ell) == baseline, (p, q, ell, n)
                            cases += 1
                            beyond_original += 2*n - 1 >= 3*ell**2
                        factor = (4*n + 2) * q
                        prev_real, real = real, factor * real - p*p * prev_real
                        prev_imag, imag = imag, factor * imag - p*p * prev_imag
    print(f'{cases} cases tested; {beyond_original} beyond the original window')
    print(f'{check_even_numerators()} even-numerator cases tested at prime 2')


if __name__ == '__main__':
    main()
