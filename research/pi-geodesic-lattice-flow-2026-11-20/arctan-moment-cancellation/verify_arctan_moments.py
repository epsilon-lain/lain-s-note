"""Check fixed-prime obstruction and exact denominator cancellation."""
from fractions import Fraction
from math import isqrt


def truncated_atan(argument, terms):
    return sum(
        (-1)**index * argument**(2 * index + 1) / (2 * index + 1)
        for index in range(terms)
    )


def primes_below(bound):
    for value in range(2, bound):
        if all(value % d for d in range(2, isqrt(value) + 1)):
            yield value


def dyadic_valuation(value):
    return (value & -value).bit_length() - 1


def main():
    left_weight = Fraction(-5706, 269)
    right_weight = Fraction(5975, 269)
    first_moment_left = 4 * (Fraction(1, 2) + Fraction(1, 3))
    first_moment_right = Fraction(16, 5) - Fraction(4, 239)
    assert left_weight + right_weight == 1
    assert left_weight * first_moment_left + right_weight * first_moment_right == 0
    # This gives right arm <= half the main left arm for every d >= 3.
    assert 20 * right_weight * Fraction(2, 5)**5 < 3 * abs(left_weight) / 2
    removed_cases = 0
    for n in (2, 10, 120, 404, 600):
        degree = 2 * n - 1
        left = 4 * (truncated_atan(Fraction(1, 2), n) + truncated_atan(Fraction(1, 3), n))
        right = 16 * truncated_atan(Fraction(1, 5), n) - 4 * truncated_atan(Fraction(1, 239), n)
        mixed = left_weight * left + right_weight * right
        assert dyadic_valuation(mixed.denominator) == degree - 3
        removed = [
            p for p in primes_below(degree + 1)
            if p > 269 and 3 * p > degree and p * p > degree
        ]
        assert all(left.denominator % p == 0 for p in removed)
        assert all((right.denominator % p != 0) == (first_moment_right.numerator % p == 0)
                   for p in removed)
        assert all(mixed.denominator % p != 0 for p in removed)
        removed_cases += len(removed)
        print(f"N={n}: {len(removed)} primes removed, denominator 2-adic exponent {degree-3}")
    print(f"Checked {removed_cases} prime removals across all listed truncations.")


if __name__ == "__main__":
    main()
