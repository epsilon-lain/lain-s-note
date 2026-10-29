#!/usr/bin/env python3
"""Exact finite checks for rational decimal words evaluated in base pi.

The arithmetic intervals are rational, using the Machin pi enclosure from
the preceding companion script.  Transcendence and injectivity are proved
in the note, rather than inferred from decimal output.
"""
from fractions import Fraction
from pi_gene_base_2026_10_30 import (
    pi_bounds, decimal_digits, weighted_prefix, format_interval,
)


def decimal_word(value):
    """Canonical nonnegative decimal: integer part, prefix, period.

    Terminating decimals use a finite prefix and an empty period, i.e.
    the zero tail.  No trailing-nine alternate representation is used.
    """
    assert value >= 0
    integer, remainder = divmod(value.numerator, value.denominator)
    digits, seen = [], {}
    while remainder and remainder not in seen:
        seen[remainder] = len(digits)
        digit, remainder = divmod(10*remainder, value.denominator)
        digits.append(digit)
    if remainder == 0:
        return integer, digits, []
    split = seen[remainder]
    return integer, digits[:split], digits[split:]


def word_polynomial(digits, base):
    out = Fraction()
    for digit in digits:
        out = out*base + digit
    return out


def word_value(word, base):
    integer, prefix, cycle = word
    power = base**len(prefix)
    value = Fraction(integer) + word_polynomial(prefix, base)/power
    if cycle:
        value += word_polynomial(cycle, base)/(power*(base**len(cycle)-1))
    return value


def value_interval(value, lower_pi, upper_pi):
    word = decimal_word(value)
    # Nonnegative fractional digits make this rational function decreasing.
    return word_value(word, upper_pi), word_value(word, lower_pi)


def quotient_interval(first, second):
    a, b = first
    c, d = second
    assert 0 < c <= d
    return a/d, b/c


def main():
    lo_pi, hi_pi = pi_bounds()
    for q in range(1, 41):
        for p in range(0, min(4*q, 15)):
            value = Fraction(p, q)
            assert word_value(decimal_word(value), Fraction(10)) == value

    samples = (Fraction(1, 3), Fraction(1, 7), Fraction(22, 7), Fraction(355, 113))
    intervals = {}
    for value in samples:
        word = decimal_word(value)
        lo, hi = value_interval(value, lo_pi, hi_pi)
        assert lo < hi
        intervals[value] = (lo, hi)
        print(f'{value}: prefix {word[1]}, period length {len(word[2])}')
        print('  T_pi =', format_interval(lo, hi, 25) + '...')
    ratio = quotient_interval(intervals[Fraction(355,113)], intervals[Fraction(22,7)])
    print('T_pi(355/113)/T_pi(22/7) =', format_interval(*ratio, 25) + '...')

    digits = decimal_digits(lo_pi, hi_pi, 30)
    for n in (1, 2, 3, 6, 12, 20):
        initial = digits[:n]
        next_digit = digits[n]
        # Directly certified images of two successive decimal truncations.
        a_lo = weighted_prefix(initial, hi_pi)
        a_hi = weighted_prefix(initial, lo_pi)
        b_lo = weighted_prefix(digits[:n+1], hi_pi)
        b_hi = weighted_prefix(digits[:n+1], lo_pi)
        ratio_lo, ratio_hi = quotient_interval((b_lo,b_hi), (a_lo,a_hi))
        assert ratio_hi-1 <= 3/lo_pi**(n+1)
        if next_digit:
            assert ratio_lo > 1
        print(f'n={n}: next digit {next_digit}; A_(n+1)/A_n =',
              format_interval(ratio_lo, ratio_hi, 20) + '...')

    # Rational-source order is not preserved at b=pi (<10).
    first = value_interval(Fraction(19,100), lo_pi, hi_pi)
    second = value_interval(Fraction(1,5), lo_pi, hi_pi)
    assert Fraction(19,100) < Fraction(1,5) and first[0] > second[1]

    # Two transcendental images can have an exactly rational quotient when
    # their original integer parts vanish.
    at_tenth = word_value(decimal_word(Fraction(1,10)), Fraction(7,2))
    at_fifth = word_value(decimal_word(Fraction(1,5)), Fraction(7,2))
    assert at_fifth == 2*at_tenth
    print('All exact finite checks passed. Infinite claims rely on the written proofs.')


if __name__ == '__main__':
    main()
