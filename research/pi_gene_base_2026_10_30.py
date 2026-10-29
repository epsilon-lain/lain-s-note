#!/usr/bin/env python3
"""Certified numeric examples for reweighting decimal digits of pi.

No third-party dependencies.  Uses Machin's identity and alternating-series
rational enclosures for pi.  Every printed digit is certified by an interval;
numeric examples do not prove any infinite repetition condition.
"""
from decimal import Decimal, localcontext
from fractions import Fraction


def atan_bounds(inv, terms):
    total = sum((Fraction((-1)**j, (2*j+1)*inv**(2*j+1))
                 for j in range(terms)), Fraction())
    err = Fraction(1, (2*terms+1)*inv**(2*terms+1))
    return (total, total+err) if terms % 2 == 0 else (total-err, total)


def pi_bounds():
    a, b = atan_bounds(5, 80)
    c, d = atan_bounds(239, 30)
    return 16*a-4*d, 16*b-4*c


def dec_floor(x):
    return x.numerator // x.denominator


def decimal_digits(lo, hi, n):
    scale = 10**n
    x, y = dec_floor(lo*scale), dec_floor(hi*scale)
    assert x == y
    assert str(x)[0] == '3'
    return [int(c) for c in str(x)[1:]]


def weighted_prefix(digits, base):
    acc = Fraction(3)
    weight = Fraction(1)
    for digit in digits:
        weight /= base
        acc += digit * weight
    return acc


def gene_interval(digits, base_lo, base_hi):
    assert 1 < base_lo <= base_hi
    n = len(digits)
    lo = weighted_prefix(digits, base_hi)
    hi = weighted_prefix(digits, base_lo) + 9 / (base_lo-1) / base_lo**n
    assert lo < hi
    return lo, hi


def format_interval(lo, hi, places):
    scale = 10**places
    lower, upper = dec_floor(lo*scale), dec_floor(hi*scale)
    assert lower == upper, 'Increase digit count to certify requested places'
    s = str(lower).zfill(places+1)
    return s[:-places] + '.' + s[-places:]


def outward_decimal_bounds(lo, hi, places):
    scale = 10**places
    lower = dec_floor(lo*scale)
    upper = -dec_floor(-hi*scale)
    def render(v):
        s = str(v).zfill(places+1)
        return s[:-places] + '.' + s[-places:]
    return render(lower), render(upper)


def greedy_from_interval(lo, hi, base, count):
    lo, hi = lo-3, hi-3
    out = []
    for _ in range(count):
        lo, hi = lo*base, hi*base
        digit = dec_floor(lo)
        assert digit == dec_floor(hi)
        out.append(digit)
        lo, hi = lo-digit, hi-digit
    return out


def pure_cycle_in_integer_base(integer_part, digits, base):
    value = 0
    for digit in digits:
        assert 0 <= digit < base
        value = value*base + digit
    return Fraction(integer_part) + Fraction(value, base**len(digits)-1)


def main():
    pi_lo, pi_hi = pi_bounds()
    assert pi_hi-pi_lo < Fraction(1, 10**110)
    digits = decimal_digits(pi_lo, pi_hi, 60)
    print('Certified pi decimal digits:', '3.' + ''.join(map(str, digits)))
    examples = [('pi', pi_lo, pi_hi), ('4', Fraction(4), Fraction(4)),
                ('10', Fraction(10), Fraction(10)), ('11', Fraction(11), Fraction(11)),
                ('12', Fraction(12), Fraction(12)), ('16', Fraction(16), Fraction(16))]
    prev = None
    intervals = {}
    for label, lo_b, hi_b in examples:
        lo, hi = gene_interval(digits, lo_b, hi_b)
        if label == '10':
            assert lo < pi_lo < pi_hi < hi
        if label in ('11', '12', '16'):
            assert greedy_from_interval(lo, hi, lo_b, 35) == digits[:35]
        if prev is not None:
            assert hi < prev
        prev = lo
        intervals[label] = (lo, hi)
        print(f'G({label}) = {format_interval(lo, hi, 25)}...')
        if label == 'pi':
            custom_lo = 3*pi_lo + lo-3
            custom_hi = 3*pi_hi + hi-3
            print('H(pi) = 3*pi + (G(pi)-3) =',
                  format_interval(custom_lo, custom_hi, 25) + '...')
    a, b = intervals['pi']
    print('G(pi)/pi =', format_interval(a/pi_hi, b/pi_lo, 25) + '...')
    a, b = intervals['11']
    print('G(11)/pi =', format_interval(a/pi_hi, b/pi_lo, 25) + '...')
    a, b = intervals['12']
    print('G(12)/pi =', format_interval(a/pi_hi, b/pi_lo, 25) + '...')
    literal_lo, literal_hi = gene_interval([1,4,1,5,9,3,6], pi_lo, pi_hi)
    print('Literal prefix 3.1415936...: G(pi) is in',
          outward_decimal_bounds(literal_lo, literal_hi, 10))
    print('Literal prefix 3.1415936...: H(pi) is in',
          outward_decimal_bounds(3*pi_lo+literal_lo-3, 3*pi_hi+literal_hi-3, 10))
    rational_11 = pure_cycle_in_integer_base(3, [1, 4, 2, 8, 5, 7], 11)
    lo_11, hi_11 = intervals['11']
    assert rational_11 > hi_11
    q = rational_11.denominator
    metric_lo, metric_hi = q*q*(rational_11-hi_11), q*q*(rational_11-lo_11)
    print(f'22/7 decimal cycle in base 11: {rational_11}; reduced q={q}')
    print('q^2 |G(11)-cycle_11(22/7)| =', format_interval(metric_lo, metric_hi, 6) + '...')
    print('Recovered first 35 decimal-pi digits from G(11), G(12), G(16).')
    print('All interval checks passed.')


if __name__ == '__main__':
    main()
