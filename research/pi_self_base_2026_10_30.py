#!/usr/bin/env python3
"""Exact interval checks for pi-base expansions and nearby rational bases.

Uses only Python's standard library.  The pi interval follows from Machin's
identity and alternating rational arctangent sums.  Decimal values printed
below are for readability; every digit and sign is established by Fractions.
"""

from decimal import Decimal, localcontext
from fractions import Fraction


def atan_bounds(inv, terms):
    s = sum((Fraction((-1) ** j, (2*j+1)*inv**(2*j+1))
             for j in range(terms)), Fraction())
    e = Fraction(1, (2*terms+1)*inv**(2*terms+1))
    return (s, s+e) if terms % 2 == 0 else (s-e, s)


def pi_bounds():
    a, b = atan_bounds(5, 90)
    c, d = atan_bounds(239, 35)
    return 16*a-4*d, 16*b-4*c


def fmt(x, places=14):
    with localcontext() as ctx:
        ctx.prec = places + 12
        y = Decimal(x.numerator) / Decimal(x.denominator)
        return f"{y:.{places}g}"


def digit_intervals(base_lo, base_hi, rem_lo, rem_hi, count):
    digits = []
    for _ in range(count):
        # The expansion of 1 starts at the endpoint; later remainders are < 1.
        assert 0 <= rem_lo <= rem_hi <= 1
        low, high = base_lo*rem_lo, base_hi*rem_hi
        a = low.numerator // low.denominator
        b = high.numerator // high.denominator
        if a != b:
            raise ArithmeticError("interval cannot certify next digit")
        digits.append(a)
        rem_lo, rem_hi = low-a, high-a
    return digits, rem_lo, rem_hi


def prefix_value(digits, base):
    return sum((Fraction(d, 1) / base**j for j, d in enumerate(digits, 1)), Fraction())


def algebraic_root(digits, bits=95):
    lo, hi = Fraction(3), Fraction(22, 7)
    assert prefix_value(digits, lo) >= 1
    assert prefix_value(digits, hi) < 1
    if prefix_value(digits, lo) == 1:
        return lo, lo
    for _ in range(bits):
        m = (lo+hi)/2
        if prefix_value(digits, m) >= 1:
            lo = m
        else:
            hi = m
    return lo, hi


def bubble_edge(k, side, lower_pi, upper_pi, bits=90):
    if side == '-':
        lo, hi = Fraction(5, 2), Fraction(22, 7)
    else:
        lo, hi = Fraction(3), Fraction(4)
    for _ in range(bits):
        mid = (lo+hi)/2
        if side == '-':
            f_lo, f_hi = lower_pi-mid-mid**(-k), upper_pi-mid-mid**(-k)
            # f strictly decreases here.
            if f_lo > 0:
                lo = mid
            elif f_hi < 0:
                hi = mid
            else:
                raise ArithmeticError('pi interval cannot decide edge')
        else:
            f_lo, f_hi = mid-upper_pi-mid**(-k), mid-lower_pi-mid**(-k)
            # f strictly increases here.
            if f_hi < 0:
                lo = mid
            elif f_lo > 0:
                hi = mid
            else:
                raise ArithmeticError('pi interval cannot decide edge')
    return lo, hi


def one_sided_counterexample_checks():
    # alpha = [3; 1, 2, 1, 3, 1, 4, ...].  Two consecutive
    # convergents rigorously bracket alpha; the limiting construction in the
    # note uses the usual convergent and Legendre inequalities for all n.
    a = [3, 1] + [k + 1 if j % 2 == 0 else 1
                  for j in range(2, 82) for k in [j // 2]]
    prev_p, p, prev_q, q = 0, 1, 1, 0
    convergents = []
    for term in a:
        prev_p, p = p, term*p + prev_p
        prev_q, q = q, term*q + prev_q
        convergents.append(Fraction(p, q))
    alpha_lo, alpha_hi = convergents[80], convergents[81]
    assert alpha_lo < alpha_hi
    for n in range(1, 41):
        approximation = convergents[n]
        q = approximation.denominator
        if n % 2 == 0:
            assert approximation < alpha_lo
            assert q*q*(alpha_lo-approximation) > Fraction(1, 3)
        else:
            assert approximation > alpha_hi
            assert q*q*(approximation-alpha_lo) < Fraction(1, a[n+1])


def main():
    pl, ph = pi_bounds()
    assert ph-pl < Fraction(1, 10**100)
    one_digits, _, _ = digit_intervals(pl, ph, Fraction(1), Fraction(1), 30)
    print('Digits of 1 in base pi (30):', ''.join(map(str, one_digits)))
    print('Non-greedy digits of pi in base pi:', str(one_digits[0]) + '.' + ''.join(map(str, one_digits[1:])))
    print('Algebraic roots of 1 = sum d_j/b^j:')
    for n in (1, 2, 3, 4, 6, 10, 15, 20):
        lo, hi = algebraic_root(one_digits[:n])
        assert hi < pl
        bound = ph / 3 if n == 1 else pl**(2-n) / 3
        assert ph-lo < bound
        print(f' n={n:2}  root={fmt((lo+hi)/2)}  pi-root={fmt((pl+ph-lo-hi)/2, 6)}')
    print('Nearby rational bases: b | sign pi-b | initial zero count | digits of |pi-b| | q^2 |pi-b|')
    for base in (Fraction(3), Fraction(157,50), Fraction(22,7), Fraction(333,106), Fraction(355,113), Fraction(103993,33102)):
        if base < pl:
            dl, dh, sign = pl-base, ph-base, '+'
        elif base > ph:
            dl, dh, sign = base-ph, base-pl, '-'
        else:
            raise ArithmeticError('uncertified sign')
        digits, _, _ = digit_intervals(base, base, dl, dh, 25)
        count = next(i for i,d in enumerate(digits) if d != 0)
        assert base**(-count-1) <= dl <= dh < base**(-count)
        print(f' {str(base):>12} | {sign} | {count:2} | {"".join(map(str,digits))} | {fmt((dl+dh)*base.denominator**2/2, 10)}')
        if sign == '-':
            original, _, _ = digit_intervals(base, base, pl-3, ph-3, 20)
            print(f'                actual pi expansion: (3.{"".join(map(str, original))}...)_{base}')
    print('Symmetric zero-gap bubble edges (certified brackets):')
    for k in (2, 5, 8, 13):
        l0,l1 = bubble_edge(k, '-', pl, ph)
        u0,u1 = bubble_edge(k, '+', pl, ph)
        print(f' k={k:2} | {fmt((l0+l1)/2, 15)} | {fmt((u0+u1)/2, 15)} | width~{fmt((u0+u1-l0-l1)/2, 8)}')
    one_sided_counterexample_checks()
    print('One-sided counterexample: first 40 convergents checked exactly.')
    print('All interval checks passed.')


if __name__ == '__main__':
    main()
