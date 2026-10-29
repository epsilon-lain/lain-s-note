"""Certified audit of decimal periods and pi approximants.

Run from this directory with Python's standard library. The imported pi interval
uses exact Machin alternating-series bounds. The general argument in the
companion note is not a Lean-formalized proof of any statement about pi.
"""

from fractions import Fraction as F
from math import gcd, lcm

from dyadic_rectangles_2026_10_29 import pi_bounds


def factors(value):
    n = value
    output = []
    d = 2
    while d*d <= n:
        if n % d == 0:
            exponent = 0
            while n % d == 0:
                exponent += 1
                n //= d
            output.append((d, exponent))
        d += 1 if d == 2 else 2
    if n > 1:
        output.append((n, 1))
    assert value == product(p**e for p, e in output)
    return output


def product(values):
    result = 1
    for v in values:
        result *= v
    return result


def order_10(q):
    """Least period and phi/lambda for the component coprime to 10."""
    rest = q
    while rest % 2 == 0:
        rest //= 2
    while rest % 5 == 0:
        rest //= 5
    if rest == 1:
        return rest, 0, 1, 1
    fact = factors(rest)
    phi = product((p-1)*p**(e-1) for p,e in fact)
    exponent = lcm(*[(p-1)*p**(e-1) for p,e in fact])
    order = exponent
    for prime, _ in factors(exponent):
        while order % prime == 0 and pow(10, order // prime, rest) == 1:
            order //= prime
    assert pow(10, order, rest) == 1
    assert all(pow(10, order // prime, rest) != 1 for prime,_ in factors(order))
    return rest, order, phi, exponent


def scaled_error(lo, hi, p, q):
    value = F(p, q)
    assert not lo <= value <= hi
    if value < lo:
        return q*q*(lo-value), q*q*(hi-value)
    return q*q*(value-hi), q*q*(value-lo)


def convergents(lo, hi, limit):
    left, right = lo, hi
    a_old, a_new, b_old, b_new = 0, 1, 1, 0
    output = []
    while True:
        digit = left.__floor__()
        assert digit == right.__floor__()
        numerator, denominator = digit*a_new+a_old, digit*b_new+b_old
        a_old, a_new, b_old, b_new = a_new, numerator, b_new, denominator
        if denominator > limit:
            return output
        assert gcd(numerator, denominator) == 1
        output.append((numerator, denominator, digit))
        left, right = 1/(right-digit), 1/(left-digit)


def divisors(n):
    output = [1]
    for prime, exponent in factors(n):
        output = [a*prime**i for a in output for i in range(exponent+1)]
    assert all(n % x == 0 for x in output)
    return sorted(output)


def best_divisor_cycle(lo, hi, length):
    candidates = []
    for q in divisors(10**length - 1):
        lower_candidate = (q*lo + F(1, 2)).__floor__()
        upper_candidate = (q*hi + F(1, 2)).__floor__()
        assert lower_candidate == upper_candidate
        p = lower_candidate
        low, high = scaled_error(lo, hi, p, q)
        candidates.append((low, high, p, q))
    candidates.sort(key=lambda v: v[0])
    assert candidates[0][1] < candidates[1][0]
    return candidates[0], len(candidates)


def construct_cf_counterexample_prefix():
    # a1=a2=2, a_n=n-1 for odd n>=3, a_n=5(n/2-1) for even n>=4.
    a_old, a_new, b_old, b_new = 1, 0, 0, 1
    for n in range(1, 41):
        digit = 2 if n <= 2 else n-1 if n % 2 else 5*(n//2-1)
        numerator, denominator = digit*a_new+a_old, digit*b_new+b_old
        a_old, a_new, b_old, b_new = a_new, numerator, b_new, denominator
        assert denominator % 10 == (2 if n % 2 else 5)
    return True


def main():
    lo, hi = pi_bounds()
    convergent_list = convergents(lo, hi, 7_000_000_000)
    convergent_set = {(p,q) for p,q,_ in convergent_list}
    exceptionally_good = []
    for p,q,_ in convergent_list:
        low, high = scaled_error(lo, hi, p, q)
        assert high < F(1,10) or low > F(1,10)
        if high < F(1,10):
            exceptionally_good.append((p,q))
    assert exceptionally_good == [
        (22,7), (355,113), (5419351,1725033), (21053343141,6701487259),
    ]
    print('All convergents with q<=7e9 and score<0.1:', exceptionally_good)
    samples = [
        (22, 7), (53, 17), (355, 113), (1146408, 364913),
        (5419351, 1725033), (21053343141, 6701487259),
        (312689, 99532),
    ]
    for p,q in samples:
        assert gcd(p,q) == 1
        rest, order, phi, exponent = order_10(q)
        low, high = scaled_error(lo, hi, p, q)
        print(f'{p}/{q}: score≈{float((low+high)/2):.12g}; q0={rest}; period={order}; phi={phi}; lambda={exponent}')
    assert (21053343141, 6701487259) in convergent_set
    assert (312689, 99532) in convergent_set
    assert factors(6701487259) == [(101,1), (4457,1), (14887,1)]
    assert order_10(6701487259) == (6701487259, 33166008, 6633201600, 829150200)
    assert order_10(1725033) == (1725033, 31824, 1145664, 31824)
    assert order_10(99532)[1] == 12284
    assert pow(10, 96, 113) == 16 and pow(10, 10, 17) == 2
    assert 3 + sum((F(1,8**j) for j in range(1,12)),F()) < F(22,7)
    assert 3 + F(1,8-1) == F(22,7)
    assert 3 + F(16,114-1) == F(355,113)
    windows = [(2, 99), (100, 999), (1724033, 1726033), (6701486259, 6701488259)]
    expected_counts = [(16,39), (137,360), (330,801), (295,801)]
    for (start,end), expected in zip(windows, expected_counts):
        candidate_q = [q for q in range(start,end+1) if gcd(q,10)==1]
        maximum = sum(order_10(q)[1] == order_10(q)[3] for q in candidate_q)
        assert (maximum,len(candidate_q)) == expected
        print(f'maximal order q={start}..{end}: {maximum}/{len(candidate_q)}')
    for length in (6,12):
        best, count = best_divisor_cycle(lo, hi, length)
        low, high, p, q = best
        assert (p,q) == (22,7)
        print(f'period dividing {length}: {count} divisors; optimum {p}/{q}; score≈{float((low+high)/2):.12g}')
    assert construct_cf_counterexample_prefix()
    print('Certified samples, minimal-period certificates, divisor scans, and 40 counterexample recurrences passed.')


if __name__ == '__main__':
    main()
