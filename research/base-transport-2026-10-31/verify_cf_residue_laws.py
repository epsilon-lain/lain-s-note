"""Exact certificates for finite-depth CF residue laws under base transport.

Uses only the Python standard library. Infinite claims follow from the
displayed polynomial positivity certificates, not from sampled bases.
"""

import importlib.util
import json
from fractions import Fraction as F
from math import comb
from pathlib import Path


HELPER = Path(__file__).resolve().parents[1] / 'pi-native-space-2026-10-31' / 'verify_native.py'
spec = importlib.util.spec_from_file_location('native_certificate', HELPER)
native = importlib.util.module_from_spec(spec)
spec.loader.exec_module(native)


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(p, q):
    r = [0] * max(len(p), len(q))
    for i, v in enumerate(p):
        r[i] += v
    for i, v in enumerate(q):
        r[i] += v
    return trim(r)


def neg(p):
    return [-v for v in p]


def mul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, v in enumerate(p):
        for j, w in enumerate(q):
            r[i+j] += v*w
    return trim(r)


def shift(p, offset):
    r = [0] * len(p)
    for n, v in enumerate(p):
        for k in range(n+1):
            r[k] += v * comb(n, k) * offset**(n-k)
    return trim(r)


def difference_numerator(left, right):
    p, q = left
    r, s = right
    return add(mul(p, s), neg(mul(r, q)))


def positive_from(p, threshold):
    q = shift(p, threshold)
    assert q[0] > 0 and all(v >= 0 for v in q)
    return q


def certified_cf(low, high, length):
    result = []
    for _ in range(length):
        digit = low.numerator // low.denominator
        assert digit == high.numerator // high.denominator
        result.append(digit)
        low, high = 1/(high-digit), 1/(low-digit)
    return result


def main():
    pi_low, pi_high = native.pi_bounds(100)
    denominator = 2**500
    pi_low = F((pi_low*denominator).__floor__(), denominator)
    pi_high = F((pi_high*denominator).__ceil__(), denominator)
    digits = native.certified_digits(F(1), F(1), pi_low, pi_high, 64)
    assert digits[:9] == [3, 0, 1, 1, 0, 2, 1, 1, 1]

    # Polynomial coefficients below are listed in ascending order.
    b8 = [0]*8+[1]
    bm1 = [-1, 1]
    ln = [1, 1, 1, 2, 0, 1, 1]
    lower = (ln, b8)
    upper = (add(mul(ln, bm1), [3]), mul(b8, bm1))
    w0 = ([4, 3], [-5, -1, 1, 3])
    w1 = ([21, 16], [-27, -5, 5, 16])
    w2 = ([8, 9, 2], [-11, -5, 1, 7, 2])
    w3 = ([6, 6, 1], [-8, -3, 1, 5, 1])
    cases = [
        ('L-W0', lower, w0, [-5, -6, -5, -7, 2], 5),
        ('W1-U', w1, upper, [54, 10, -10, -59, 49, -12, 1], 5),
        ('W2-U', w2, upper, [22, 10, -2, -25, 13], 14),
        ('L-W3', lower, w3, [-8, -11, -10, -13, 1], 14),
    ]
    positivity = []
    for label, left, right, expected, threshold in cases:
        assert difference_numerator(left, right) == expected
        shifted = positive_from(expected, threshold)
        positive_from(left[1], threshold)
        positive_from(right[1], threshold)
        positivity.append({
            'difference': label,
            'threshold': threshold,
            'numerator_coefficients_ascending': expected,
            'shifted_positive_coefficients_ascending': shifted,
        })

    prefix_counts = {str(r): 0 for r in range(3)}
    next_counts = {str(r): 0 for r in range(3)}
    samples = []
    bases = list(range(5, 101)) + [1000]
    for b in bases:
        low = sum((F(d, b**i) for i, d in enumerate(digits)), F(0))
        high = low + F(3, b**(len(digits)-1)*(b-1))
        actual = certified_cf(low, high, 8)
        r = b % 3
        a = b*b-b
        if r == 0:
            predicted = [3, a, 1, (b-3)//3, 2, 3, 1]
            next_predicted = (b-2)//6
        elif r == 1:
            predicted = [3, a, 1, (b-4)//3, 1, 3, 2]
            next_predicted = b//6
        else:
            predicted = [3, a, 1, (b-2)//3, 9]
            next_predicted = (b+1)//6
        assert actual[:len(predicted)] == predicted
        prefix_counts[str(r)] += 1
        if b >= 14:
            assert actual[len(predicted)] == next_predicted
            next_counts[str(r)] += 1
        if b in [5, 6, 7, 10, 14, 15, 16, 1000]:
            samples.append({'base': b, 'certified_cf_first_8': actual})

    return {
        'native_digit_count_certified': len(digits),
        'native_prefix_18': ''.join(map(str, digits[:18])),
        'pi_enclosure': 'Machin alternating sums, 100 terms; outward dyadic enclosure at 500 bits',
        'polynomial_certificates': positivity,
        'sampled_bases': '5 through 100, and 1000',
        'sampled_base_count': len(bases),
        'mod3_prefix_counts': prefix_counts,
        'mod6_next_partial_quotient_counts': next_counts,
        'samples': samples,
        'scope': 'Exact finite checks plus symbolic polynomial positivity certificates; no P2848 conclusion.',
    }


if __name__ == '__main__':
    print(json.dumps(main(), ensure_ascii=False, indent=2))
