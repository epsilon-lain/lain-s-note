"""Exact finite checks; see PRIME_SAMPLING_AND_GENERATORS.md for proofs."""
from fractions import Fraction as F
from math import gcd, isqrt
from pathlib import Path
import json
from verify_transport import certify_digits, evaluate_word, interval_record


def prime_numbers(limit):
    sieve = bytearray([1]) * (limit + 1)
    sieve[0:2] = b'\x00\x00'
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            for k in range(p*p, limit+1, p):
                sieve[k] = 0
    return [p for p in range(2, limit+1) if sieve[p]]


def error_interval(rational, low, high):
    assert rational < low or rational > high
    return (low-rational, high-rational) if rational < low else (rational-high, rational-low)


def phase(x):
    return x - x.__floor__()


def main():
    coefficients = [3, 4, 2]
    for k in range(3, 81):
        coefficients.append(k-1 if k % 2 else 3*((k-2)//2))
    ps, qs = [], []
    pprev, plast, qprev, qlast = 0, 1, 1, 0
    for k, a in enumerate(coefficients):
        p, q = a*plast+pprev, a*qlast+qprev
        assert gcd(p, q) == 1
        if k:
            assert q >= (4 if k % 2 else 9)
            assert q % (2 if k % 2 else 3) == 0
        ps.append(p); qs.append(q)
        pprev, plast, qprev, qlast = plast, p, qlast, q
    low, high = sorted([F(ps[-1], qs[-1]), F(ps[-2], qs[-2])])
    assert F(29, 9) < low < high < F(13, 4)
    for k in range(1, 79):
        lo, hi = error_interval(F(ps[k], qs[k]), low, high)
        assert qs[k]**2*hi < F(1, coefficients[k+1])
    prime_scores = []
    primes = prime_numbers(2000)
    for q in primes:
        p = (q*low+F(1, 2)).__floor__()
        assert p == (q*high+F(1, 2)).__floor__()
        lo, hi = error_interval(F(p, q), low, high)
        assert q*q*lo >= F(1, 2)
        prime_scores.append((q*q*lo, q*q*hi, q, p))
    min_prime = min(prime_scores)
    L, U, digits = certify_digits(64)
    pi_examples = []
    errors = {}
    for q in (2, 3, 5, 6):
        p = (q*L+F(1, 2)).__floor__()
        assert p == (q*U+F(1, 2)).__floor__()
        lo, hi = error_interval(F(p, q), L, U)
        errors[q] = (lo, hi)
        pi_examples.append({'p': p, 'q': q, 'absolute_error': interval_record(lo, hi)})
    assert errors[6][1] < min(errors[q][0] for q in (2, 3, 5))
    native_values = []
    for b in (3, 5, 7, 11):
        lo = evaluate_word(digits, F(b))
        hi = lo + F(3, (b-1)*b**64)
        if b == 3:
            assert lo-1 >= F(1, 27) and hi-1 <= F(1, 6)
        else:
            assert hi <= F(3, 4)
        native_values.append({'base': b, 'F': interval_record(lo, hi)})
    composition_checks = 0
    for x in [F(k, 37) for k in range(37)] + [L-3, U-3]:
        for a in (2, 3, 5, 7, 11):
            for b in (2, 3, 5, 7, 11):
                assert phase(a*phase(b*x)) == phase(a*b*x)
                assert phase(a*phase(b*x)) == phase(b*phase(a*x))
                composition_checks += 1
    result = {'status': 'exact finite checks passed; P2848 unresolved',
              'counterexample': {'coefficients_first_12': coefficients[:12],
                                'denominators_first_12': qs[:12],
                                'composite_denominators_checked': 80,
                                'convergent_error_bounds_checked': 78,
                                'prime_limit': 2000, 'primes_checked': len(primes),
                                'minimum_prime_sample': {'p': min_prime[3], 'q': min_prime[2],
                                    'normalized_error': interval_record(min_prime[0], min_prime[1])}},
              'pi_finite_counterexample': pi_examples,
              'native_digits_certified': 64, 'prime_base_values': native_values,
              'commuting_compositions_checked': composition_checks,
              'scope': 'General counterexample is not pi. Infinite claims follow from the note proofs, not finite tests.'}
    Path(__file__).with_name('prime_sampling_results.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
