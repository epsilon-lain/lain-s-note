"""Independent checks for the finite irrational-numerator search certificate.

Uses Fraction arithmetic and independently summed arctangent partial sums for
the small exhaustive check. Does not import the certificate generator.
The finite checks support the accompanying proofs, not an infinite conclusion.
"""
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys


def atan_bounds(n, terms):
    partial = sum(
        (Fraction((-1)**k, (2*k+1)*n**(2*k+1)) for k in range(terms)),
        Fraction(),
    )
    adjacent = partial + Fraction((-1)**terms, (2*terms+1)*n**(2*terms+1))
    return min(partial, adjacent), max(partial, adjacent)


def binary_atan_bounds(n, terms):
    """Exact adjacent partial sums, without termwise fixed-point rounding.

    Each recursive block is normalized to start with power zero:
    F(a,b) = sum((-1/n**2)**(k-a)/(2*k+1), k=a,...,b-1).
    Splitting at m gives F(a,m) + (-1/n**2)**(m-a)*F(m,b).
    The final arctangent partial sum is F(0,terms)/n.
    """
    assert n >= 2 and terms > 0
    powers = {}

    def signed_power(length):
        if length not in powers:
            powers[length] = Fraction((-1)**length, n**(2*length))
        return powers[length]

    def block(start, stop):
        if stop-start == 1:
            return Fraction(1, 2*start+1)
        middle = (start+stop)//2
        return block(start, middle) + signed_power(middle-start)*block(middle, stop)

    partial = block(0, terms)/n
    adjacent = partial + Fraction((-1)**terms, (2*terms+1)*n**(2*terms+1))
    return min(partial, adjacent), max(partial, adjacent)


def abs_interval(lower, upper):
    assert lower * upper > 0
    return (lower, upper) if lower > 0 else (-upper, -lower)


def check_displayed_intervals(row, lower, upper, scale):
    p, q = int(row['p']), int(row['q'])
    beta_lo, beta_hi = Fraction(q*lower-p*scale, scale), Fraction(q*upper-p*scale, scale)
    assert beta_lo*beta_hi > 0
    assert row['beta_sign'] == (1 if beta_lo > 0 else -1)
    assert row['q_digits'] == len(str(q))
    for field in ['beta_interval', 'beta_scientific_interval']:
        shown_lo, shown_hi = map(Fraction, row[field])
        assert shown_lo <= beta_lo < beta_hi <= shown_hi
    shown_lo, shown_hi = map(Fraction, row['pi_minus_p_over_q_scientific_interval'])
    assert shown_lo <= beta_lo/q < beta_hi/q <= shown_hi
    shown_lo, shown_hi = map(Fraction, row['q_abs_beta_interval'])
    abs_lo, abs_hi = abs_interval(beta_lo, beta_hi)
    assert shown_lo <= q*abs_lo < q*abs_hi <= shown_hi


def valuation(n, prime):
    count = 0
    while n % prime == 0:
        n //= prime
        count += 1
    return count


def largest_prime_power(prime, limit):
    power, exponent = 1, 0
    while power * prime <= limit:
        power *= prime
        exponent += 1
    return power, exponent


def check_seed_bridge(p, q):
    assert q > 0 and math.gcd(p, q) == 1
    u, v = 304*q-105*p, 3*p-8*q
    rawp, rawq = 8*u+304*v, 3*u+105*v
    assert (rawp, rawq) == (72*p, 72*q)
    assert math.gcd(rawp, rawq) == 72
    return u, v


def verify(result_path, max_q=10000, tail_terms=1000):
    raw = result_path.read_bytes()
    result = json.loads(raw)
    lower = int(result['pi_scaled_lower'])
    upper = int(result['pi_scaled_upper'])
    digits = result['decimal_scale_digits']
    scale = 10**digits
    limit = 10**result['maximum_denominator_power_of_ten']
    assert lower < upper
    assert upper-lower == result['pi_width_units']

    for base in [2, 5, 239]:
        for count in [1, 2, 3, 7, 20, 70]:
            assert binary_atan_bounds(base, count) == atan_bounds(base, count)

    # Certify every supplied digit using independent exact rational sums.
    # Unlike the generator, this performs no termwise outward rounding.
    terms5 = result['arctan_counts']['atan5_terms']
    terms239 = result['arctan_counts']['atan239_terms']
    exact5_lo, exact5_hi = binary_atan_bounds(5, terms5)
    exact239_lo, exact239_hi = binary_atan_bounds(239, terms239)
    exact_pi_lo = 16*exact5_lo-4*exact239_hi
    exact_pi_hi = 16*exact5_hi-4*exact239_lo
    assert Fraction(lower, scale) <= exact_pi_lo < exact_pi_hi <= Fraction(upper, scale)

    # Rebuild the complete interval-stable prefix using rational intervals.
    lo, hi = Fraction(lower, scale), Fraction(upper, scale)
    coefficients = []
    while lo.numerator//lo.denominator == hi.numerator//hi.denominator:
        a = lo.numerator//lo.denominator
        coefficients.append(a)
        lo, hi = lo-a, hi-a
        if not lo or not hi:
            break
        lo, hi = 1/hi, 1/lo
    assert coefficients == result['certified_cf_coefficients']
    assert len(coefficients) == result['certified_cf_length']

    p2, p1, q2, q1 = 0, 1, 1, 0
    best = positive_best = None
    records, positive_records = [], []
    checked = 0
    expected_best_index = result['strict_quality_records'][-1]['convergent_index']
    local_pairs = {}
    for index, a in enumerate(coefficients):
        p, q = a*p1+p2, a*q1+q2
        assert math.gcd(p, q) == 1
        if index:
            assert abs(p*q1-p1*q) == 1
            assert q > q1  # For pi, a_1=7, so no equal initial denominators.
        p2, p1, q2, q1 = p1, p, q1, q
        if expected_best_index-2 <= index <= expected_best_index:
            local_pairs[index] = p, q
        if q > limit:
            stop = result['first_convergent_above_denominator_limit']
            assert index == stop['convergent_index']
            assert str(q) == stop['q']
            break
        signed_lo, signed_hi = q*lower-p*scale, q*upper-p*scale
        abs_lo, abs_hi = abs_interval(signed_lo, signed_hi)
        kl, ku = q*abs_lo, q*abs_hi
        if best is None or ku < best[0]:
            best = kl, ku, p, q, index
            records.append((p, q, index, kl, ku))
        else:
            assert kl > best[1]
        if signed_lo > 0:
            if positive_best is None or ku < positive_best[0]:
                positive_best = kl, ku, p, q, index
                positive_records.append((p, q, index, kl, ku))
            else:
                assert kl > positive_best[1]
        checked += 1
    else:
        raise AssertionError('No certified first convergent above the limit')
    assert checked == result['quality_checked_convergents']
    assert 2*best[1] < scale and 2*positive_best[1] < scale

    for rebuilt, saved in (
        (records, result['strict_quality_records']),
        (positive_records, result['positive_beta_quality_records']),
    ):
        assert len(rebuilt) == len(saved)
        for (p, q, i, kl, ku), row in zip(rebuilt, saved):
            assert (str(p), str(q), i) == (row['p'], row['q'], row['convergent_index'])
            shown_lo, shown_hi = map(Fraction, row['q_abs_beta_interval'])
            assert shown_lo <= Fraction(kl, scale) < Fraction(ku, scale) <= shown_hi
            check_seed_bridge(p, q)
            check_displayed_intervals(row, lower, upper, scale)

    # Exact local signed cancellation at the winning convergent.
    local_rows = result['best_local_cancellation_certificate']
    assert [row['convergent_index'] for row in local_rows] == [best[4]-2, best[4]-1, best[4]]
    for row in local_rows:
        index = row['convergent_index']
        assert (int(row['p']), int(row['q'])) == local_pairs[index]
        assert row['partial_quotient'] == coefficients[index]
        check_displayed_intervals(row, lower, upper, scale)
    local_coefficient = coefficients[best[4]]
    oldp, oldq = local_pairs[best[4]-2]
    prevp, prevq = local_pairs[best[4]-1]
    newp, newq = local_pairs[best[4]]
    assert (newp, newq) == (local_coefficient*prevp+oldp, local_coefficient*prevq+oldq)
    assert [row['beta_sign'] for row in local_rows] == [1, -1, 1]
    assert 21000*best[0] > scale and 20000*best[1] < scale
    seedu, seedv = check_seed_bridge(best[2], best[3])
    saved_seed = result['best_leibniz_seed_certificate']
    assert (int(saved_seed['u']), int(saved_seed['v'])) == (seedu, seedv)
    assert saved_seed['common_factor'] == 72
    for row in result['cancellation_examples']:
        check_displayed_intervals(row, lower, upper, scale)

    # Independent exact-rational Machin enclosure, then exhaustive denominators.
    a, b = atan_bounds(5, 70)
    c, d = atan_bounds(239, 25)
    pi_lo, pi_hi = 16*a-4*d, 16*b-4*c
    assert pi_lo < Fraction(lower, scale) < Fraction(upper, scale) < pi_hi
    brute_best, brute_records = None, []
    for q in range(1, max_q+1):
        lower_nearest = q*pi_lo + Fraction(1, 2)
        upper_nearest = q*pi_hi + Fraction(1, 2)
        p = lower_nearest.numerator//lower_nearest.denominator
        assert p == upper_nearest.numerator//upper_nearest.denominator
        sl, su = q*pi_lo-p, q*pi_hi-p
        kl, ku = abs_interval(q*sl, q*su)
        if brute_best is None or ku < brute_best[0]:
            brute_best = kl, ku, p, q
            brute_records.append((p, q))
        else:
            assert kl > brute_best[1]
    eligible = [(p, q) for p, q, *_ in records if q <= max_q]
    assert brute_records == eligible

    # Independent exact denominator checks for the original positive-tail family.
    partial = Fraction()
    for n in range(1, tail_terms+1):
        partial += Fraction(4, 4*n-3) - Fraction(4, 4*n-1)
        p, q = partial.numerator, partial.denominator
        D = 4*n-1
        power3, k = largest_prime_power(3, D)
        power5, ell = largest_prime_power(5, D)
        assert valuation(q, 3) == k and valuation(q, 5) == ell
        assert math.gcd(p, q) == 1
        assert q % (power3*power5) == 0
        assert 15*q > D*D
        assert D+2 <= 3*power3
        assert Fraction(2*q, D+2) >= Fraction(2, 3)
        assert Fraction(2*q*q, D+2) >= 2
        assert Fraction(2*q, D+2) > Fraction(2*D*D, 15*(D+2))
        assert Fraction(2*q*q, D+2) > Fraction(2*D**4, 225*(D+2))
        check_seed_bridge(p, q)

    assert check_seed_bridge(355, 113) == (-2923, 161)
    assert check_seed_bridge(333, 106) == (-2741, 151)
    assert 3*pi_lo-8 > 0 and 105*pi_lo-304 > 0
    assert 113*pi_hi-355 < 0 and 106*pi_lo-333 > 0
    assert (22+333, 7+106) == (355, 113)
    for q in range(1, 101):
        for p in range(1, 101):
            if math.gcd(p, q) == 1:
                check_seed_bridge(p, q)

    return {
        'status': 'passed',
        'input_file': result_path.name,
        'input_sha256': hashlib.sha256(raw).hexdigest(),
        'full_pi_interval_independently_certified': True,
        'full_pi_interval_method': 'Exact Fraction binary sums and adjacent alternating partial sums; generated interval contains the independently constructed Machin interval.',
        'full_pi_interval_exact_atan_terms': {'at_1_over_5': terms5, 'at_1_over_239': terms239},
        'binary_sum_vs_direct_fraction_checks': 18,
        'independently_rebuilt_cf_count': len(coefficients),
        'independently_checked_convergents': checked,
        'winner_convergent_index': best[4],
        'winner_denominator_digits': len(str(best[3])),
        'global_quality_record_count': len(records),
        'positive_quality_record_count': len(positive_records),
        'scientific_and_fixed_interval_displays_checked': True,
        'winner_local_cancellation': {'convergent_index': best[4], 'coefficient': local_coefficient, 'identity': 'beta_n = a_n*beta_(n-1) + beta_(n-2)', 'signs': [1, -1, 1]},
        'winner_quality_thresholds': '1/21000 < K < 1/20000, certified by exact integer comparisons.',
        'global_coverage_condition': 'All convergents up to Q checked and best K strictly below 1/2; Legendre excludes a better nonconvergent.',
        'small_exhaustive_max_q': max_q,
        'small_exhaustive_records': brute_records,
        'small_pi_enclosure_method': 'Independent Fraction sums of adjacent alternating arctangent partial sums: 70 terms at 1/5, 25 at 1/239.',
        'leibniz_family_exact_cases': tail_terms,
        'leibniz_checks': ['v3', 'v5', 'coprime reduced pair', 'quadratic denominator lower bound', 'derived beta and K lower bounds'],
        'seed_bridge_checks': 'All saved records, all tested Leibniz pairs, all primitive p,q in [1,100]^2, exact named examples and residual signs.',
        'scope': 'Finite certificate audit and finite checks of general identities; no proof of an infinite subsequence with K tending to zero.',
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', type=Path, default=Path(__file__).with_name('results_q_10pow5000.json'))
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('independent_verification.json'))
    parser.add_argument('--max-q', type=int, default=10000)
    parser.add_argument('--tail-terms', type=int, default=1000)
    args = parser.parse_args()
    if args.max_q < 1 or args.tail_terms < 1:
        parser.error('Check ranges must be positive')
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    report = verify(args.results, args.max_q, args.tail_terms)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
