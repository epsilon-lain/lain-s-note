"""Exact finite checks for MULTI_BAND_DEFICIT_CERTIFICATE.md.

No floating-point test is used for a divisibility, valuation, or quality floor.
The global asymptotic theorem is proved in the note, not inferred from this run.
"""
from fractions import Fraction
from math import gcd, lcm
from pathlib import Path
import json


def v2(x):
    assert x > 0
    return (x & -x).bit_length() - 1


def odd(x):
    return x >> v2(x)


def vp(x, p):
    assert x > 0
    value = 0
    while x % p == 0:
        x //= p
        value += 1
    return value


def power_two(k):
    return Fraction(2**k) if k >= 0 else Fraction(1, 2**(-k))


def cube_root_floor(n):
    b = 1
    while (b + 1)**3 <= n:
        b += 1
    return b


def primes_to(n):
    flags = bytearray(b'\1') * (n + 1)
    flags[0:2] = b'\0\0'
    for p in range(2, n + 1):
        if flags[p]:
            for k in range(p * p, n + 1, p):
                flags[k] = 0
    return [p for p in range(3, n + 1) if flags[p]]


def recurrence_to(limit):
    aa, bb, ff, uu = [1, 1], [0, 1], [1, 1], [1, 2]
    for n in range(1, limit):
        aa.append((2*n + 1)*aa[n] + n*n*aa[n-1])
        bb.append((2*n + 1)*bb[n] + n*n*bb[n-1])
        ff.append(ff[-1]*(n+1))
        numerator = 2*(2*n + 1)*uu[n] + 4*n*uu[n-1]
        assert numerator % (n+1) == 0
        uu.append(numerator // (n+1))
        assert Fraction(aa[-1], ff[-1]) == Fraction(uu[-1], 2**(n+1))
    return aa, bb, ff, uu


def refined_gain_divisor(c, n, ll):
    """Earlier proved two-neighbor cap and whole-prime deletion sieve."""
    cap_lcm = gcd(ll[n-1] // gcd(ll[n-1], n), ll[n+1] // (n+1))
    gg = gcd(c, odd(cap_lcm))
    hh = gcd(c // gg, gg)
    while True:
        next_h = gcd(gg, hh*hh)
        if next_h == hh:
            return gg // hh
        hh = next_h


def check(limit=1000):
    aa, bb, ff, uu = recurrence_to(limit)
    primes = primes_to(limit)
    ll = [1]
    for n in range(1, limit + 2):
        ll.append(lcm(ll[-1], n))
    largest_band = max(16, cube_root_floor(limit))
    prefix = [1]
    for b in range(1, largest_band + 1):
        assert 1 <= bb[b] <= aa[b]
        prefix.append(prefix[-1]*bb[b])

    certificate_checks = valuation_checks = positive_support_checks = 0
    improved_pairs = 0
    base_survivors, improved_survivors = [], []
    samples = []
    all_band_choices = set()
    for n in range(2, limit + 1):
        actual_g = gcd(aa[n], 4*bb[n])
        ratio = Fraction(actual_g, ff[n])
        e, d = odd(ratio.numerator), odd(ratio.denominator)
        binary = v2(ratio.numerator) - v2(ratio.denominator)
        assert ratio == power_two(binary)*Fraction(e, d)
        c = odd(uu[n])
        simple_gain = gcd(c, odd(ll[n-1]))
        sharp_gain = refined_gain_divisor(c, n, ll)
        assert simple_gain % e == 0
        assert sharp_gain % e == 0
        s = v2(uu[n])
        assert binary <= s - n
        old_floor = 2*power_two(2*(n-s))/((2*n+1)*sharp_gain**2)
        exact_analytic_floor = Fraction(2, 2*n+1)/ratio**2
        assert old_floor <= exact_analytic_floor
        if old_floor < 1:
            base_survivors.append(n)
        primary_band = cube_root_floor(n)
        band_choices = sorted({1, 2, 3, 8, 16, primary_band})
        all_band_choices.update(band_choices)
        for bands in band_choices:
            # These integer tests exactly encode p > max(sqrt(n), n/(B+1)).
            selected = [p for p in primes
                        if p <= n and p*p > n and p*(bands+1) > n]
            h = 1
            for p in selected:
                h *= p
            pp = gcd(c, h)
            f = gcd(h, prefix[bands])
            m = h // gcd(h, pp*f)
            assert d % m == 0
            assert (d*pp*f) % h == 0
            certificate_checks += 1
            e_high = 1
            for p in selected:
                exponent = vp(actual_g, p) - vp(ff[n], p)
                assert exponent <= 1
                if pp % p and f % p:
                    assert exponent == -1
                    valuation_checks += 1
                if exponent > 0:
                    assert pp % p == 0
                    e_high *= p**exponent
                    positive_support_checks += 1
            assert pp % e_high == 0
            # All other primes contribute only to the low-prime part.
            assert e % e_high == 0
            low = e // e_high
            for p in primes:
                if low % p == 0:
                    assert p not in selected
                    k = vp(low, p)
                    assert p**k <= n - 1
                    low //= p**k
            assert low == 1
            new_floor = old_floor*m*m
            assert old_floor <= new_floor <= exact_analytic_floor
            looser_exact_binary_floor = (power_two(1-2*binary)*h*h /
                                        ((2*n+1)*e*e*pp*pp*f*f))
            assert looser_exact_binary_floor <= exact_analytic_floor
            if m > 1:
                improved_pairs += 1
            if bands == primary_band:
                if new_floor < 1:
                    improved_survivors.append(n)
                if n in (21, 72, 107, 500, 1000):
                    samples.append({
                        'n': n, 'bands': bands,
                        'H': str(h), 'P': str(pp), 'F': str(f),
                        'M': str(m), 'E': str(e),
                        'gain_divisor_Ghat': str(sharp_gain),
                        'strict_floor_improvement_M_squared': str(m*m),
                        'old_floor_at_least_one': old_floor >= 1,
                        'new_floor_at_least_one': new_floor >= 1,
                        'all_generic_selected_prime_deficits_are_exact': True,
                    })
    return {
        'status': 'passed',
        'maximum_n': limit,
        'band_choices_tested': sorted(all_band_choices),
        'finite_divisibility_certificates_checked': certificate_checks,
        'generic_selected_prime_valuations_checked': valuation_checks,
        'positive_selected_prime_support_checks': positive_support_checks,
        'pairs_with_strict_M_squared_improvement': improved_pairs,
        'primary_band_rule': 'floor(cuberoot(n))',
        'gain_only_sieve_survivors_below_quality_one': base_survivors,
        'improved_sieve_survivors_below_quality_one': improved_survivors,
        'samples': samples,
        'scope': 'Finite exact checks only. The sufficient global bound on '
                 'log(P) is not proved, and no assertion about the bad '
                 'approximability of pi follows from this check.',
    }


def main():
    result = check()
    destination = Path(__file__).with_name('multi_band_verification.json')
    destination.write_text(json.dumps(result, indent=2) + '\n')
    # Keep the console report compact; complete samples are in the JSON file.
    brief = {k: v for k, v in result.items() if k != 'samples'}
    print(json.dumps(brief, indent=2))


if __name__ == '__main__':
    main()
