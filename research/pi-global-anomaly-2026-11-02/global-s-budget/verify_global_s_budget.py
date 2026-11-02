"""Exact finite audit of the upper-half first-layer prime budget.

All universal conclusions are established in the companion proof. The tests
below use only Python standard-library integers and rational numbers.
"""
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path
import argparse
import json


def primes_upto(limit):
    flags = bytearray(b'\1') * (limit + 1)
    flags[:2] = b'\0\0'
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p*p:limit+1:p] = b'\0' * ((limit-p*p)//p+1)
    return [p for p in range(3, limit+1) if flags[p]]


def odd(x):
    assert x > 0
    return x // (x & -x)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--limit', type=int, default=1000)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    assert args.limit >= 21
    aa, bb = [1, 1], [0, 1]
    for n in range(1, args.limit):
        aa.append((2*n+1)*aa[n] + n*n*aa[n-1])
        bb.append((2*n+1)*bb[n] + n*n*bb[n-1])
    ps = primes_upto(args.limit)
    nf, lcm_prev = 1, 1
    classes = {'nonexceptional_deficit': 0,
               'exceptional_no_positive_excess': 0,
               'exceptional_positive_excess': 0}
    certificate_checks = 0
    dangerous = []
    for n in range(1, args.limit+1):
        nf *= n
        if n >= 2:
            lcm_prev = lcm_prev * (n-1) // gcd(lcm_prev, n-1)
        g = gcd(aa[n], 4*bb[n])
        r = Fraction(g, nf)
        E, D = odd(r.numerator), odd(r.denominator)
        H, P, P2 = 1, 1, 1
        exceptional = []
        lifted = []
        for p in ps:
            if p >= n:
                break
            if 2*p <= n:
                continue
            H *= p
            is_exception = aa[n-p] % p == 0
            if is_exception:
                P *= p
                exceptional.append(p)
            valuation = 0
            pp = p
            while g % pp == 0:
                valuation += 1
                pp *= p
            assert (valuation in (1, 2)) if is_exception else (valuation == 0)
            if valuation == 2:
                P2 *= p
                lifted.append(p)
            key = ('nonexceptional_deficit' if not is_exception else
                   'exceptional_positive_excess' if valuation == 2 else
                   'exceptional_no_positive_excess')
            classes[key] += 1
        if n == 21:
            assert P == 17
        if n >= 2:
            left = E * odd(n) * H*H
            right = D * P*P * odd(lcm_prev)
            sharp_right = D * P*P2 * odd(lcm_prev)
            assert right % left == 0
            assert sharp_right % left == 0
            assert P % P2 == 0
            binary = ((r.numerator & -r.numerator).bit_length()-1
                      -((r.denominator & -r.denominator).bit_length()-1))
            old_floor = Fraction(2, 2*n+1)/(r*r)
            budget_floor = (Fraction(2)**(1-2*binary)*odd(n)**2*H**4
                            / ((2*n+1)*P**4*odd(lcm_prev)**2))
            sharper_floor = (Fraction(2)**(1-2*binary)*odd(n)**2*H**4
                             / ((2*n+1)*P**2*P2**2*odd(lcm_prev)**2))
            assert budget_floor <= sharper_floor <= old_floor
            certificate_checks += 1
        # A rigorous finite S(n)<n/6 certificate: e>8/3, while the
        # exact comparison below gives P^6<(8/3)^n<e^n.
        assert 3**n * P**6 < 8**n
        # The submitted finite maximum is independently checked here
        # through this script's declared limit; do not extrapolate it.
        assert P**21 <= 17**n
        if P > 1:
            dangerous.append({'n': n, 'exceptional_primes': exceptional,
                              'lifted_primes': lifted,
                              'exceptional_product': str(P)})
    result = {
        'status': 'passed',
        'arithmetic': 'Python standard-library exact integers and fractions',
        'maximum_index': args.limit,
        'exact_first_layer_divisibility_checks': certificate_checks,
        'exact_two_layer_divisibility_checks': certificate_checks,
        'exact_finite_quality_budget_floors_checked': certificate_checks,
        'top_half_classification_counts': classes,
        'finite_S_less_than_n_over_6_certificates': args.limit,
        'finite_maximum_certificate': {
            'attained_at_n': 21,
            'exceptional_prime': 17,
            'exact_comparison_for_every_index': 'P_n^21 <= 17^n',
            'meaning': 'S(n)/n <= ln(17)/21 on the stated finite range'},
        'nonempty_exceptional_sets': dangerous,
        'scope': 'The finite range does not prove the global S(n) bound. '
                 'The exact aggregate divisibility and its PNT consequence '
                 'are proved in the companion note, not extrapolated.'}
    out = args.output or Path(__file__).with_name('global_s_budget_verification.json')
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items()
                      if k != 'nonempty_exceptional_sets'}, indent=2))


if __name__ == '__main__':
    main()
