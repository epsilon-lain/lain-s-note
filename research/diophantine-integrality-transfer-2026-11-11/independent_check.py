"""Independent elementary checks; no imported constructor or paper theorem.

The scalar map is evaluated as homogeneous integer forms, rather than via
successive rational beta/contact evaluations. These finite checks supplement
the proofs in INDEPENDENT_REVIEW.md; they do not verify the integrality paper.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, isqrt, lcm, prod
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
FORMS = [(51, 2504), (51, 2503), (51, 103), (1632, 80471)]


def valuation(n, prime):
    assert n != 0
    n = abs(n)
    value = 0
    while n % prime == 0:
        n //= prime
        value += 1
    return value


def factors(n):
    n = abs(n)
    result = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            n //= d
            result[d] = result.get(d, 0) + 1
        d += 1
    if n > 1:
        result[n] = result.get(n, 0) + 1
    return result


def run():
    prime = 149
    assert all(prime % d for d in range(2, isqrt(prime) + 1))
    roots = [0, 1, 2401, -343 * pow(32, -1, prime)]
    bad_squares = [((r - 2504) * pow(51, -1, prime)) % prime for r in roots]
    assert bad_squares == [59, 97, 109, 108]
    assert [pow(r, 74, prime) for r in bad_squares] == [148] * 4
    for x in range(prime):
        beta = (2504 + 51 * x * x) % prime
        assert all(v % prime for v in
                   [beta, beta - 1, beta - 2401, 32 * beta + 343])

    determinant_records = []
    exceptions = set()
    for i, j in combinations(range(4), 2):
        a, b = FORMS[i]
        c, d = FORMS[j]
        determinant = a * d - b * c
        fs = factors(determinant)
        exceptions.update(fs)
        determinant_records.append({'i': i, 'j': j,
                                    'determinant': determinant,
                                    'factorization': {str(p): e for p, e in fs.items()}})
    assert exceptions == {2, 3, 5, 7, 17}

    gcd_cases = 0
    for q in range(1, 51):
        for p in range(-50, 51):
            if gcd(p, q) != 1:
                continue
            vals = [a * p * p + b * q * q for a, b in FORMS]
            for row in determinant_records:
                g = gcd(vals[row['i']], vals[row['j']])
                assert row['determinant'] % g == 0
                gcd_cases += 1
            assert gcd(vals[0], vals[1]) == gcd(51, q)

    q_values = sorted(set(range(1, 76)) |
                      {149, 298, 149**2, 2*149**2, 3*17*149, 8*3*17*149**2})
    scalar_cases = 0
    valuation_cases = {'integral_at_149': 0, 'denominator_at_149': 0}
    for q in q_values:
        for p in range(-100, 101):
            if gcd(p, q) != 1:
                continue
            contact_product = prod(a*p*p + b*q*q for a, b in FORMS)
            denominator = q**8 + 149**2 * contact_product
            expected_gcd = (2**min(5, 2*valuation(q, 2)) *
                            3**(4*(q % 3 == 0)) *
                            17**(4*(q % 17 == 0)) *
                            149**(2*(q % 149 == 0)))
            for N in [1, 2, 7]:
                numerator = 149**(2*N+1) * contact_product
                assert gcd(numerator, denominator) == expected_gcd
                expected_v = 2*N - 1 if q % 149 == 0 else 2*N + 1
                assert valuation(numerator, 149) - valuation(denominator, 149) == expected_v
                kind = 'denominator_at_149' if q % 149 == 0 else 'integral_at_149'
                valuation_cases[kind] += 1
                scalar_cases += 1

    recovery_cases = 0
    for a in [1, 2, 3, 5, 6, 7, 10, 13, 17, 30, 51]:
        assert all(e == 1 for e in factors(a).values())
        for q in range(1, 33):
            for p in range(-30, 31):
                if gcd(p, q) != 1:
                    continue
                beta = 2504 + a*F(p, q)**2
                assert (beta.denominator == 1) == (q == 1)
                recovery_cases += 1
    assert F(1, 2).denominator > 1 and (2504 + 4*F(1, 2)**2).denominator == 1

    quartic_checks = 0
    for N in range(1, 129):
        x, s, t = F(1), F(1, N), F(N)
        loss = (s*x)**2 + (s*t-1)**2
        L = lcm(x.denominator, s.denominator, t.denominator)
        assert loss == F(1, N*N)
        assert L**4 * loss == N*N
        quartic_checks += 1

        # General quadratic-system wrapper: f=(x^2-2,y-x), b=f(x,y).
        x = y = F(3, 2)
        b0, b1 = x*x-2, y-x
        residuals = [x*x-2-b0, y-x-b1, s*b0, s*b1, s*t-1]
        L = lcm(*(z.denominator for z in [x, y, b0, b1, s, t]))
        integers = [L*L*r for r in residuals]
        assert all(a.denominator == 1 for a in integers)
        loss = sum(r*r for r in residuals)
        assert loss == F(1, 16*N*N)
        assert L**4*loss == sum(a*a for a in integers)
        assert (L**4*loss).denominator == 1 and L**4*loss >= 1
        quartic_checks += 1

    # A positive transfer control: sqrt(2), with exact integer/rational
    # certificates only. Adjacent Pell forms replace a real-value evaluator.
    pell = [(1, 1)]  # (B,A), so B^2-2A^2 is +1 or -1.
    while pell[-1][0] + pell[-1][1] < 1000:
        B, A = pell[-1]
        pell.append((B + 2*A, B + A))
    linear_certificates = 0
    polynomial_certificates = 0
    for q in range(1, 101):
        k = next(j for j, (B, A) in enumerate(pell[:-1]) if A+B >= 2*q)
        B0, A0 = pell[k]
        B1, A1 = pell[k+1]
        assert abs(B0*B0-2*A0*A0) == abs(B1*B1-2*A1*A1) == 1
        assert abs(A0*B1-A1*B0) == 1
        assert F(q, A0+B0) <= F(1, 2) and F(q, A1+B1) <= F(1, 2)
        c = F(1, 2*q*max(A0, A1))
        for p in range(q, 2*q+1):
            if gcd(p, q) != 1:
                continue
            x = F(p, q)
            assert x*x != 2
            if x*x < 2:
                assert (x+c)**2 <= 2
            else:
                assert x-c > 0 and (x-c)**2 >= 2
            linear_certificates += 1
            # Exact root case of P=X^2-2: no half-error allowance needed.
            cpoly = F(1, 4*q*q)
            assert p*p-2*q*q != 0
            if x*x < 2:
                assert (x+cpoly)**2 <= 2
            else:
                assert x-cpoly > 0 and (x-cpoly)**2 >= 2
            polynomial_certificates += 1

    out = {'checks': 'PASS',
           'method': 'Homogeneous integer forms, modular residue exhaustion, direct integer gcd, rational SOS normalization.',
           'nonsquare_residues_mod_149': bad_squares,
           'all_149_residues_checked': 149,
           'pair_determinants': determinant_records,
           'pair_gcd_instances_checked': gcd_cases,
           'scalar_map_instances_checked': scalar_cases,
           'valuation_regimes': valuation_cases,
           'squarefree_recovery_instances_checked': recovery_cases,
           'quartic_scaling_and_height_instances_checked': quartic_checks,
           'positive_sqrt2_linear_form_certificates_checked': linear_certificates,
           'positive_sqrt2_polynomial_gap_certificates_checked': polynomial_certificates,
           'scope': 'Elementary lemmas only. Finite checks supplement symbolic proofs; no verification of the claimed definition of Z or any pi approximation theorem.'}
    (ROOT/'independent_results.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    run()
