"""Finite exact checks for BOUNDARY_KERNEL_AND_MODULAR_STRUCTURE.md.

Uses the existing native certificate helper; no third-party dependencies.
Infinite summability/non-generation statements are proved in the note.
"""

import importlib.util
import json
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from pathlib import Path

helper_path = Path(__file__).resolve().parent.parent / 'pi-native-space-2026-10-31' / 'verify_native.py'
spec = importlib.util.spec_from_file_location('native_certificates', helper_path)
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)


def outward_decimal(value, rounding):
    with localcontext() as ctx:
        ctx.prec = 24
        ctx.rounding = rounding
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def main():
    original_low, original_high = helper.pi_bounds(160)
    grid = 2**600
    low = F((original_low * grid).__floor__(), grid)
    high = F((original_high * grid).__ceil__(), grid)
    assert low <= original_low < original_high <= high
    degree = 192
    digits = helper.certified_digits(F(1), F(1), low, high, degree)
    counts = [1]
    for n in range(1, degree + 1):
        counts.append(sum(digits[k-1] * counts[n-k] for k in range(1, n+1)))

    quotient = [counts[n] // (n+1) for n in range(degree+1)]
    product = [quotient[0]]
    for n in range(1, degree+1):
        product.append(quotient[n] - sum(digits[k-1]*quotient[n-k] for k in range(1, n+1)))
    for n in range(degree+1):
        recovered = sum(product[k]*counts[n-k] for k in range(n+1))
        assert recovered == quotient[n]
    assert product[0] == 1

    scaled_coefficients = []
    for n in [16, 32, 64, 128, 192]:
        endpoints = [F(n*n*product[n]) / low**n, F(n*n*product[n]) / high**n]
        lower, upper = min(endpoints), max(endpoints)
        scaled_coefficients.append({
            'n': n,
            'quantity': 'n^2 * f_n * pi^(-n)',
            'lower': outward_decimal(lower, ROUND_FLOOR),
            'upper': outward_decimal(upper, ROUND_CEILING),
            'bounds': 'exact rational interval, displayed with outward Decimal rounding'
        })

    block_cases = 0
    example = None
    for modulus in range(2, 14):
        for start in range(1, 41):
            if counts[start] % modulus:
                continue
            stop = start
            while stop+1 <= degree and counts[stop+1] % modulus == 0:
                stop += 1
            for n in range(start, stop+1):
                prediction = -sum(counts[j]*digits[n-j-1] for j in range(1, start))
                assert (digits[n-1] - prediction) % modulus == 0
                block_cases += 1
            if modulus == 4 and start == 3:
                example = {'modulus': 4, 'start': start, 'stop': stop,
                           'count_block': counts[start:stop+1],
                           'forced_digit_block': digits[start-1:stop]}

    return {
        'native_digit_certificate': {'count': degree, 'source': 'Machin rational bounds, 160 terms, then dyadic outward enclosure'},
        'harmonic_quotient_first_16': quotient[:16],
        'zero_relation_first_16': product[:16],
        'formal_D_quotient_and_C_recovery_orders': degree,
        'scaled_product_coefficients': scaled_coefficients,
        'divisible_block_identity_cases': block_cases,
        'mod4_example': example,
        'scope': 'Finite exact certificates only. Norm convergence, non-finite-generation, and modular nonrecurrence require the accompanying general proofs.'
    }


if __name__ == '__main__':
    print(json.dumps(main(), ensure_ascii=False, indent=2))
