"""Exact certificates for pi carry survival, using only the standard library.

Run: python3 verify_survival_intervals.py
No decimal expansion or continued-fraction list of pi is input.  The program
certifies a native prefix, searches Stern--Brocot intervals, then independently
checks the final certificates using the homogeneous polynomial formula.
All assertions needed for certification use check(), so -O does not disable them.
"""
from fractions import Fraction as F
from math import gcd, isqrt
from pathlib import Path
import hashlib
import json

CAP = 10**12
BITS = 1024
COUNT = 96


def check(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def atan_bounds(v: int, terms: int) -> tuple[F, F]:
    total = sum((F((-1)**j, (2*j+1)*v**(2*j+1))
                 for j in range(terms)), F(0))
    next_term = F((-1)**terms, (2*terms+1)*v**(2*terms+1))
    return min(total, total+next_term), max(total, total+next_term)


def native_digits(kind: str) -> tuple[list[int], dict]:
    unit = 1 << BITS
    if kind == 'pi':
        a, b = atan_bounds(5, 250)
        c, d = atan_bounds(239, 90)
        lower, upper = 16*a-4*d, 16*b-4*c
        low = lower.numerator*unit//lower.denominator
        high = -((-upper.numerator*unit)//upper.denominator)
    elif kind == 'sqrt10':
        low = isqrt(10*unit*unit)
        high = low+1
    else:
        raise ValueError(kind)
    rem_low = rem_high = unit
    digits: list[int] = []
    for index in range(1, COUNT+1):
        first = low*rem_low//unit
        last = -((-high*rem_high)//unit)
        check(first//unit == last//unit, f'ambiguous digit {index}')
        digit = first//unit
        digits.append(digit)
        rem_low, rem_high = first-digit*unit, last-digit*unit
        check(0 < rem_low <= rem_high < unit, f'remainder {index}')
    return digits, {'kind': kind, 'count': COUNT, 'dyadic_bits': BITS,
                    'base_interval_numerators': [str(low), str(high)],
                    'digits': ''.join(map(str, digits)),
                    'digit_bytes_sha256': hashlib.sha256(bytes(digits)).hexdigest()}


def carry(p: int, q: int, digits: list[int]) -> tuple[int, int]:
    numerator = scale = 1
    for digit in digits:
        scale *= q
        numerator = p*numerator-digit*scale
    return numerator, scale


def direct(p: int, q: int, digits: list[int]) -> tuple[int, int]:
    degree = len(digits)
    numerator = p**degree-sum(d*q**j*p**(degree-j)
                               for j, d in enumerate(digits, 1))
    return numerator, q**degree


def minimum_survivor(digits: list[int]) -> tuple[tuple[int, int], tuple[int, int], int]:
    left, right = (3, 1), (4, 1)
    check(carry(*left, digits)[0] < -2, 'left outside')
    check(carry(*right, digits)[0] > 2, 'right outside')
    calls = 0
    while calls < 1_000_000:
        p, q = left[0]+right[0], left[1]+right[1]
        numerator, scale = carry(p, q, digits)
        calls += 1
        if numerator <= -2*scale:
            left = p, q
        elif numerator >= 2*scale:
            right = p, q
        else:
            return left, right, calls
    raise RuntimeError('search budget exceeded; no certificate emitted')


def verify_row(row: list[int], digits: list[int], h: int) -> None:
    baseline, a, b, c, d = row
    depth = baseline+h
    prefix = digits[:depth]
    check(len(prefix) == depth, 'insufficient digits')
    check(b > 0 and d > 0 and b*c-a*d == 1, 'Farey adjacency')
    check(3*b <= a and c <= 4*d, 'bracket outside [3,4]')
    low, low_scale = direct(a, b, prefix)
    high, high_scale = direct(c, d, prefix)
    mid, mid_scale = direct(a+c, b+d, prefix)
    check(low <= -2*low_scale and high >= 2*high_scale, 'outside signs')
    check(abs(mid) < 2*mid_scale, 'mediant does not survive')
    # Every survivor has denominator >= b+d and base < c/d.
    check(c**baseline < (b+d)**2*d**baseline, 'height inequality')


def baseline_depth(p: int, q: int) -> int:
    left, right, count = 1, q*q, 0
    while left < right:
        left *= p
        right *= q
        count += 1
    return count


def first_escape(p: int, q: int, digits: list[int]) -> int:
    numerator = scale = 1
    for depth, digit in enumerate(digits, 1):
        scale *= q
        numerator = p*numerator-digit*scale
        if abs(numerator) >= 2*scale:
            return depth
    raise ArithmeticError('escape not reached in certified prefix')


def build_certificate(digits: list[int], h: int) -> dict:
    largest_baseline = 0
    while 3**largest_baseline < CAP*CAP:
        largest_baseline += 1
    rows, calls = [], 0
    for baseline in range(largest_baseline+1):
        left, right, used = minimum_survivor(digits[:baseline+h])
        row = [baseline, *left, *right]
        verify_row(row, digits, h)
        rows.append(row)
        calls += used
    return {'denominator_cap': CAP, 'H_upper_bound': h,
            'baseline_max': largest_baseline,
            'native_digits_used_for_global_bound': largest_baseline+h,
            'certificate_rows': len(rows), 'search_evaluations': calls,
            'row_fields': ['B', 'left_p', 'left_q', 'right_p', 'right_q'],
            'rows': rows}


def independent_small_search(digits: list[int]) -> dict:
    depths = range(4, 13)
    observed = {depth: None for depth in depths}
    count = 0
    for q in range(1, 301):
        for p in range(3*q, 4*q+1):
            if gcd(p, q) != 1:
                continue
            count += 1
            for depth in depths:
                numerator, scale = direct(p, q, digits[:depth])
                if abs(numerator) < 2*scale and observed[depth] is None:
                    observed[depth] = (p, q)
    for depth in depths:
        left, right, _ = minimum_survivor(digits[:depth])
        candidate = (left[0]+right[0], left[1]+right[1])
        check(observed[depth] == candidate, 'exhaustive minimum mismatch')
    return {'q_cap': 300, 'reduced_bases_checked': count,
            'depths': list(depths),
            'minimum_survivors': {str(k): list(v) for k, v in observed.items()}}


def main() -> None:
    pi_digits, pi_info = native_digits('pi')
    control_digits, control_info = native_digits('sqrt10')
    pi_bound = build_certificate(pi_digits, 6)
    control_bound = build_certificate(control_digits, 5)
    witness = (355, 113)
    b = baseline_depth(*witness)
    e = first_escape(*witness, pi_digits)
    check((b, e, e-b) == (9, 15, 6), 'sharpness witness')
    result = {'status': 'all_exact_checks_passed',
              'scope': 'Finite denominator theorem only; no unbounded-depth or P2848 conclusion.',
              'pi_native': pi_info, 'pi_global_bound': pi_bound,
              'sharpness_witness': {'p': 355, 'q': 113, 'B': b, 'E': e, 'H': e-b},
              'sqrt10_native': control_info, 'sqrt10_finite_control': control_bound,
              'independent_exhaustive_check': independent_small_search(pi_digits)}
    destination = Path(__file__).with_name('survival_interval_results.json')
    # All certificate integers are retained; no floating-point output is used.
    destination.write_text(json.dumps(result, ensure_ascii=False, separators=(',', ':'))+'\n', encoding='utf-8')
    print(json.dumps({'status': result['status'],
                      'pi': {k: v for k, v in pi_bound.items() if k != 'rows'},
                      'control_rows': control_bound['certificate_rows'],
                      'exhaustive': result['independent_exhaustive_check'],
                      'output': str(destination)}, indent=2))


if __name__ == '__main__':
    main()
