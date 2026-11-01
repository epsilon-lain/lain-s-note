"""Independent finite review of bulk-compatible-transport/probe_bulk.py.

Uses a polynomial-convolution entry implementation and reduced-row pivot
path, not the probe's Taylor-allocation formula. Exact transport blocks are
checked over Q. This is one fixed H=20 certificate, with no external imports.
"""
from collections import defaultdict
from fractions import Fraction
from math import comb, factorial, prod
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
P = 65537
HEIGHT = 20
ROWS = [(j, s, (b1, b2)) for j in range(2) for b2 in range(3)
        for b1 in range(6) for s in range(HEIGHT)
        if s + 4*b1 + 10*b2 < HEIGHT]
COLS = [(h, (a1, a2)) for a2 in range(5) for a1 in range(11)
        for h in range(HEIGHT+1) if h + 2*a1 + 5*a2 <= HEIGHT]
FACT = [factorial(k) % P for k in range(HEIGHT+1)]
INVFACT = [pow(k, -1, P) for k in FACT]


def convolution(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] = (out[i+j] + x*y) % P
    return out


def entry(row, column, centers):
    j, s, beta = row
    h, alpha = column
    if any(b > a for b, a in zip(beta, alpha)):
        return 0
    poly = [1]
    for center, a, b in zip(centers, alpha, beta):
        n, base = a-b, j*center % P
        factor = [comb(a, b)*comb(n, k)*pow(base, n-k, P) % P
                  for k in range(n+1)]
        poly = convolution(poly, factor)
    value = sum(poly[k]*pow(h, s-k, P)*INVFACT[s-k]
                for k in range(min(s, len(poly)-1)+1))
    return value*FACT[s] % P


def lex_pivots(matrix):
    a = [row[:] for row in matrix]
    selected, k = [], 0
    for c in range(len(a[0])):
        pivot = next((r for r in range(k, len(a)) if a[r][c] % P), None)
        if pivot is None:
            continue
        a[k], a[pivot] = a[pivot], a[k]
        inverse = pow(a[k][c], -1, P)
        a[k] = [x*inverse % P for x in a[k]]
        # Fully reduced elimination differs from the probe's echelon path.
        for r in range(len(a)):
            if r == k:
                continue
            factor = a[r][c]
            if factor:
                a[r] = [(x-factor*y) % P for x, y in zip(a[r], a[k])]
        selected.append(c)
        k += 1
        if k == len(a):
            break
    return selected


def determinant_mod(matrix):
    a = [row[:] for row in matrix]
    answer = 1
    for k in range(len(a)):
        pivot = next((r for r in range(k, len(a)) if a[r][k] % P), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            answer = -answer
        diagonal = a[k][k] % P
        answer = answer*diagonal % P
        inverse = pow(diagonal, -1, P)
        for r in range(k+1, len(a)):
            quotient = a[r][k]*inverse % P
            if quotient:
                for c in range(k+1, len(a)):
                    a[r][c] = (a[r][c]-quotient*a[k][c]) % P
    return answer % P


def determinant_q(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    answer = Fraction(1)
    for k in range(len(a)):
        pivot = next((r for r in range(k, len(a)) if a[r][k]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            answer = -answer
        diagonal = a[k][k]
        answer *= diagonal
        for r in range(k + 1, len(a)):
            quotient = a[r][k] / diagonal
            for c in range(k + 1, len(a)):
                a[r][c] -= quotient * a[k][c]
    return answer


def main():
    assert all(P % d for d in range(2, 257))
    assert 256 * 256 % P == P - 1
    buckets = defaultdict(list)
    for j, s, beta in ROWS:
        if j == 0:
            buckets[s + sum(beta)].append(beta)
    # Independently require both exact support inequalities for every target.
    packed = []
    for n, source in buckets.items():
        candidates = [(2*b[0] + j, b[1]) for j in range(2) for b in source]
        if all(sum(a) <= n and 2*a[0] + 5*a[1] <= 20 for a in candidates):
            packed.append(n)
    packed.sort()
    target = []
    blocks = []
    coefficient = Fraction(1)
    gamma = [0, 0]
    for j, s, beta in ROWS:
        n = s + sum(beta)
        a = (2*beta[0] + j, beta[1]) if n in packed else beta
        target.append((j, n - sum(a), a))
    for n in packed:
        source = buckets[n]
        for j in range(2):
            images = [(2*b[0] + j, b[1]) for b in source]
            matrix = [[prod(comb(a_i, b_i) if a_i >= b_i else 0
                            for a_i, b_i in zip(a, b)) for a in images]
                      for b in source]
            integer = determinant_q(matrix)
            expected = 2 ** sum(b[0] for b in source)
            assert integer == expected
            coefficient *= integer
            for i in range(2):
                gamma[i] += sum(a[i] for a in images) - sum(b[i] for b in source)
            blocks.append({'N': n, 'center': j, 'size': len(source),
                           'integer_determinant': str(integer)})
    assert packed == list(range(8, 20))
    assert len(set(target)) == len(ROWS) == 156
    assert all(s >= 0 and 2*a[0] + 5*a[1] <= 20 for j, s, a in target)
    outside = sum(s + 4*a[0] + 10*a[1] >= 20 for j, s, a in target)
    changed = sum(a != b for a, b in zip(ROWS, target))
    assert outside == 33
    period = [[entry(row, column, (1, 1)) for column in COLS]
              for row in target]
    selected = lex_pivots(period)
    assert len(selected) == 156
    centers = (44*256*pow(7, -1, P) % P,
               710*256*pow(113, -1, P) % P)
    source = [[entry(row, COLS[c], centers) for c in selected]
              for row in ROWS]
    period_residue = determinant_mod([[row[c] for c in selected]
                                     for row in period])
    source_residue = determinant_mod(source)
    assert period_residue == 9333 and source_residue == 44706
    factorial_ratio = Fraction(prod(factorial(s) for _, s, _ in ROWS),
                               prod(factorial(s) for _, s, _ in target))
    result = {
        'status': 'passed', 'prime': P, 'image_of_i': 256,
        'source_rows': len(ROWS), 'ambient_columns': len(COLS),
        'packed_buckets': packed, 'changed_rows': changed,
        'packed_rows': sum(2*len(buckets[n]) for n in packed),
        'outside_rows': outside, 'period_target_rank_at_formal_period_1': len(selected),
        'selected_column_indices': selected,
        'selected_raw_columns': [{'h': COLS[c][0], 'alpha': list(COLS[c][1])}
                                 for c in selected],
        'period_target_factorial_scaled_determinant_mod_p': period_residue,
        'rational_source_factorial_scaled_determinant_mod_p': source_residue,
        'normalized_transport_integer_coefficient': str(coefficient),
        'normalized_transport_monomial': '2^60 * (-zeta_1)^93',
        'error_exponent': gamma,
        'scaled_to_normalized_transport_factor': str(factorial_ratio),
        'factorial_scaled_transport_integer_coefficient': str(coefficient*factorial_ratio),
        'exact_rational_block_checks': blocks,
        'scope': 'One fixed complete H=20 packet on common raw columns. '
                 'Transcendence promotes the nonzero formal-period polynomial '
                 'to the actual period, but supplies neither magnitude nor an '
                 'asymptotic family of compatible full-rank packed packets. '
                 'A nonzero individual target-set term does not prove that '
                 'the full eta^33 zeta_1^93 coefficient survives summation.'}
    (HERE / 'bulk_review_verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in
                      ('selected_column_indices', 'selected_raw_columns',
                       'exact_rational_block_checks')}, indent=2))


if __name__ == '__main__':
    main()
