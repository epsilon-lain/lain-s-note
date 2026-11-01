"""Exact certificate for a complete bulk-transport packet, using Taylor sums.

Standard-library only. This is one finite certificate, not an asymptotic
interpolation theorem or a lower bound on a period determinant's magnitude.
The independent verifier uses polynomial convolution and reduced elimination.
"""
from collections import defaultdict
from fractions import Fraction
from math import comb, factorial, prod
from pathlib import Path
import json

P = 65537
H = 20
ROWS = [(j, s, (b1, b2))
        for j in range(2) for b2 in range(3) for b1 in range(6)
        for s in range(H) if s + 4*b1 + 10*b2 < H]
COLS = [(h, (a1, a2))
        for a2 in range(5) for a1 in range(11) for h in range(H+1)
        if h + 2*a1 + 5*a2 <= H]


def entry(row, column, centers):
    """s! times the ordinary z^s u^beta Taylor coefficient, mod P."""
    j, s, beta = row
    h, alpha = column
    if any(b > a for a, b in zip(alpha, beta)):
        return 0
    value = 0
    for k1 in range(min(s, alpha[0]-beta[0])+1):
        for k2 in range(min(s-k1, alpha[1]-beta[1])+1):
            power = s-k1-k2
            value += (comb(alpha[0], beta[0])*comb(alpha[1], beta[1])
                      * comb(alpha[0]-beta[0], k1)
                      * comb(alpha[1]-beta[1], k2)
                      * (factorial(s)//factorial(power))*pow(h, power, P)
                      * pow(j*centers[0] % P, alpha[0]-beta[0]-k1, P)
                      * pow(j*centers[1] % P, alpha[1]-beta[1]-k2, P))
    return value % P


def pivots(matrix):
    a = [r[:] for r in matrix]
    out = []
    row = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(row, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        inverse = pow(a[row][c], -1, P)
        a[row] = [(x*inverse) % P for x in a[row]]
        for i in range(row+1, len(a)):
            factor = a[i][c]
            if factor:
                a[i] = [(x-factor*y) % P for x, y in zip(a[i], a[row])]
        out.append(c)
        row += 1
        if row == len(a):
            break
    return out


def det_mod(matrix):
    a = [r[:] for r in matrix]
    out = 1
    for c in range(len(a)):
        pivot = next((i for i in range(c, len(a)) if a[i][c]), None)
        if pivot is None:
            return 0
        if pivot != c:
            a[c], a[pivot] = a[pivot], a[c]
            out = -out
        diagonal = a[c][c]
        out = out*diagonal % P
        inverse = pow(diagonal, -1, P)
        for i in range(c+1, len(a)):
            factor = a[i][c]*inverse % P
            if factor:
                for k in range(c+1, len(a)):
                    a[i][k] = (a[i][k]-factor*a[c][k]) % P
    return out % P


def det_exact(matrix):
    a = [[Fraction(x) for x in r] for r in matrix]
    out = Fraction(1)
    for c in range(len(a)):
        pivot = next((i for i in range(c, len(a)) if a[i][c]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != c:
            a[c], a[pivot] = a[pivot], a[c]
            out = -out
        value = a[c][c]
        out *= value
        for i in range(c+1, len(a)):
            factor = a[i][c]/value
            for k in range(c+1, len(a)):
                a[i][k] -= factor*a[c][k]
    return out


def main():
    assert all(P % d for d in range(2, 257))
    assert 256*256 % P == P-1
    buckets = defaultdict(list)
    for j, s, beta in ROWS:
        if j == 0:
            buckets[s+sum(beta)].append(beta)
    packed = []
    for n, source in buckets.items():
        images = [(2*b[0]+j, b[1]) for j in range(2) for b in source]
        if all(sum(a) <= n and 2*a[0]+5*a[1] <= H for a in images):
            packed.append(n)
    packed.sort()
    assert packed == list(range(8, 20))
    target = []
    for j, s, beta in ROWS:
        n = s+sum(beta)
        a = (2*beta[0]+j, beta[1]) if n in packed else beta
        target.append((j, n-sum(a), a))
    assert len(set(target)) == len(ROWS) == 156
    assert len(COLS) == 242
    assert all(s >= 0 and 2*a[0]+5*a[1] <= H for j, s, a in target)
    assert all(r[0] == t[0] and r[1]+sum(r[2]) == t[1]+sum(t[2])
               for r, t in zip(ROWS, target))
    gamma = tuple(sum(t[2][i]-r[2][i] for r, t in zip(ROWS, target))
                  for i in range(2))
    coefficient = Fraction(1)
    blocks = []
    for n in packed:
        source = buckets[n]
        for j in range(2):
            images = [(2*b[0]+j, b[1]) for b in source]
            matrix = [[prod(comb(x, y) if x >= y else 0
                            for x, y in zip(a, b)) for a in images]
                      for b in source]
            determinant = det_exact(matrix)
            assert determinant == 2**sum(b[0] for b in source)
            coefficient *= determinant
            blocks.append({'N': n, 'center': j, 'size': len(source),
                           'integer_determinant': str(determinant)})
    assert coefficient == 2**60 and gamma == (93, 0)
    changed = sum(r != t for r, t in zip(ROWS, target))
    outside = len(set(target)-set(ROWS))
    packed_count = sum(s+sum(b) in packed for j, s, b in ROWS)
    assert (packed_count, changed, outside) == (66, 51, 33)
    period = [[entry(r, c, (1, 1)) for c in COLS] for r in target]
    selected = pivots(period)
    assert len(selected) == len(ROWS)
    rational_centers = (44*256*pow(7, -1, P) % P,
                        710*256*pow(113, -1, P) % P)
    source = [[entry(r, COLS[i], rational_centers) for i in selected]
              for r in ROWS]
    period_residue = det_mod([[line[i] for i in selected] for line in period])
    source_residue = det_mod(source)
    assert (period_residue, source_residue) == (9333, 44706)
    source_factorials = prod(factorial(s) for j, s, b in ROWS)
    target_factorials = prod(factorial(s) for j, s, b in target)
    ratio = Fraction(source_factorials, target_factorials)
    result = {
        'status': 'passed', 'prime': P, 'image_of_i': 256,
        'source_rows': len(ROWS), 'ambient_columns': len(COLS),
        'packed_buckets': packed, 'rows_in_packed_buckets': packed_count,
        'changed_rows': changed, 'outside_rows': outside,
        'source_rows_ordered': ROWS, 'target_rows_ordered': target,
        'selected_column_indices': selected,
        'selected_raw_columns': [{'h': COLS[c][0], 'alpha': COLS[c][1]}
                                 for c in selected],
        'period_target_factorial_scaled_determinant_mod_p': period_residue,
        'rational_source_factorial_scaled_determinant_mod_p': source_residue,
        'period_target_normalized_determinant_mod_p':
            period_residue*pow(target_factorials % P, -1, P) % P,
        'rational_source_normalized_determinant_mod_p':
            source_residue*pow(source_factorials % P, -1, P) % P,
        'normalized_transport_integer_coefficient': str(coefficient),
        'error_exponent': gamma,
        'transport_monomial_convention': 'integer_coefficient * (-zeta)^error_exponent',
        'factorial_scaled_transport_factor': str(ratio),
        'factorial_scaled_transport_integer_coefficient': str(coefficient*ratio),
        'exact_rational_block_checks': blocks,
        'scope': 'One complete finite packet. A nonzero formal-period polynomial '
                 'has nonzero value at 2*pi*i by known transcendence. Neither '
                 'magnitude, a full aggregate coefficient, nor an asymptotic '
                 'family of full-rank bulk packets is certified.'}
    Path(__file__).with_name('bulk_transport_verification.json').write_text(
        json.dumps(result, indent=2)+'\n')
    omitted = {'source_rows_ordered', 'target_rows_ordered', 'selected_column_indices',
               'selected_raw_columns', 'exact_rational_block_checks'}
    print(json.dumps({k: v for k, v in result.items() if k not in omitted}, indent=2))


if __name__ == '__main__':
    main()
