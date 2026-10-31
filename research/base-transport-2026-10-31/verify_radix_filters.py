"""Exact constructions from certified pi-base digits; no new pi approximation."""
from fractions import Fraction as F
from functools import reduce
from math import gcd, lcm
from pathlib import Path
import importlib.util
import json

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('transport', HERE/'verify_transport.py')
transport = importlib.util.module_from_spec(spec)
spec.loader.exec_module(transport)


def solve_filter(digits, N, L, R):
    """Solve M*b=0 with last coefficient 1 over Q; clear denominators."""
    matrix = [[F(digits[N+h+j]) for j in range(L)]
              + [-F(digits[N+h+L])] for h in range(R)]
    pivot_row = 0
    pivots = []
    for col in range(L):
        found = next((r for r in range(pivot_row, R) if matrix[r][col]), None)
        if found is None:
            continue
        matrix[pivot_row], matrix[found] = matrix[found], matrix[pivot_row]
        divisor = matrix[pivot_row][col]
        matrix[pivot_row] = [x/divisor for x in matrix[pivot_row]]
        for r in range(R):
            if r != pivot_row and matrix[r][col]:
                multiplier = matrix[r][col]
                matrix[r] = [x-multiplier*y for x, y in zip(matrix[r], matrix[pivot_row])]
        pivots.append((pivot_row, col))
        pivot_row += 1
    if any(all(x == 0 for x in row[:L]) and row[L] != 0 for row in matrix):
        return None
    answer = [F(0)]*L+[F(1)]
    for row, col in pivots:
        answer[col] = matrix[row][L]
    scale = lcm(*(x.denominator for x in answer))
    integers = [int(x*scale) for x in answer]
    common = reduce(gcd, map(abs, integers))
    integers = [x//common for x in integers]
    assert integers[-1] > 0
    assert all(sum(integers[j]*digits[N+h+j] for j in range(L+1)) == 0
               for h in range(R))
    return integers


def maximal_filter(digits, N, L):
    previous = solve_filter(digits, N, L, 0)
    for R in range(1, len(digits)-N-L+1):
        current = solve_filter(digits, N, L, R)
        if current is None:
            return R-1, previous
        previous = current
    raise AssertionError('Prefix too short to witness maximality')


def carry(p, q, digits):
    result = [1]
    for n, digit in enumerate(digits, 1):
        result.append(p*result[-1]-digit*q**n)
    return result


def check_filter(p, q, digits, N, b, R):
    L = len(b)-1
    assert gcd(p, q) == gcd(b[-1], q) == 1
    t = F(p, q)
    horizon = N+L+R
    A = carry(p, q, digits[:horizon])
    def J(n):
        return sum(b[j]*F(A[n+j], q**(n+j)) for j in range(L+1))
    Z = sum(b[j]*q**(L-j)*A[N+j] for j in range(L+1))
    assert Z == q**(N+L)*J(N)
    assert Z % q == (b[-1]*p**(N+L)) % q
    assert gcd(Z, q) == 1 and Z != 0
    for h in range(R):
        n = N+h
        forcing = sum(b[j]*digits[n+j] for j in range(L+1))
        assert forcing == 0
        assert J(n+1) == t*J(n)-forcing
    assert J(N+R) == t**R*J(N)
    B = sum(map(abs, b))
    survived = all(abs(value) < 2*q**n for n, value in enumerate(A))
    trigger = p**R >= 2*B*q**(N+L+R)
    if survived:
        assert not trigger
        assert t**R < 2*B*q**(N+L)
    return {'base': str(t), 'N': N, 'L': L, 'R': R, 'coefficients': b,
            'coefficient_height_l1': B, 'Z': str(Z),
            'horizon': horizon, 'survived': survived,
            'contradiction_trigger': trigger}


def main():
    _, _, digits = transport.certify_digits(128)
    assert digits[:5] == [3, 0, 1, 1, 0]
    profiles = []
    records = []
    bases = [(22, 7), (333, 106), (355, 113), (19, 6), (25, 8), (13, 4)]
    for N in range(3):
        for L in range(9):
            R, b = maximal_filter(digits, N, L)
            profiles.append({'N': N, 'L': L, 'maximal_R': R, 'coefficients': b})
            for p, q in bases:
                if gcd(b[-1], q) == 1:
                    records.append(check_filter(p, q, digits, N, b, R))

    b = [1, 3, -3]
    example = check_filter(22, 7, digits, 0, b, 2)
    assert example['Z'] == '4'
    assert [sum(b[j]*digits[n+j] for j in range(3)) for n in range(3)] == [0, 0, 4]
    matrix = [[digits[i+j] for j in range(3)] for i in range(3)]
    determinant = sum(matrix[0][j]*(matrix[1][(j+1)%3]*matrix[2][(j+2)%3]
                                   -matrix[1][(j+2)%3]*matrix[2][(j+1)%3])
                      for j in range(3))
    assert determinant == -4
    assert F(22, 7)**2*F(4, 49) == F(1936, 2401)

    A = carry(22, 7, digits[:4])
    repeat_Z = (22-7)*A[2]-7**2*7
    assert repeat_Z == -13
    assert F(A[4], 7**4)-F(A[2], 7**2) == F(-377, 2401)

    expected = {(0, 0): 0, (1, 0): 1, (2, 0): 0,
                (0, 1): 1, (1, 1): 0, (0, 2): 2}
    for profile in profiles:
        key = profile['N'], profile['L']
        if key in expected:
            assert profile['maximal_R'] == expected[key]
    # Independent integer check of the cost contradiction in Section 6.
    cost_cases = 0
    for h in range(9):
        threshold = max(3, (4**(h+1)+53)//54)
        for q in (threshold, threshold+1, threshold+10):
            for p in range(3*q, 4*q+1):
                if gcd(p, q) != 1:
                    continue
                baseline = 0
                while p**baseline < q**(baseline+2):
                    baseline += 1
                for D in range(3, baseline+h+1):
                    R = baseline+h-D
                    assert p**R < 2*q**(D+R)  # already fails with B=1
                    cost_cases += 1
                if cost_cases > 5000:
                    break
            if cost_cases > 5000:
                break
        if cost_cases > 5000:
            break

    result = {'status': 'exact finite checks passed; P2848 unresolved',
              'certified_native_digits': len(digits), 'prefix': digits[:12],
              'method': 'Machin rational intervals; Fraction elimination; integer comparisons',
              'filter_profiles': profiles, 'verified_integer_objects': len(records),
              'objects': records, 'quadratic_example': example,
              'three_step_matrix': matrix, 'determinant': determinant,
              'repeated_block_Z': repeat_Z,
              'fixed_extra_depth_cost_checks': cost_cases}
    (HERE/'radix_filter_results.json').write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps({k: result[k] for k in ['status', 'certified_native_digits',
        'verified_integer_objects', 'quadratic_example', 'determinant',
        'fixed_extra_depth_cost_checks']}, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
