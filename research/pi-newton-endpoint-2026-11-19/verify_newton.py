"""Verify modular rank bounds for the two-coordinate Newton basis in the pi jet matrix."""
from functools import lru_cache
from math import comb, factorial

K, W1, W2, V, INV_THETA = 2, 1, 3, 2, 4
P1, Q1, P2, Q2 = 22, 7, 355, 113
PRIME, SQRT_MINUS_ONE = 65537, 256


def row_packet(height):
    if height <= 0:
        return []
    rows = []
    for j in range(K):
        for b1 in range((height - 1) // (W1 * INV_THETA) + 1):
            for b2 in range((height - 1 - W1 * INV_THETA * b1) // (W2 * INV_THETA) + 1):
                time_budget = height - 1 - INV_THETA * (W1 * b1 + W2 * b2)
                for s in range(time_budget // V + 1):
                    rows.append((j, s, b1, b2))
    return rows


def columns(height):
    # The last-coordinate degree is primary, in descending order.
    return [(h, a1, a2)
            for a2 in reversed(range(height // W2 + 1))
            for a1 in range((height - W2 * a2) // W1 + 1)
            for h in range(height - W2 * a2 - W1 * a1 + 1)]


@lru_cache(None)
def polynomial(j, degree, coordinate, modulus, sqrt_minus_one, shifted):
    # Coefficients of f(T) before T is replaced by q*z+v.
    p = (P1, P2)[coordinate - 1]
    coefficients = [1]
    for idx in range(degree):
        node = idx % K if coordinate == 2 and shifted else 0
        constant = 2 * sqrt_minus_one * (j - node) * p
        coefficients = [((coefficients[d] * constant if d < len(coefficients) else 0)
                         + (coefficients[d - 1] if d else 0)) % modulus
                        for d in range(len(coefficients) + 1)]
    return coefficients


@lru_cache(None)
def jet(j, degree, b, coordinate, modulus, sqrt_minus_one, shifted):
    # [v^b] f(q*z+v) in ascending z degree.
    poly = polynomial(j, degree, coordinate, modulus, sqrt_minus_one, shifted)
    q = (Q1, Q2)[coordinate - 1]
    return tuple(c * comb(d, b) * pow(q, d - b, modulus) % modulus
                 for d, c in enumerate(poly) if d >= b)


def entry(row, col, modulus, sqrt_minus_one, shifted=True):
    j, s, b1, b2 = row
    h, a1, a2 = col
    j1 = jet(j, a1, b1, 1, modulus, sqrt_minus_one, shifted)
    j2 = jet(j, a2, b2, 2, modulus, sqrt_minus_one, shifted)
    # s! [z^s] e^{h z} j1(z) j2(z), with z=log(1+t).
    return sum(x * y * (factorial(s) // factorial(s - d1 - d2))
               * pow(h, s - d1 - d2, modulus)
               for d1, x in enumerate(j1[:s+1])
               for d2, y in enumerate(j2[:s+1-d1])) % modulus


def greedy_basis(rows, cols, shifted=True):
    echelon, selected = {}, []
    for col in cols:
        vec = [entry(row, col, PRIME, SQRT_MINUS_ONE, shifted) for row in rows]
        for pivot in sorted(echelon):
            coefficient = vec[pivot]
            if coefficient:
                vec = [(x - coefficient * y) % PRIME for x, y in zip(vec, echelon[pivot])]
        pivot = next((idx for idx, x in enumerate(vec) if x), None)
        if pivot is not None:
            inv = pow(vec[pivot], -1, PRIME)
            echelon[pivot] = [(x * inv) % PRIME for x in vec]
            selected.append(col)
            if len(selected) == len(rows):
                break
    return selected


def lifted_i(prime, root, precision):
    level = 1
    while level < precision:
        level = min(2 * level, precision)
        modulus = prime ** level
        root = (root - (root*root + 1) * pow(2 * root, -1, modulus)) % modulus
    assert (root*root + 1) % (prime**precision) == 0
    return root


def residue_valuation(x, prime, precision):
    if not x:
        return precision
    v = 0
    while x % prime == 0:
        x //= prime
        v += 1
    return v


def determinant_valuation(matrix, prime, precision):
    """Determinant valuation from row operations in Z/p^N; require value < N."""
    modulus = prime**precision
    value = 0
    for k in range(len(matrix)):
        pivot_row = min(range(k, len(matrix)),
                        key=lambda i: residue_valuation(matrix[i][k], prime, precision))
        matrix[k], matrix[pivot_row] = matrix[pivot_row], matrix[k]
        pivot = matrix[k][k]
        v = residue_valuation(pivot, prime, precision)
        value += v
        assert value < precision, "raise the precision"
        prime_power = prime**v
        inverse = pow(pivot // prime_power, -1, modulus)
        for i in range(k + 1, len(matrix)):
            x = matrix[i][k]
            if x:
                assert x % prime_power == 0
                multiplier = (x // prime_power) * inverse % modulus
                for j in range(k, len(matrix)):
                    matrix[i][j] = (matrix[i][j] - multiplier * matrix[k][j]) % modulus
                assert matrix[i][k] == 0
    return value


def experiment(height):
    rows = row_packet(height)
    chosen = greedy_basis(rows, columns(height))
    assert len(rows) == len(chosen)
    nominal_power = sum(col[2] for col in chosen) - sum(row[3] for row in rows)
    counts = [(k, sum(a2 >= K*k for _, _, a2 in chosen), len(row_packet(height-V*k)),
               sum(b2 >= k for _, _, _, b2 in rows)) for k in range(1, height//(K*W2)+1)]
    assert all(high <= bound for _, high, bound, _ in counts)
    gain = sum(max(high - already_used, 0) for _, high, _, already_used in counts)
    # The split prime 113 has conjugate roots of -1, ±15 modulo 113.
    # Evaluate both to obtain the exact rational 113-content of the Gaussian determinant.
    precision = nominal_power + 100
    modulus = Q2**precision
    values, raw_values = [], []
    for base_root in (15, -15):
        root = lifted_i(Q2, base_root, precision)
        for shifted, target in ((True, values), (False, raw_values)):
            matrix = [[entry(row, col, modulus, root, shifted) for col in chosen] for row in rows]
            target.append(determinant_valuation(matrix, Q2, precision))
    saved = min(values)
    # Conversion from time z to t is unit triangular, row factors s! and q1 are 113-units.
    assert gain <= saved
    actual_den_power = max(0, nominal_power - saved)
    print(f"H={height}; rows={len(rows)}; layers={counts}")
    print(f"nominal={nominal_power}; guaranteed={gain}; local_valuations={values}; "
          f"true_denominator_113_exponent={actual_den_power}")
    raw_den_power = max(0, nominal_power - min(raw_values))
    print(f"same exponents, raw monomials: local valuations={raw_values}; "
          f"exact denominator exponent={raw_den_power}")
    # Compare the ascending last-degree order. Every nonzero determinant has
    # C2>=0 since a nonzero entry requires a2>=b2. Here the minimum is attained.
    cheap_order = sorted(columns(height), key=lambda c: (c[2], c[1], c[0]))
    cheap_raw = greedy_basis(rows, cheap_order, shifted=False)
    cheap_newton = greedy_basis(rows, cheap_order)
    assert cheap_raw == cheap_newton and len(cheap_raw) == len(rows)
    cheap_counts = [sum(a2 == level for _, _, a2 in cheap_raw) for level in (0, 1)]
    assert sum(cheap_counts) == len(rows)
    assert sum(c[2] for c in cheap_raw) - sum(row[3] for row in rows) == 0
    # In every nonzero permutation all a2=b2. Then every last-coordinate
    # factor is monic and contributes its leading transverse coefficient 1.
    print(f"cheapest raw and Newton: degrees 0/1 counts={cheap_counts}; "
          "nominal=0; exact denominator exponent=0; independent of r2")
    return nominal_power, gain, values, actual_den_power, raw_values, raw_den_power


if __name__ == '__main__':
    assert (SQRT_MINUS_ONE**2 + 1) % PRIME == 0
    assert all(PRIME % d for d in range(2, 257))
    assert experiment(14) == (78, 28, [31, 31], 47, [25, 25], 53)
    assert experiment(18) == (172, 71, [74, 74], 98, [52, 52], 120)
