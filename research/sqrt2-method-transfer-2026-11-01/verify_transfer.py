"""Exact finite checks for two attempted transfers of the uploaded pi method.

The accompanying rank obstruction is proved for all sufficiently large H.
These checks verify concrete matrices; they do not certify the original paper.
Only Python's standard library is required.
"""
from fractions import Fraction
from math import comb, factorial, isqrt
from pathlib import Path
import json


def prime_check(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def rank_mod(matrix, modulus):
    a = [[v % modulus for v in row] for row in matrix]
    if not a:
        return 0
    row_count, column_count = len(a), len(a[0])
    rank = 0
    for column in range(column_count):
        pivot = next((r for r in range(rank, row_count) if a[r][column]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inverse = pow(a[rank][column], -1, modulus)
        a[rank] = [(value * inverse) % modulus for value in a[rank]]
        pivot_row = a[rank]
        for r in range(rank + 1, row_count):
            multiplier = a[r][column]
            if multiplier:
                a[r][column:] = [
                    (x - multiplier * y) % modulus
                    for x, y in zip(a[r][column:], pivot_row[column:])
                ]
        rank += 1
        if rank == row_count:
            break
    return rank


def polynomial_restriction_matrix(H, w):
    """Columns Z^(2h)*(Z-1)^k, after removing nonzero sqrt(2)^k.

    The local square-root substitution is Y=1+t,
    X=u+sqrt(2)*(sqrt(1+t)-1). On u=0, Z=sqrt(1+t).
    """
    columns = [(h, k) for k in range(H // w + 1)
               for h in range(H - w*k + 1)]
    matrix = [[0] * len(columns) for _ in range(2*H + 1)]
    for c, (h, k) in enumerate(columns):
        for j in range(k + 1):
            matrix[2*h + j][c] = comb(k, j) * (-1)**(k-j)
    return columns, matrix


def logarithmic_jet_matrix(columns, jet_count, modulus):
    """Derivatives at z=0 of exp(h*z)*z^k, with z=ln(1+t).

    The formal change t=exp(z)-1 is invertible and preserves the rank of
    the first jet_count univariate jets. Entries are integers.
    """
    matrix = []
    for s in range(jet_count):
        row = []
        for h, k in columns:
            if s < k:
                row.append(0)
            else:
                falling = factorial(s) // factorial(s-k)
                row.append((falling * pow(h, s-k, modulus)) % modulus)
        matrix.append(row)
    return matrix


def exp_bounds(x, N=8):
    """Positive rational x: exact Taylor lower bound and geometric tail upper."""
    total = sum((x**k / factorial(k) for k in range(N+1)), Fraction())
    next_term = x**(N+1) / factorial(N+1)
    ratio = x / (N+2)
    assert 0 < ratio < 1
    return total, total + next_term/(1-ratio)


def main():
    modulus = 1_000_003
    assert prime_check(modulus)
    theta, w0, v0, K, m = Fraction(1, 10), Fraction(1), Fraction(1, 5), 1, 1
    assert K * theta**m < 1
    assert K * w0/v0 * theta**m == Fraction(1, 2) < 1
    cases = []
    for H, w in [(20, 1), (40, 2)]:
        columns, polynomial_matrix = polynomial_restriction_matrix(H, w)
        t_jets = 5*H
        full_target = sum(5*(H - 10*w*b)
                          for b in range((H-1)//(10*w) + 1))
        observed = rank_mod(polynomial_matrix, modulus)
        predicted = 2*H - w + 2
        # The accompanying proof gives the matching characteristic-zero upper
        # bound. Thus the modular lower bound certifies the exact rational rank.
        assert observed == predicted <= 2*H + 1 < t_jets
        assert len(columns) > full_target
        log_rank = rank_mod(logarithmic_jet_matrix(columns, t_jets, modulus), modulus)
        assert log_rank == t_jets
        cases.append({
            'H': H, 'w1': w, 'domain_monomials': len(columns),
            'full_target_jet_count': full_target, 'pure_t_target_count': t_jets,
            'sqrt_restriction_exact_rank': observed,
            'logarithmic_pure_t_rank': log_rank,
            'modular_prime': modulus,
            'scope': 'The log check concerns only the pure-t projection, not full multivariate surjectivity.',
        })

    # A large weight count illustrates the universal bound without a large matrix.
    H, w = 2000, 100
    monomials = sum(H-w*k+1 for k in range(H//w+1))
    target = sum(5*(H-10*w*b) for b in range((H-1)//(10*w)+1))
    assert monomials > target and 2*H+1 < 5*H

    # The repeated (a,d)=(0,0) coefficient rows for columns 1,Y are
    # (1,1) and (1,exp(-sqrt(2))). Their determinant is exp(-sqrt(2))-1.
    assert Fraction(7, 5)**2 < 2 < Fraction(3, 2)**2
    lower, _ = exp_bounds(Fraction(7, 5))
    _, upper = exp_bounds(Fraction(3, 2))
    assert 4 < lower and upper < 5
    # Monotonicity now gives 1/5 < exp(-sqrt(2)) < 1/4.
    determinant_bounds = [str(Fraction(1, 5)-1), str(Fraction(1, 4)-1)]

    # An admissible Section 4 parameter family also enters the forbidden
    # v0/w0 < 1/2 range for the square-root replacement.
    nu, theta4, A, B, C, eta = (Fraction(3), Fraction(9, 10),
                              Fraction(47, 50), Fraction(377, 400),
                              Fraction(21, 20), Fraction(1, 1000))
    assert 0 < theta4 < A < B < 1
    assert C > 1 and B < 1/C and B < C*theta4 < 1
    assert nu*(A*(1-eta)-theta4) - (1-theta4) > 0
    section4_cases = []
    for dimension in (30, 100):
        power = C**dimension
        centers = power.numerator // power.denominator
        ratio = 2*centers*theta4**dimension
        assert centers*theta4**dimension < 1
        assert centers*theta4**dimension/ratio == Fraction(1, 2)
        assert ratio < Fraction(1, 2)
        section4_cases.append({'m': dimension, 'K': centers,
                               'v0_over_w0': str(ratio),
                               'below_algebraic_curve_necessary_ratio': True})

    report = {
        'status': 'passed',
        'purpose': 'Audit direct transfers of the actual interpolation/translation/collision pipeline.',
        'paper_numeric_hypotheses': {'m': m, 'K': K, 'theta': str(theta),
                                     'w0': str(w0), 'v0': str(v0),
                                     'K_theta_power_m': '1/10',
                                     'K_w0_over_v0_theta_power_m': '1/2'},
        'finite_matrix_checks': cases,
        'large_weight_count_only': {'H': H, 'w1': w, 'domain_monomials': monomials,
                                   'full_target_jet_count': target,
                                   'pure_t_target_count': 5*H,
                                   'sqrt_rank_upper_bound': 2*H+1},
        'repeated_degree_coefficient_determinant': {
            'expression': 'exp(-sqrt(2))-1',
            'strict_rational_interval': determinant_bounds,
            'zero': False,
            'reason': 'After compensating for a nonperiodic translation, coefficient rows depend on the center j.',
        },
        'uniform_obstruction': 'For every integer w1>=1 and every integer H>=w1, pure-t image rank is 2H-w1+2, while the requested pure-t target has dimension 5H.',
        'section4_parameter_examples': {
            'nu': str(nu), 'theta': str(theta4), 'A': str(A),
            'B': str(B), 'C': str(C), 'eta': str(eta),
            'cases': section4_cases,
            'scope': 'These examples check Lemma 4.1 numerical conditions and interpolation ratios, not all later determinant inequalities.',
            'asymptotic_obstruction': 'v0/w0 = 2*floor(C^m)*theta^m <= 2*(C*theta)^m tends to zero. A square-root-curve pure-t projection requires v0/w0 >= 1/2 when wi >= w0/2.',
        },
        'not_claimed': [
            'An error in the original logarithmic interpolation theorem.',
            'Impossibility of every modified migration to sqrt(2).',
            'A completed proof of mu(sqrt(2))=2 using the original pipeline.',
        ],
    }
    path = Path(__file__).with_name('verification.json')
    path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
