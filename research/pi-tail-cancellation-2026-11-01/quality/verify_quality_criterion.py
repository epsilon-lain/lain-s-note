#!/usr/bin/env python3
"""Exact checks for the formal anisotropic criterion, not a pi proof."""
from fractions import Fraction as F
from math import comb, prod
from pathlib import Path
import json


def capacity(H, coefficients):
    a, b = coefficients
    return sum((H - a * x) // b + 1 for x in range(H // a + 1))


def main():
    # ln(2)/4 = (1/2) sum_{k>=0} 1/((2k+1) 3^(2k+1)).
    # The first two terms give a strict lower bound. For k>=1 replace
    # 2k+1 by 3; strictness follows from k>=2.
    c_lower = F(1, 2) * (F(1, 3) + F(1, 81))
    c_upper = F(1, 6) + F(1, 6) * F(1, 27) / (1 - F(1, 9))
    assert c_lower == F(14, 81) and c_upper == F(25, 144)

    m, theta, eta = 500, F(1, 2), F(1, 100)
    B, A, E = F(17, 25), F(13, 5), F(1, 100)
    C = F(139, 100)
    K = (C ** m).numerator // (C ** m).denominator
    w0 = B ** (-m)
    v0 = 2 * K * theta ** m * w0
    assert K * theta ** m < 1
    assert K * w0 * theta ** m / v0 == F(1, 2)
    assert B < C * theta < 1 and B * C < 1
    u_sum, u_max = F(m), F(1)
    R = u_max - theta * u_sum / (m + 2)
    nu_uniform = [F(4)] * m
    nu_heterogeneous = [F(3)] * (m // 2) + [F(5)] * (m // 2)
    assert sum(nu_uniform) == sum(nu_heterogeneous) == 4 * m
    mean = F(sum(nu_uniform), m)
    assert mean == 4
    P_uniform, P_heterogeneous = prod(nu_uniform), prod(nu_heterogeneous)
    assert P_uniform == mean ** m
    assert P_heterogeneous == P_uniform * F(15, 16) ** (m // 2)
    assert P_heterogeneous < mean ** m
    D = theta * sum(nu_uniform) / (m + 2)
    L_uniform = eta**2 * K * theta**m * P_uniform / ((m + 1) * v0 * A**m)
    L_heterogeneous = eta**2 * K * theta**m * P_heterogeneous / ((m + 1) * v0 * A**m)
    assert L_uniform == eta**2 * F(68, 65)**500 / (2 * (m + 1))
    high_margin = A * (1 - eta) - R - D - E
    assert high_margin == F(4391, 62750) > 0
    assert c_lower * L_uniform > R + E
    assert c_upper * L_heterogeneous < R + E

    # Same quality mean 3 and same original weights (1,3).
    # The heterogeneous ratios (2,4) create more low types asymptotically
    # than the formal homogeneous ratios (3,3).
    lattice = []
    for H in (72, 144, 288, 576):
        values = []
        for coefficients in ((2, 12), (3, 9)):
            n = capacity(H, coefficients)
            a, b = coefficients
            area = F(H * H, 2 * a * b)
            # Unit squares cover the simplex and fit in its expanded simplex.
            upper = F((H + a + b)**2, 2 * a * b)
            assert area <= n <= upper
            values.append(n)
        assert values[0] > values[1]
        lattice.append({'H': H, 'heterogeneous_types': values[0],
                        'homogeneous_types': values[1]})

    # Equation (4) and both sides of (14), checked without logarithms.
    # These positive rational errors are symbolic surrogates, not pi errors.
    errors, weights, H, K_small = (F(1, 100), F(1, 1000)), (1, 3), 12, 3
    scalar_cases = 0
    anchor_cases = 0
    for a0 in range(H + 1):
        for a1 in range(H // 3 + 1):
            a = (a0, a1)
            if sum(w*x for w, x in zip(weights, a)) > H:
                continue
            for b0 in range(a0 + 1):
                for b1 in range(a1 + 1):
                    beta = (b0, b1)
                    d = tuple(x-y for x, y in zip(a, beta))
                    base = prod(e**z for e, z in zip(errors, d))
                    for j in range(1, K_small):
                        scalar = prod(F(comb(x, y)) * (2*j*e)**(x-y)
                                      for x, y, e in zip(a, beta, errors))
                        prefactor = prod(comb(x, y) for x, y in zip(a, beta)) * (2*j)**sum(d)
                        assert scalar == base * prefactor
                        assert 1 <= prefactor <= (4*(K_small-1))**sum(a)
                        assert base <= scalar <= base*(4*(K_small-1))**sum(a)
                        scalar_cases += 1
                    if beta == a:
                        assert prod(F(comb(x, y)) for x, y in zip(a, beta)) == 1
                        anchor_cases += 1

    # Exact normalization of the shell mean-quality coefficient, m fixed.
    shell = []
    for dimension in (8, 16, 40, 80):
        C_s, theta_s, B_s = F(13, 10), F(1, 2), F(3, 5)
        power = C_s ** dimension
        K_s = power.numerator // power.denominator
        t = C_s * (1 + F(1, 8*dimension))
        T = theta_s * t
        J = F(K_s-1, K_s) * F(dimension, dimension+2) * t**dimension / (t**dimension-1) * (T-theta_s)
        assert J > 0
        ratios = [F(2) + F(i+1, dimension*100) for i in range(dimension)]
        ratio_mean = sum(ratios) / dimension
        coordinate_cost = sum(nu*J/dimension for nu in ratios)
        assert coordinate_cost == ratio_mean * J
        shell.append({'m': dimension, 'J': str(J), 'J_display': float(J),
                      'quality_mean': str(ratio_mean),
                      'coordinate_sum_equals_mean_cost': True})

    out = {
        'status': 'passed',
        'scope': 'Synthetic parameter and structural checks only; no actual pi approximation tuple or interpolation threshold is certified.',
        'c_interval': [str(c_lower), str(c_upper)],
        'synthetic_comparison': {
            'm': m, 'mean_quality': str(mean),
            'arithmetic_threshold': str(R), 'source_row_quality_mean': str(D),
            'high_type_margin': str(high_margin),
            'uniform_collision_quantity': str(L_uniform),
            'uniform_collision_quantity_display': float(L_uniform),
            'heterogeneous_collision_quantity': str(L_heterogeneous),
            'heterogeneous_collision_quantity_display': float(L_heterogeneous),
            'uniform_passes_both_inequalities': True,
            'heterogeneous_fails_collision_inequality': True,
            'equal_mean_product_penalty': '(15/16)^250',
            'background_geometric_conditions_checked': True,
        },
        'finite_lattice_capacities': lattice,
        'scalar_prefactor_cases': scalar_cases,
        'zero_center_anchor_cases': anchor_cases,
        'shell_mean_quality_normalizations': shell,
        'not_tested': ['Actual approximation tuples for pi', 'Theorem 2.1',
                       'A lower bound on test determinant magnitudes',
                       'Cancellation between determinant terms',
                       'Bounded partial quotients of pi'],
    }
    path = Path(__file__).with_name('quality_verification.json')
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps({'status': out['status'], 'synthetic_uniform_passes': True,
                      'synthetic_heterogeneous_collision_fails': True,
                      'lattice_counts_checked': len(lattice)*2,
                      'scalar_prefactor_cases': scalar_cases,
                      'zero_center_anchor_cases': anchor_cases,
                      'shell_normalizations_checked': len(shell)}, indent=2))


if __name__ == '__main__':
    main()
