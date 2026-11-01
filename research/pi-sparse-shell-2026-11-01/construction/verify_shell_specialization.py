"""Exact rational parameter checks and finite modular rank certificates.

These computations check examples, not the conditional interpolation theorem.
The general proof is in SPARSE_SHELL_SPECIALIZATION.md.
"""
from fractions import Fraction as F
from math import comb, isqrt
from pathlib import Path
import json

PRIME = 1000003


def rational(x):
    return str(x)


def parameter_check(m):
    theta, B, C = F(1, 2), F(3, 5), F(13, 10)
    K = int(C ** m)
    w0 = B ** (-m)
    v0 = 2 * K * theta ** m * w0
    t = C * (1 + F(1, 8 * m))
    T = theta * t
    thetap = T * (1 + F(1, 8 * m))
    v0p = F(9, 10) * v0
    eta = F(K - 1) / (t ** m - 1)
    linear = 2 * F(K - 1, K) * F(m, m + 2) * t ** m / (t ** m - 1) * (T - theta)
    bracket = (1 - eta) ** 2 + eta ** 2 * t ** m + eta * (1 - eta) * (m + 2 - m / t)
    # c=ln(2)/4 < 7/40. Thus this is an exact rational upper bound.
    collision_upper = F(7, 20) * bracket / (K * v0 * (m + 2))
    cost_upper = linear + collision_upper
    deterministic_upper = 2 * F(K - 1, K) * (T * (m + 1) - theta * m) / (m + 2) + F(7, 40) / v0
    b = F(m, m + 2) * theta
    checks = {
        "original_first_volume": K * theta ** m < 1,
        "original_second_volume_half": K * w0 * theta ** m / v0 == F(1, 2),
        "critical_geometry": theta < B < C * theta < 1 and C * B < 1,
        "strict_shell_capacity": t ** m > K,
        "strict_packet_containment_margins": theta < T < thetap < 1 and 0 < v0p < v0,
        "single_center_volume": w0 * thetap ** m / v0p < 1,
        "power_bound": t ** m / K < F(4, 3),
        "linear_below_third": linear < F(1, 3),
        "deterministic_below_29_over_60": deterministic_upper < F(29, 60),
        "deterministic_below_comparison": deterministic_upper < 1 - b,
        "cost_below_two_fifths": cost_upper < F(2, 5),
        "cost_below_comparison": cost_upper < 1 - b,
    }
    assert all(checks.values()), (m, checks)
    return {
        "m": m,
        "K": K,
        "all_exact_checks_pass": True,
        "checks": checks,
        "linear": rational(linear),
        "collision_rational_upper": str(collision_upper),
        "collision_rational_upper_float": float(collision_upper),
        "cost_rational_upper": str(cost_upper),
        "cost_rational_upper_float": float(cost_upper),
        "deterministic_rational_upper": str(deterministic_upper),
        "deterministic_rational_upper_float": float(deterministic_upper),
        "comparison_1_minus_b": rational(1 - b),
        "single_center_volume_float": float(w0 * thetap ** m / v0p),
    }


def multiply(a, b, degree):
    out = [0] * (degree + 1)
    for i, x in enumerate(a):
        if not x:
            continue
        for j, y in enumerate(b[:degree + 1 - i]):
            out[i + j] = (out[i + j] + x * y) % PRIME
    return out


def rank_and_pivots(matrix):
    a = [r[:] for r in matrix]
    n, rank, pivots = len(a), 0, []
    for col in range(len(a[0])):
        p = next((i for i in range(rank, n) if a[i][col]), None)
        if p is None:
            continue
        a[rank], a[p] = a[p], a[rank]
        scale = pow(a[rank][col], -1, PRIME)
        a[rank] = [x * scale % PRIME for x in a[rank]]
        for i in range(rank + 1, n):
            k = a[i][col]
            if k:
                a[i] = [(x - k * y) % PRIME for x, y in zip(a[i], a[rank])]
        pivots.append(col)
        rank += 1
        if rank == n:
            break
    return rank, pivots


def determinant_mod(matrix):
    a = [r[:] for r in matrix]
    ans = 1
    for col in range(len(a)):
        p = next((i for i in range(col, len(a)) if a[i][col]), None)
        if p is None:
            return 0
        if p != col:
            a[p], a[col] = a[col], a[p]
            ans = -ans
        pivot = a[col][col]
        ans = ans * pivot % PRIME
        inv = pow(pivot, -1, PRIME)
        for i in range(col + 1, len(a)):
            fac = a[i][col] * inv % PRIME
            if fac:
                for k in range(col + 1, len(a)):
                    a[i][k] = (a[i][k] - fac * a[col][k]) % PRIME
    return ans % PRIME


def one_dimensional_finite_check(H):
    """An explicit m=1 cross-center shell at v0<w0.

    The finite injection a=j*d_s+beta is special to this 1-D example.
    It occupies a complete shell, with j=0 fixed, and is not used as a
    substitute for the general annulus transport proof.
    """
    K, theta, w0, v0, w1 = 2, F(1, 5), F(1), F(1, 2), F(1)
    cols = [(h, alpha) for alpha in range(H + 1) for h in range(H - alpha + 1)]
    rows = []
    for s in range(2 * H):
        cap = theta * (H - v0 * s) / w1
        d = -((-cap.numerator) // cap.denominator)
        for j in range(K):
            for beta in range(d):
                rows.append((j, s, beta, j * d + beta))
    assert all(a >= beta and w1*a <= H for j, s, beta, a in rows)
    assert K*theta < 1 and K*w0*theta/v0 < 1
    assert all(a == beta for j, s, beta, a in rows if j == 0)
    assert len({(s, a) for j, s, beta, a in rows}) == len(rows)
    degree = max(s for j, s, beta, a in rows)
    assert degree < PRIME
    logarithm = [0] + [((-1) ** (k + 1) * pow(k, -1, PRIME)) % PRIME for k in range(1, degree + 1)]
    matrices = {}
    for W in (0, 1):
        powers = {}
        for j in range(K):
            base = logarithm[:]
            base[0] = j * W
            pp = [[1] + [0] * degree]
            for k in range(H):
                pp.append(multiply(pp[-1], base, degree))
            powers[j] = pp
        matrix = []
        for j, s, beta, a in rows:
            line = []
            for h, alpha in cols:
                if alpha < a:
                    line.append(0)
                else:
                    z = powers[j][alpha - a]
                    line.append(comb(alpha, a) * sum(comb(h, k) * z[s - k] for k in range(min(h, s) + 1)) % PRIME)
            matrix.append(line)
        matrices[W] = matrix
    rank, pivots = rank_and_pivots(matrices[0])
    assert rank == len(rows), (H, rank, len(rows), len(cols))
    determinant_at_zero = determinant_mod([[row[p] for p in pivots] for row in matrices[0]])
    assert determinant_at_zero != 0
    return {
        "H": H,
        "m": 1,
        "K": K,
        "w0": str(w0),
        "v0": str(v0),
        "theta": str(theta),
        "rows": len(rows),
        "columns": len(cols),
        "all_s_a_distinct_across_centers": True,
        "anchor_fixed": True,
        "rank_at_zero_mod_prime": rank,
        "prime": PRIME,
        "minor_at_zero_mod_prime": determinant_at_zero,
        "same_minor_at_one_mod_prime": determinant_mod([[row[p] for p in pivots] for row in matrices[1]]),
        "pivot_column_indices": pivots,
        "scope": "Finite exact modular certificate, independent of unknown interpolation weight thresholds.",
    }


def main():
    assert all(PRIME % d for d in range(2, isqrt(PRIME) + 1))
    # The elementary series bound ln(2)<25/36 implies c<7/40.
    assert F(25, 144) < F(7, 40)
    result = {
        "status": "passed",
        "parameter_checks": [parameter_check(m) for m in (8, 10, 16, 24, 40, 80)],
        "finite_specialization_checks": [one_dimensional_finite_check(H) for H in (15, 20)],
        "logical_scope": "Finite checks do not prove the general theorem; its proof is conditional on source Theorem 2.1 and uses known pi transcendence.",
    }
    path = Path(__file__).with_name("shell_specialization_verification.json")
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"all_checks_pass": True, "output": str(path), "parameter_cases": len(result["parameter_checks"]), "finite_ranks": [r["rank_at_zero_mod_prime"] for r in result["finite_specialization_checks"]]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
