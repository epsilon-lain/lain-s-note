#!/usr/bin/env python3
"""Exact finite verification of the stage-gap test for the fixed V_n grids."""

from math import comb, gcd, isqrt

from verify_round2_pi import pi_bounds


def grid_denominator(stage):
    return ((4 * stage + 1) * comb(4 * stage, 2 * stage)) >> stage.bit_count()


def check_scale(grids, X, B, K, gap_cap, pi_lower, pi_upper, scale):
    """Enumerate all viable pairs at X in every grid beneath the proven cap."""
    max_grid = gap_cap * B * B * X * X
    assert grids[-1] > max_grid
    outer_radius = B * K * gap_cap * scale
    pairs = []
    for n, V in enumerate(grids, start=2):
        if V > max_grid:
            break
        pi_times_V_lower = V * pi_lower
        pi_times_V_upper = V * pi_upper
        k_min = -((outer_radius - pi_times_V_lower) // scale)
        k_max = (pi_times_V_upper + outer_radius) // scale
        accepted = []
        for k in range(k_min, k_max + 1):
            common = gcd(k, V)
            if V // common > B * X:
                continue
            left = pi_times_V_lower - k * scale
            right = pi_times_V_upper - k * scale
            error_threshold = common * K * scale
            error_lower = max(left, -right, 0)
            if X * error_lower > error_threshold:
                continue
            if X * max(-left, right) > error_threshold:
                raise ArithmeticError("The rational pi enclosure is too coarse.")
            for previous_k in reversed(accepted):
                if k - previous_k > gap_cap:
                    break
                pairs.append((n, previous_k, k))
            accepted.append(k)
    return pairs


def main():
    scale = 10**1600
    pi_lower, pi_upper = (bound * scale for bound in pi_bounds(scale))
    assert pi_lower.denominator == pi_upper.denominator == 1
    pi_lower, pi_upper = int(pi_lower), int(pi_upper)
    grids = [grid_denominator(n) for n in range(2, 513)]
    assert all(next_V > V for V, next_V in zip(grids, grids[1:]))
    for B, K, gap_cap in [(1, 1, 1), (4, 4, 4), (8, 8, 4)]:
        for exponent in range(2, 10):
            stage = 2**exponent
            V = grids[stage - 2]
            V_previous = grids[stage - 3]
            assert V * (2 * stage - 1) == V_previous * (4 * stage - 1) * (4 * stage + 1)
            assert V > 8 * stage * V_previous
            X = isqrt(V // (2 * gap_cap * B * B))
            pairs = check_scale(grids, X, B, K, gap_cap, pi_lower, pi_upper, scale)
            expected = [(2, 987, 990)] if stage == 4 and B != 1 else []
            assert pairs == expected, (B, K, gap_cap, stage, pairs)
            print(f"B={B}, K={K}, M={gap_cap}, m={stage}, X has {X.bit_length()} bits, pairs={pairs}")


if __name__ == "__main__":
    main()
