#!/usr/bin/env python3
"""Verify balanced and shifted interval bounds at convergents of 1/pi. Requires mpmath.

The phase in W(t) is 2*n^2 radians. Index 0 is the convergent 0/1.
Run: python verify_resonant_windows.py
"""
from mpmath import iv

ENERGY_BOUNDS = {
    4: ("148.37294893724", "148.37294893725"),
    79: ("81.85721868", "81.85721872"),
    307: ("55.20558193474", "55.20558193475"),
    431: ("10500.70009373578", "10500.70009373579"),
}
RESIDUAL_BOUNDS = {
    4: ("2.37638757690902309045e-20", "2.37638757690902309047e-20"),
    431: ("1.07895873288330e-1417", "1.07895873288332e-1417"),
}


def convergents(alpha, last_index):
    """Fix each partial quotient by matching its two interval endpoints."""
    value = alpha
    prev_num, num, prev_den, den = 0, 1, 1, 0
    for index in range(last_index + 1):
        quotient = int(value.a)
        assert int(value.b) == quotient and value.a >= 0
        prev_num, num = num, quotient * num + prev_num
        prev_den, den = den, quotient * den + prev_den
        yield index, num, den
        value = 1 / (value - quotient)


def energy_interval(epsilon, den):
    """Bound the integral of the zero mode and all other modes uniformly."""
    if den % 2:
        zero_mode_fourth_power = 1
    else:
        # For den = 8 mod 16 the parity-twisted Gauss mean is the negative
        # of the untwisted Gauss mean, whose normalized modulus is sqrt(2).
        assert den % 16 == 8
        zero_mode_fourth_power = iv.mpf(1) / 4
    main = iv.pi * iv.atan2(1, 3) * zero_mode_fourth_power / (2 * epsilon)
    decay = iv.pi / (20 * epsilon)
    # |k_m/k_0| <= 3. Consecutive squares starting at 1 differ by >= 3.
    tail = 6 * iv.exp(-decay) / (1 - iv.exp(-3 * decay))
    assert tail.b < 1
    return main * (1 + iv.mpf([-1, 1]) * tail)**4


def odd_dual_sum(time_multiple, signed_error, num, den, cutoff):
    """Normalized Poisson sum at s*T, for odd den and balanced T=|error|."""
    assert den % 2
    assert signed_error.a > 0 or signed_error.b < 0
    sign = 1 if signed_error.a > 0 else -1
    scale = abs(signed_error)
    exponent = iv.pi**2 / (4 * den**2 * scale * (time_multiple - sign * iv.j))
    inverse = pow(4 * num, -1, den)
    total = iv.mpc(1, 0)
    for mode in range(1, cutoff + 1):
        # Half frequencies of odd index come from the weight (-1)^n/2.
        # The pair +/-mode thus has coefficient magnitude 1 (odd) or 2 (even).
        harmonic = (mode + (den if mode % 2 else 0)) // 2
        residue = (inverse * harmonic**2) % den
        phase = -2 * iv.pi * residue / den
        pair_weight = 1 if mode % 2 else 2
        total += pair_weight * iv.exp(iv.j * phase - mode**2 * exponent)
    decay = exponent.real
    # n >= cutoff+1, so all further differences of squares are >= 2*cutoff+3.
    tail = 2 * iv.exp(-(cutoff + 1)**2 * decay) / (1 - iv.exp(-(2 * cutoff + 3) * decay))
    component_error = iv.mpf([-1, 1]) * tail
    return total + iv.mpc(component_error, component_error)


def residual_interval(signed_error, num, den, time_multiple=1):
    sign = 1 if signed_error.a > 0 else -1
    at_one = odd_dual_sum(time_multiple, signed_error, num, den, cutoff=2)
    at_two = odd_dual_sum(2 * time_multiple, signed_error, num, den, cutoff=2)
    assert abs(at_one).a > 0
    normalized_prefactor = (time_multiple - sign * iv.j) / (2 * time_multiple - sign * iv.j)
    z = normalized_prefactor * (at_two / at_one)**2
    return abs(5 * z**2 - 6 * z + 2)


def check_interval(label, result, bounds):
    lower, upper = bounds
    assert iv.mpf(lower) < result < iv.mpf(upper), label
    print(f"{lower} < {label} < {upper}")


def main():
    iv.dps = 2600
    alpha = 1 / iv.pi
    for index, num, den in convergents(alpha, max(ENERGY_BOUNDS)):
        if index not in ENERGY_BOUNDS:
            continue
        epsilon = den * abs(den * alpha - num)
        assert epsilon.a > 0
        check_interval(f"E(T_{index})", energy_interval(epsilon, den), ENERGY_BOUNDS[index])
        if index in RESIDUAL_BOUNDS:
            signed_error = 2 - 2 * iv.pi * num / den
            result = residual_interval(signed_error, num, den)
            check_interval(f"|R(T_{index})|", result, RESIDUAL_BOUNDS[index])
            if index == 4:
                shifted_scale = 1 - iv.mpf("5.94e-20")
                shifted = residual_interval(signed_error, num, den, shifted_scale)
                check_interval("|R((1-5.94e-20) T_4)|", shifted, ("3.9971e-22", "3.9973e-22"))


if __name__ == "__main__":
    main()
