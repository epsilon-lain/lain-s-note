"""Exact finite certificates for the native pi-base model.

All decisions use Fraction arithmetic. No decimal expansion of pi is used.
The accompanying README supplies the infinite statements and their proofs.
Run: python3 verify_native.py > results.json
"""

import json
from fractions import Fraction as F
from math import gcd


def atan_bounds(denominator, terms):
    value = sum((F((-1)**k, (2*k+1)*denominator**(2*k+1))
                 for k in range(terms)), F(0))
    adjacent = value + F((-1)**terms, (2*terms+1)*denominator**(2*terms+1))
    return min(value, adjacent), max(value, adjacent)


def pi_bounds(terms=100):
    low5, high5 = atan_bounds(5, terms)
    low239, high239 = atan_bounds(239, terms)
    return 16*low5-4*high239, 16*high5-4*low239


def certified_digits(low, high, base_low, base_high, count):
    digits = []
    for _ in range(count):
        assert 0 <= low <= high <= 1
        low, high = base_low*low, base_high*high
        digit = low.numerator // low.denominator
        assert digit == high.numerator // high.denominator, 'Increase interval precision'
        low, high = low-digit, high-digit
        assert 0 <= low <= high < 1
        assert digit in (0, 1, 2, 3)
        digits.append(digit)
    return digits


def main():
    base_low, base_high = pi_bounds()
    assert 3 < base_low < base_high < 4
    boundary = certified_digits(F(1), F(1), base_low, base_high, 48)
    assert boundary[:6] == [3, 0, 1, 1, 0, 2]
    four_tail = certified_digits(4-base_high, 4-base_low, base_low, base_high, 32)

    returns = [1]
    for length in range(1, 41):
        returns.append(sum(boundary[k-1]*returns[length-k]
                           for k in range(1, length+1)))
    assert returns[:7] == [1, 3, 9, 28, 88, 276, 867]
    assert F(returns[4], returns[3]) == F(22, 7)

    # Independently count labelled paths by their current graph state.
    states = {0: 1}
    for length in range(1, 31):
        next_states = {}
        for state, count in states.items():
            digit = boundary[state]
            next_states[0] = next_states.get(0, 0) + digit*count
            next_states[state+1] = next_states.get(state+1, 0) + count
        states = next_states
        assert states[0] == returns[length]
        assert sum(states.values()) == sum(returns[:length+1])

    # Coefficients of (1-D(z))*A(z) through the certified length.
    for length in range(1, 41):
        coefficient = returns[length] - sum(
            boundary[k-1]*returns[length-k] for k in range(1, length+1)
        )
        assert coefficient == 0

    radius = F(1, 2)
    tail_bound = sum(F(boundary[k-1])*radius**k for k in range(2, 7))
    tail_bound += 3*radius**7/(1-radius)
    assert tail_bound == F(17, 64)
    assert tail_bound < 3*radius-1
    reciprocal_bound = 1/(3*radius-1-tail_bound)
    remainder_bound = reciprocal_bound + 2
    assert reciprocal_bound == F(64, 15)
    assert remainder_bound == F(94, 15)

    ratios = []
    for index in range(1, 13):
        common = gcd(returns[index+1], returns[index])
        ratios.append({"n": index, "raw_numerator": returns[index+1],
                       "raw_denominator": returns[index], "gcd": common,
                       "reduced_ratio": str(F(returns[index+1], returns[index]))})

    return {
        "pi_certificate": "Machin identity; 100 alternating terms per arctangent; exact Fraction bounds",
        "boundary_digits": ''.join(map(str, boundary)),
        "four_units_pi_base": '10.' + ''.join(map(str, four_tail)) + '...',
        "return_counts_A_0_through_A_20": returns[:21],
        "graph_vs_renewal_checked_lengths": 30,
        "formal_inverse_checked_coefficients": 40,
        "rouche_circle_radius": str(radius),
        "nonlinear_tail_bound": str(tail_bound),
        "analytic_remainder_bound": str(remainder_bound),
        "adjacent_count_ratios": ratios,
        "scope": "Finite exact certificates; no claim to prove P2848 or infinite randomness."
    }


if __name__ == '__main__':
    print(json.dumps(main(), ensure_ascii=False, indent=2))
