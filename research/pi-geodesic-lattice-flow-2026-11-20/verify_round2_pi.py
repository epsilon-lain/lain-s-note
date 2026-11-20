#!/usr/bin/env python3
"""Certify the round-two pi prefix and its active nodes with exact arithmetic."""

from fractions import Fraction


def arctan_bounds(base: int, scale: int) -> tuple[int, int]:
    """Enclose scale * arctan(1/base) using directed rounding."""
    lower = upper = 0
    index = 0
    power = base
    while True:
        floor_term, remainder = divmod(scale, (2 * index + 1) * power)
        ceil_term = floor_term + bool(remainder)
        if floor_term == 0:
            # The alternating tail lies between zero and its first term.
            if index % 2 == 0:
                upper += ceil_term
            else:
                lower -= ceil_term
            return lower, upper
        if index % 2 == 0:
            lower += floor_term
            upper += ceil_term
        else:
            lower -= ceil_term
            upper -= floor_term
        power *= base**2
        index += 1


def pi_bounds(scale: int) -> tuple[Fraction, Fraction]:
    lower5, upper5 = arctan_bounds(5, scale)
    lower239, upper239 = arctan_bounds(239, scale)
    # Machin's identity: pi = 16 arctan(1/5) - 4 arctan(1/239).
    lower = Fraction(16 * lower5 - 4 * upper239, scale)
    upper = Fraction(16 * upper5 - 4 * lower239, scale)
    return lower, upper


def common_cf(lower: Fraction, upper: Fraction) -> list[int]:
    terms = []
    while True:
        digit = lower.numerator // lower.denominator
        if digit != upper.numerator // upper.denominator or lower == digit:
            return terms
        terms.append(digit)
        lower, upper = 1 / (upper - digit), 1 / (lower - digit)


def absolute_bounds(lower: Fraction, upper: Fraction) -> tuple[Fraction, Fraction]:
    if lower >= 0:
        return lower, upper
    if upper <= 0:
        return -upper, -lower
    return Fraction(0), max(-lower, upper)


def main() -> None:
    lower, upper = pi_bounds(10**3010)
    terms = common_cf(lower, upper)
    assert len(terms) == 2929
    assert max(terms) == 20776 and terms.index(max(terms)) == 431

    # Finer enclosure certifies the inequalities even at the final denominator.
    lower, upper = pi_bounds(10**3040)
    assert common_cf(lower, upper)[:len(terms)] == terms
    half = Fraction(1, 2)
    records = []
    active = []
    p, q, last_p, last_q = 1, 0, 0, 1
    for index, digit in enumerate(terms):
        p, last_p = digit * p + last_p, p
        q, last_q = digit * q + last_q, q
        value_lower, value_upper = absolute_bounds(q * (q * lower - p), q * (q * upper - p))
        if not records or value_upper < records[-1][2]:
            records.append((index, q, value_lower, value_upper))
        else:
            assert value_lower > records[-1][3], "Record ordering is unresolved"
        if value_upper < half:
            active.append((index, p, q))
        else:
            assert value_lower >= half, "Active-node classification is unresolved"

    assert [record[0] for record in records] == [0, 1, 3, 306, 430]
    best_index, best_q, best_lower, best_upper = records[-1]
    assert len(str(best_q)) == 216 and len(str(q)) == 1503
    # Verify exact rounded bounds before displaying short decimals.
    assert Fraction("0.0000481305612462893") < best_lower
    assert best_upper < Fraction("0.0000481305612462895")

    largest_ratio = Fraction(0)
    largest_pair = None
    for (i, p, q), (j, next_p, next_q) in zip(active, active[1:]):
        assert abs(p * next_q - next_p * q) == 1, (i, j)
        ratio = Fraction(next_q, q)
        if ratio > largest_ratio:
            largest_ratio, largest_pair = ratio, (i, j)

    print(f"Certified partial quotients: {len(terms)}; maximum: a_431 = {max(terms)}")
    print(f"Defect record indices: {[record[0] for record in records]}")
    print(f"Finite minimum at n={best_index}: 0.0000481305612462893 < defect < 0.0000481305612462895")
    print(f"Active nodes, with q²|pi - p/q| < 1/2: {len(active)}")
    print(f"Consecutive active-node pairs with absolute determinant 1: {len(active) - 1}")
    print(f"Largest consecutive active denominator ratio: {float(largest_ratio):.12f}, indices {largest_pair}")


if __name__ == "__main__":
    main()
