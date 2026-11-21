#!/usr/bin/env python3
"""Certify historical normalized-error records from a rational enclosure of pi.

All comparisons use integers. The decimal display serves only to present
the relative sizes. A convergent is indexed from n=0 with a_0=3.
"""

from decimal import Decimal, localcontext


def arctangent_bounds(inverse, scale):
    """Return integer lower and upper bounds on scale * arctan(1/inverse)."""
    lower = upper = 0
    index = 0
    power = inverse
    sign = 1
    while True:
        floor, remainder = divmod(scale, (2 * index + 1) * power)
        ceiling = floor + bool(remainder)
        if floor == 0:
            # The alternating tail has the sign of this first omitted
            # term and magnitude at most its ceiling.
            return (lower, upper + ceiling) if sign > 0 else (lower - ceiling, upper)
        if sign > 0:
            lower += floor
            upper += ceiling
        else:
            lower -= ceiling
            upper -= floor
        index += 1
        power *= inverse * inverse
        sign = -sign


def pi_bounds(scale):
    """Use Machin's identity to bound scale * pi."""
    lo5, hi5 = arctangent_bounds(5, scale)
    lo239, hi239 = arctangent_bounds(239, scale)
    return 16 * lo5 - 4 * hi239, 16 * hi5 - 4 * lo239


def common_convergents(lower, upper, scale):
    """Yield the convergents prescribed by every irrational in the interval."""
    left_num, left_den = lower, scale
    right_num, right_den = upper, scale
    p_old, p = 0, 1
    q_old, q = 1, 0
    index = 0
    while True:
        left_digit, left_remainder = divmod(left_num, left_den)
        right_digit, right_remainder = divmod(right_num, right_den)
        if left_digit != right_digit:
            return
        p_old, p = p, left_digit * p + p_old
        q_old, q = q, left_digit * q + q_old
        yield index, p, q
        if left_remainder == 0:
            return
        left_num, left_den, right_num, right_den = (
            right_den, right_remainder, left_den, left_remainder
        )
        index += 1


def error_interval(p, q, lower, upper, scale):
    residual_lower = q * (q * lower - p * scale)
    residual_upper = q * (q * upper - p * scale)
    if residual_lower <= 0 <= residual_upper:
        return 0, max(-residual_lower, residual_upper)
    return sorted((abs(residual_lower), abs(residual_upper)))


def certified_records(lower, upper, scale):
    best_interval = None
    records = []
    count = 0
    for index, p, q in common_convergents(lower, upper, scale):
        count += 1
        lo, hi = error_interval(p, q, lower, upper, scale)
        if best_interval is None or hi < best_interval[0]:
            best_interval = (lo, hi)
            records.append((index, p, q, lo, hi))
        elif lo < best_interval[1]:
            raise ValueError(f"The error intervals for record comparison overlap at n={index}.")
    return count, records


def main():
    scale = 10**2100
    lower, upper = pi_bounds(scale)
    count, records = certified_records(lower, upper, scale)
    assert count == 2048
    assert [r[0] for r in records] == [0, 1, 3, 306, 430]
    assert records[2][1:3] == (355, 113)
    print(f"Certified {count} consecutive convergents from n=0.")
    with localcontext() as context:
        context.prec = 19
        for index, p, q, lo, hi in records:
            # This display is approximate. The interval and strict
            # comparison above use exact integers.
            estimated_error = Decimal(lo + hi) / (2 * scale)
            print(f"n={index}  E≈{estimated_error}\n  p={p}\n  q={q}")


if __name__ == "__main__":
    main()
