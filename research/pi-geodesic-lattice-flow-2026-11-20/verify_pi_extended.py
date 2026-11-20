#!/usr/bin/env python3
"""Certify a finite range of pi approximants using integer arithmetic."""

from decimal import Decimal, localcontext


def arctangent_bounds(scale, base):
    """Return rational endpoints with common denominator scale."""
    power = scale // base
    total = 0
    count = 0
    while power:
        term = power // (2 * count + 1)
        total += term if count % 2 == 0 else -term
        count += 1
        power //= base * base

    # The rounding errors are each smaller than one unit. The first
    # omitted term bounds the alternating tail by less than one unit.
    return total - count // 2 - 1, total + (count + 1) // 2 + 1


def common_prefix(lower, upper, scale):
    left_num, left_den = lower, scale
    right_num, right_den = upper, scale
    digits = []
    while left_den and right_den:
        digit = left_num // left_den
        if digit != right_num // right_den:
            break
        digits.append(digit)
        left_rem = left_num - digit * left_den
        right_rem = right_num - digit * right_den
        left_num, right_num = right_den, left_den
        left_den, right_den = right_rem, left_rem
    return digits


def defect_bounds(p, q, lower, upper, scale):
    """Bound scale times the weighted error q * abs(q * pi - p)."""
    left = q * (q * lower - p * scale)
    right = q * (q * upper - p * scale)
    low = 0 if left <= 0 <= right else min(abs(left), abs(right))
    return low, max(abs(left), abs(right))


def in_critical_band(digits, index, p, q, q_prev, lower, upper, scale):
    """Return whether the convergent defect is strictly below 1/2."""
    tail_digits = digits[index + 1 : index + 13]
    if tail_digits and tail_digits[0] >= 2:
        return True
    if tail_digits:
        tail_num_prev, tail_num = 0, 1
        tail_den_prev, tail_den = 1, 0
        for digit in tail_digits:
            tail_num_prev, tail_num = tail_num, digit * tail_num + tail_num_prev
            tail_den_prev, tail_den = tail_den, digit * tail_den + tail_den_prev

        # The complete quotient lies strictly between the last
        # convergent of the tail and its mediant with the predecessor.
        # Defect < 1/2 iff that quotient plus q_prev/q exceeds 2.
        base = q_prev - 2 * q
        first_sign = q * tail_num + base * tail_den
        second_sign = first_sign + q * tail_num_prev + base * tail_den_prev
        if min(first_sign, second_sign) >= 0:
            return True
        if max(first_sign, second_sign) <= 0:
            return False

    low, high = defect_bounds(p, q, lower, upper, scale)
    if 2 * high < scale:
        return True
    if 2 * low >= scale:
        return False
    raise ArithmeticError(f"Cannot classify critical-band convergent {index}")


def improves_ratio(numerator, denominator, best):
    return best is None or numerator * best[1] > best[0] * denominator


def decimal_ratio(numerator, denominator):
    with localcontext() as context:
        context.prec = 30
        return Decimal(numerator) / Decimal(denominator)


def main():
    scale = 10**30000
    low5, high5 = arctangent_bounds(scale, 5)
    low239, high239 = arctangent_bounds(scale, 239)
    lower, upper = 16 * low5 - 4 * high239, 16 * high5 - 4 * low239
    digits = common_prefix(lower, upper, scale)
    assert len(digits) == 29166

    p_prev, p = 1, digits[0]
    q_prev, q = 0, 1
    best_low = best_high = scale
    records = []
    previous_active = None
    max_active_ratio = None
    max_active_step = None
    active_count = 0

    for index in range(len(digits)):
        if index:
            digit = digits[index]
            p_prev, p = p, digit * p + p_prev
            q_prev, q = q, digit * q + q_prev

        if index + 1 == len(digits) or best_high * (digits[index + 1] + 2) > scale:
            low, high = defect_bounds(p, q, lower, upper, scale)
            if high < best_low:
                best_low, best_high = low, high
                records.append((index, q, low, high))
            else:
                assert low > best_high, index

        if in_critical_band(digits, index, p, q, q_prev, lower, upper, scale):
            active_count += 1
            if previous_active is not None:
                previous_index, previous_q = previous_active
                if improves_ratio(q, previous_q, max_active_ratio):
                    max_active_ratio = (q, previous_q, previous_index, index)
            previous_active = index, q
            if index + 1 < len(digits):
                next_q = digits[index + 1] * q + q_prev
                if improves_ratio(next_q, q, max_active_step):
                    max_active_step = next_q, q, index, index + 1

    assert [record[0] for record in records] == [0, 1, 3, 306, 430, 28420]
    index, record_q, record_low, record_high = records[-1]
    assert digits[index + 1] == 78629
    assert record_low * 10**32 > 1271785737275745660358446939 * scale
    assert record_high * 10**32 < 1271785737275745660358446940 * scale
    assert 10**14592 < record_q < 2 * 10**14592

    # The first unclassified convergent denominator is at least q+q_prev.
    # A distinct fraction beating the record must be a convergent by
    # Legendre's theorem, after reducing to lowest terms if necessary.
    assert q + q_prev > 10**14997
    final_low, final_high = defect_bounds(p, q, lower, upper, scale)

    print(f"Machin interval width: {upper - lower}/10^30000")
    print(f"Common continued-fraction terms: {len(digits)}")
    print("Record indices:", [record[0] for record in records])
    print(f"Last record is at {index}, with next partial quotient {digits[index + 1]}")
    print("Last record defect:", decimal_ratio(record_low, scale))
    print("Final convergent defect:", decimal_ratio(final_low, scale),
          decimal_ratio(final_high, scale))
    print("Last record remains minimal up to Q=10^14997")
    print(f"Convergents in the strict half-band: {active_count}")
    for label, (num, den, first, second) in [
        ("Successive half-band nodes", max_active_ratio),
        ("Half-band node to next convergent", max_active_step),
    ]:
        print(f"Largest ratio ({label}) has indices {first}->{second}:",
              decimal_ratio(num, den))


if __name__ == "__main__":
    main()
