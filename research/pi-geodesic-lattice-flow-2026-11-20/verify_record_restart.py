"""Reproduce the true successor of historical seed 355/113 in determinant coordinates."""
from fractions import Fraction
from math import lcm

from verify_historical_records import common_convergents, pi_bounds


def first_record(lower, upper, p, q, previous_p, previous_q, digit):
    if abs(p * previous_q - q * previous_p) != 1:
        raise ValueError("the two vectors must form a unimodular basis")

    def read_tail(value):
        return (previous_p - value * previous_q) / (value * q - p) - digit

    if (lower - Fraction(p, q)) * (upper - Fraction(p, q)) <= 0:
        raise ValueError("the interval must avoid the seed pole")
    tail_low, tail_high = sorted(map(read_tail, (lower, upper)))
    if not 0 < tail_low < tail_high < 1:
        raise ValueError("the seed and next partial quotient are not certified")
    common_denominator = lcm(tail_low.denominator, tail_high.denominator)
    x = Fraction(previous_q, q)
    tail_convergents = common_convergents(int(tail_low * common_denominator),
                                        int(tail_high * common_denominator),
                                        common_denominator)
    for index, k, v in tail_convergents:
        center = Fraction(k, v)
        radius = 1 / (v * (v * (digit + x) + k))
        if tail_low >= center + radius or tail_high <= center - radius:
            continue
        if tail_low <= center - radius or tail_high >= center + radius:
            raise ValueError("increase the precision before deciding the first record")
        u = digit * v + k
        return index, k, v, u * p + v * previous_p, u * q + v * previous_q
    raise ValueError("no record in the certified tail prefix")


if __name__ == "__main__":
    scale = 10**420
    lower, upper = (Fraction(x, scale) for x in pi_bounds(scale))
    index, k, v, p, q = first_record(lower, upper, 355, 113, 333, 106, 292)
    pi_convergents = common_convergents(int(lower * scale), int(upper * scale), scale)
    expected = next((cp, cq) for i, cp, cq in pi_convergents if i == 306)
    assert index == 302 and (p, q) == expected
    print(f"Confirmed the first strict successor after {index + 1} tail convergents.")
    print(f"v={v}\nk={k}\np={p}\nq={q}")
