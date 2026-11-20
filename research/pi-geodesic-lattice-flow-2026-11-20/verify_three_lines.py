#!/usr/bin/env python3
"""Certify a single curve's three priority lines per dyadic denominator band."""
from fractions import Fraction
from itertools import islice
from math import gcd

HALF = Fraction(1, 2)
FRAME_LIMIT = Fraction(5, 4)


def arctan_bounds(base, terms):
    partial = sum(
        Fraction((-1) ** j, (2 * j + 1) * base ** (2 * j + 1))
        for j in range(terms)
    )
    next_partial = partial + Fraction(
        (-1) ** terms, (2 * terms + 1) * base ** (2 * terms + 1)
    )
    return min(partial, next_partial), max(partial, next_partial)


def pi_bounds():
    a, b = arctan_bounds(5, 80)
    c, d = arctan_bounds(239, 25)
    return 16 * a - 4 * d, 16 * b - 4 * c


def convergents(rational):
    a, b = rational.as_integer_ratio()
    p_prev, p = 0, 1
    q_prev, q = 1, 0
    while b:
        digit, remainder = divmod(a, b)
        p_prev, p = p, digit * p + p_prev
        q_prev, q = q, digit * q + q_prev
        yield p, q
        a, b = b, remainder


def certified_frame(pool, low, high, scale, ratio=2,
                    threshold=HALF, frame_limit=FRAME_LIMIT):
    for a, b in pool:
        if threshold * b >= scale * frame_limit:
            break
        error = max(abs(b * low - a), abs(b * high - a))
        if max(threshold * b / scale + scale * error,
               threshold * b / (ratio * scale) + ratio * scale * error) < frame_limit:
            assert gcd(a, b) == 1
            return a, b
    raise AssertionError(f"No certified frame at scale {scale}")


def endpoints(a, b, k, first, last):
    """Return first denominator, last progression index."""
    residue = (-k * pow(a, -1, b)) % b if b > 1 else 0
    start = first + (residue - first) % b
    return start, (last - start) // b


def index_before_zero(a, b, k, start, slope):
    root_index = (k / (slope * b - a) - start) / b
    return root_index.numerator // root_index.denominator


def line_minimum_candidates(a, b, k, first, last, slope):
    start, count = endpoints(a, b, k, first, last)
    if count < 0:
        return set()
    # q * |q*slope - p| = q * |(slope-a/b)*q - k/b|.
    candidates = {0, count}
    if slope * b != a:
        before = index_before_zero(a, b, k, start, slope)
        candidates.update((before, before + 1))
    results = set()
    for index in candidates:
        if 0 <= index <= count:
            q = start + index * b
            p, remainder = divmod(a * q + k, b)
            assert not remainder
            results.add((p, q))
    return results


def half_threshold_error(low, high, p, q):
    left, right = q * low - p, q * high - p
    return (
        q * min(abs(left), abs(right)) if left * right > 0 else Fraction(0),
        q * max(abs(left), abs(right)),
    )


def brute_candidates(slope, scale):
    good = set()
    for q in range(scale, 2 * scale):
        p = round(q * slope)
        if q * abs(q * slope - p) <= HALF:
            good.add((p, q))
    return good


def enumerate_three_lines(a, b, slope, scale):
    good = set()
    for k in (-1, 0, 1):
        first, count = endpoints(a, b, k, scale, 2 * scale - 1)
        for index in range(count + 1):
            q = first + index * b
            p = (a * q + k) // b
            if q * abs(q * slope - p) <= HALF:
                good.add((p, q))
    return good


def validate_clock(limit):
    jobs = []
    for diagonal in range(limit):
        for i in range(diagonal + 1):
            m = diagonal - i
            k = (2 + m // 2) * (1 if m % 2 == 0 else -1)
            jobs.append((i, k))
            assert len(jobs) - 1 == diagonal * (diagonal + 1) // 2 + i
    assert len(jobs) == len(set(jobs))
    return len(jobs)


def accelerated_candidates(pool, low, high, cap):
    threshold, ratio = Fraction(1, 280), 256
    assert (ratio + 1) * threshold < 1
    assert cap >= 113
    assert half_threshold_error(low, high, 355, 113)[1] < threshold
    scale = 1
    bands = 0
    candidates = set()
    while scale <= cap:
        a, b = certified_frame(pool, low, high, scale, ratio, threshold, 1)
        first, count = endpoints(a, b, 0, scale, min(cap, ratio * scale - 1))
        if count >= 0:
            candidates.add((a * (first // b), first))
        scale *= ratio
        bands += 1
    return bands, candidates


def main():
    low, high = pi_bounds()
    slope = (low + high) / 2
    pool = list(islice(convergents(slope), 140))
    checked = 0
    candidates = set()
    for power in range(60):
        scale = 2**power
        a, b = certified_frame(pool, low, high, scale)
        for k in (-1, 0, 1):
            at_low = line_minimum_candidates(a, b, k, scale, 2 * scale - 1, low)
            at_high = line_minimum_candidates(a, b, k, scale, 2 * scale - 1, high)
            assert at_low == at_high
            # The signed drift is affine in the true slope and has constant sign.
            assert (b * low - a) * (b * high - a) > 0
            # q at a zero crossing moves monotonically as the slope varies.
            start, count = endpoints(a, b, k, scale, 2 * scale - 1)
            if count >= 0:
                assert index_before_zero(a, b, k, start, low) == (
                    index_before_zero(a, b, k, start, high)
                )
            candidates.update(at_low)
        if power < 15:
            direct = brute_candidates(slope, scale)
            on_lines = enumerate_three_lines(a, b, slope, scale)
            assert direct == on_lines
            # Enclosure cannot change a nearest-integer or a threshold decision.
            for q in range(scale, 2 * scale):
                p = round(q * slope)
                assert round(q * low) == round(q * high) == p
                left, right = half_threshold_error(low, high, p, q)
                assert right < HALF or left > HALF
            checked += scale

    winner = min(candidates, key=lambda pq: half_threshold_error(low, high, *pq)[1])
    lo_value, hi_value = half_threshold_error(low, high, *winner)
    assert hi_value < HALF
    assert all(hi_value < half_threshold_error(low, high, *pq)[0]
               for pq in candidates if pq != winner)
    print(f"Independent all-denominator comparison through q={checked}: exact match")
    print(f"Total three-line work through q={2**60 - 1}: "
          f"180 lines, {len(candidates)} distinct candidate pairs")
    print(f"Certified minimizer: {winner[0]}/{winner[1]}, "
          f"normalized error {float((lo_value + hi_value) / 2):.15g}")
    bands, fast = accelerated_candidates(pool, low, high, 2**60 - 1)
    assert winner in fast
    assert all(hi_value < half_threshold_error(low, high, *pq)[0]
               for pq in fast if pq != winner)
    print(f"Adaptive clock at the same cap: {bands} bands, {len(fast)} tested values")
    print(f"Deferred single-clock jobs checked: {validate_clock(64)}")


if __name__ == "__main__":
    main()
