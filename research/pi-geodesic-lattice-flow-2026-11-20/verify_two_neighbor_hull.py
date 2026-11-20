"""Certify the local selector for pi against all past and finite future candidates."""
from argparse import ArgumentParser
from fractions import Fraction


def atan_bounds(base, scale):
    lower = upper = 0
    power = base
    j = 0
    while (2 * j + 1) * power <= scale:
        denominator = (2 * j + 1) * power
        low_term, remainder = divmod(scale, denominator)
        high_term = low_term + bool(remainder)
        if j % 2 == 0:
            lower += low_term
            upper += high_term
        else:
            lower -= high_term
            upper -= low_term
        j += 1
        power *= base**2
    # The next omitted alternating term has size less than one grid unit.
    if j % 2 == 0:
        upper += 1
    else:
        lower -= 1
    return lower, upper


def pi_bounds(method, scale):
    if method == "machin":
        lo5, hi5 = atan_bounds(5, scale)
        lo239, hi239 = atan_bounds(239, scale)
        return 16 * lo5 - 4 * hi239, 16 * hi5 - 4 * lo239
    lo2, hi2 = atan_bounds(2, scale)
    lo3, hi3 = atan_bounds(3, scale)
    return 4 * (lo2 + lo3), 4 * (hi2 + hi3)


def convergents(lo, hi, scale, count):
    ln, ld = lo, scale
    hn, hd = hi, scale
    p_prev, q_prev, p, q = 0, 1, 1, 0
    points = [(1, 0)]  # Index 0 represents the vertical vector.
    for _ in range(count):
        digit = ln // ld
        assert digit == hn // hd, "increase --digits to certify the prefix"
        p_prev, p = p, digit * p + p_prev
        q_prev, q = q, digit * q + q_prev
        points.append((p, q))
        ln, ld, hn, hd = hd, hn % hd, ld, ln % ld
    return points


def error_bounds(p, q, lo, hi, scale):
    left, right = q * lo - p * scale, q * hi - p * scale
    if right < 0:
        return -right, -left
    assert left > 0, "increase --digits to separate a convergent error from zero"
    return left, right


def crossing(points, errors, i, j):
    # All crossing scales are divided by the common factor scale**2.
    delta_q = points[j][1] ** 2 - points[i][1] ** 2
    ei_lo, ei_hi = errors[i]
    ej_lo, ej_hi = errors[j]
    assert ei_lo > ej_hi, "increase --digits to order the errors"
    return (Fraction(delta_q, ei_hi**2 - ej_lo**2),
            Fraction(delta_q, ei_lo**2 - ej_hi**2))


def certified_extreme(candidates, maximize):
    compare = max if maximize else min
    index, (lower, upper) = compare(candidates, key=lambda item: item[1][0])
    for rival, interval in candidates:
        if rival == index:
            continue
        if maximize:
            assert lower > interval[1], "past maximum not separated"
        else:
            assert upper < interval[0], "future minimum not separated"
    return index, (lower, upper)


def bezout_candidates(p, q, lo, hi, scale):
    candidates = []
    for determinant in (-1, 1):
        q0 = (-determinant * pow(p, -1, q)) % q if q > 1 else 0
        p0 = (determinant + p * q0) // q
        centers = [Fraction(q0 * endpoint - p0 * scale, p * scale - q * endpoint)
                   for endpoint in (lo, hi)]
        left, right = sorted(centers)
        j = left.numerator // left.denominator
        assert j == right.numerator // right.denominator, "center not separated"
        for offset in (j, j + 1):
            next_q = q0 + offset * q
            if next_q > q:
                candidates.append((p0 + offset * p, next_q))
    return candidates


def certify(method, scale, count):
    lo, hi = pi_bounds(method, scale)
    points = convergents(lo, hi, scale, count + 12)
    errors = [error_bounds(p, q, lo, hi, scale) for p, q in points]
    active_count = skipped_count = past_two = future_two = 0
    examples = {}
    for i in range(1, count + 1):
        past = [(j, crossing(points, errors, j, i)) for j in range(i)]
        future = [(j, crossing(points, errors, i, j)) for j in range(i + 1, len(points))]
        lower_idx, lower = certified_extreme(past, maximize=True)
        upper_idx, upper = certified_extreme(future, maximize=False)
        unseen_lower = Fraction(points[-1][1]**2 - points[i][1]**2, errors[i][1]**2)
        assert unseen_lower > upper[1]
        time_start = Fraction(1, scale**2)
        assert lower[0] > time_start
        if lower[1] < upper[0]:
            active = True
        else:
            assert lower[0] > upper[1], "global selector margin not separated"
            active = False
        assert i - lower_idx <= 2 and upper_idx - i <= 2
        p, q = points[i]
        previous_q = points[i - 1][1]
        threshold = Fraction(q + 2 * previous_q, 2 * (q**2 + q * previous_q + previous_q**2))
        local_lo, local_hi = (Fraction(endpoint, scale) for endpoint in errors[i])
        if active:
            assert local_hi < threshold
            active_count += 1
            candidates = bezout_candidates(p, q, lo, hi, scale)
            candidate_errors = [error_bounds(*candidate, lo, hi, scale)
                                for candidate in candidates]
            crossings = [(j, crossing([(p, q), candidate], [errors[i], error], 0, 1))
                         for j, (candidate, error) in enumerate(zip(candidates, candidate_errors))]
            successor, _ = certified_extreme(crossings, maximize=False)
            assert candidates[successor] == points[upper_idx]
        else:
            assert local_lo > threshold
            skipped_count += 1
        past_two += i - lower_idx == 2
        future_two += upper_idx - i == 2
        if q in {7, 106, 113, 33102, 33215, 99532, 265381}:
            examples[q] = (p, points[lower_idx][1], points[upper_idx][1], active)
    i = points.index((208341, 66317))
    lo_eta, hi_eta = (Fraction(endpoint * 66317, scale) for endpoint in errors[i])
    assert Fraction("0.53811") < lo_eta and hi_eta < Fraction("0.53813")
    assert 33215 * 312689 - 99532 * 104348 == -1
    return (active_count, skipped_count, past_two, future_two, examples)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, default=128)
    parser.add_argument("--digits", type=int, default=300)
    args = parser.parse_args()
    if args.count < 1 or args.digits < 8:
        parser.error("expected a positive count and at least eight digits")
    scale = 10**args.digits
    machin = certify("machin", scale, args.count)
    two_angle = certify("two-angle", scale, args.count)
    assert machin == two_angle
    active, skipped, past_two, future_two, examples = machin
    print(f"Count {args.count}: {active} active, {skipped} skipped.")
    print(f"Second preceding point wins {past_two} times; second following point wins {future_two} times.")
    for q, (p, before, after, active) in examples.items():
        print(f"{p}/{q}: past q={before}, future q={after}, active={active}")
