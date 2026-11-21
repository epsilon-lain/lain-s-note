"""Exact certificates for the historical and 2025 continued-fraction records."""
from fractions import Fraction as F

from verify_two_window_profiles import abs_interval, pi_interval, sinc_polynomial


RECORD_DIGIT = 845_188_380_288
RECORD_CONVERGENT = 383_059_989_912


def convergents_from_interval(lo, hi, length):
    p_prev, p = 0, 1
    q_prev, q = 1, 0
    for _ in range(length):
        digit = lo.__floor__()
        assert digit == hi.__floor__() and digit < lo
        p_prev, p = p, digit * p + p_prev
        q_prev, q = q, digit * q + q_prev
        yield digit, p, q
        lo, hi = 1 / (hi - digit), 1 / (lo - digit)


def verify_half_height(q, w_lo, w_hi, pass_h, fail_h, background, pi_lo):
    """Use the principal lobe of sinc and a bound for all later lobes."""
    t_pass = pass_h * q**2 * w_hi
    t_fail = fail_h * q**2 * w_lo
    assert 0 < t_pass < t_fail < 2
    sinc_lower = sinc_polynomial(t_pass, 13)
    sinc_upper = sinc_polynomial(t_fail, 12)
    denominator_lower = 1 - w_hi**2 / 6
    assert sinc_lower > 0 and sinc_lower**2 > F(1, 2)
    assert denominator_lower > 0 and sinc_upper > 0
    tail = F(background, fail_h**2)
    assert (sinc_upper / denominator_lower)**2 + tail < F(1, 2)
    assert 1 / (pi_lo * denominator_lower)**2 + tail < F(1, 2)


def verify_first_record():
    raw_lo, raw_hi = pi_interval(360)
    scale = 10**440
    lo = F((raw_lo * scale).__floor__(), scale)
    hi = F((raw_hi * scale).__ceil__(), scale)
    convergents = list(convergents_from_interval(lo, hi, 309))
    digits = [digit for digit, _, _ in convergents]
    assert digits[4] == 292 and digits[307:309] == [436, 8]
    assert max(digits[:4] + digits[5:307]) == 161
    _, p, q = convergents[306]
    _, p_prev, q_prev = convergents[305]
    assert q * p_prev - p * q_prev == 1 and F(p, q) < lo
    assert len(str(p)) == len(str(q)) == 155
    lower = q**2 * (lo - F(p, q))
    upper = q**2 * (hi - F(p, q))
    assert F(2288738616993029, 10**18) < lower < upper < F(2288738616993030, 10**18)

    # Compare all earlier convergents via the next-partial-quotient bound.
    old_lo, old_hi = abs_interval(lo - F(355, 113), hi - F(355, 113))
    assert upper < 113**2 * old_lo < 113**2 * old_hi < F(1, 163)
    old_w_lo, old_w_hi = (113 * old_lo / 355, 113 * old_hi / 355)
    assert 355**2 * old_w_hi / lo < F(1, 8)
    assert 8 * hi**2 / 3 < 27
    verify_half_height(355, old_w_lo, old_w_hi, 130, 131, 27, lo)

    window = 2**513
    w_lo, w_hi = q * lo / p - 1, q * hi / p - 1
    assert window <= p < 2 * window and window**2 * w_hi / lo < F(1, 800)
    assert (F(800, 199)**2 * hi**2 / 12) < 14
    verify_half_height(window, w_lo, w_hi, 605, 606, 14, lo)
    verify_half_height(p, w_lo, w_hi, 193, 194, 27, lo)

    # In the later window the nearest fraction is an intermediate convergent.
    window = 2**521
    k = (2 * window - 1 - p_prev) // p
    candidate_denominator = k * p + p_prev
    assert k == 288 and window <= candidate_denominator < 2 * window
    distance_lo = 1 / F(p * candidate_denominator) - (F(q, p) - 1 / hi)
    distance_hi = 1 / F(p * candidate_denominator) - (F(q, p) - 1 / lo)
    assert F(2456903536030, 10**11) < window**2 * distance_lo
    assert window**2 * distance_hi < F(2456903536031, 10**11)
    # Fractions on the opposite side obey the determinant lower bound.
    assert distance_hi < 1 / F(2 * window * p) + F(q, p) - 1 / lo
    assert candidate_denominator > window  # determinant at least two is worse
    print("First replacement of 355/113, convergent n=306:")
    print(f"p = {p}\nq = {q}")
    print("Last half-height: old anchored 130, new anchored 193, new dyadic 605")
    print("Q=2^521: k=288, nearest normalized distance in (24.56903536030,24.56903536031)")


def verify_giant_record_consequences():
    a = RECORD_DIGIT
    lo, hi = pi_interval(24)
    # For Q=r, only Q^2*w and an upper bound for w enter the proof.
    # Pull the worst bounds from 1/(a+2) < r^2*abs(1/pi-x) < 1/a.
    verify_half_height(1, lo / (a + 2), hi / a, 374_373_210_116, 374_373_210_118, 27, lo)
    assert a % 4 == 0
    assert a - 2 > 2**39
    # The latter gives 37 fully enclosed consecutive dyadic denominator windows.
    depth = F(a, 18)
    assert a // 18 == 46_954_910_016 and 2 / depth + 1 / (2 * depth**2) < F(43, 10**12)
    assert F(a, 16) == 52_824_273_768
    growth = F(809, 500)
    assert growth**2 < growth + 1 and growth**10000 > 10**2089
    assert 2089 * ((RECORD_CONVERGENT - 1) // 10000) > 80_000_000_000
    print("Giant record, source digit assumed: last anchored half-height 374373210116 or 374373210117")
    print("A later dyadic window has D > A/18 and wave < 4.3e-11/H^2")


if __name__ == "__main__":
    verify_first_record()
    verify_giant_record_consequences()
