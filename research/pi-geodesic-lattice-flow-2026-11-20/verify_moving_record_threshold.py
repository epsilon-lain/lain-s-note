"""Certify a moving-observer threshold using the published giant CF digit.

Place next to verify_two_window_profiles.py and run with Python 3.
"""
from fractions import Fraction as F

from verify_two_window_profiles import pi_interval, sinc_polynomial


def verify_moving_record_threshold():
    digit = 845_188_380_288
    mobility = F(1, 50)
    last = 391_276_977_722
    # c solves sinc(pi*c)^2 = 1/2. Bound it by monotonicity on the first lobe.
    c_lo = F(44294647068945233, 10**17)
    c_hi = F(44294647068945235, 10**17)
    pi_lo, pi_hi = pi_interval(24)
    sinc_lo = sinc_polynomial(c_lo * pi_hi, 15)
    sinc_hi = sinc_polynomial(c_hi * pi_lo, 16)
    assert 0 < sinc_lo < 1 and 2*sinc_lo**2 > 1
    assert 0 < sinc_hi < 1 and 2*sinc_hi**2 < 1

    # D=r^2*dist(1/pi, S_r) lies strictly between 1/(digit+2) and 1/digit.
    # Every H up to 'last' reaches a half-peak (inner arcs).
    assert F(last, digit) < mobility + c_lo
    # For H>=last+1, the LHS increases and the last term of the RHS decreases.
    assert F(last+1, digit+2) > mobility + c_hi + F(31, (last+1)**2)
    print("A =", digit, "at Q=r, mobility = 1/50")
    print("Last H admitting a mobile half-height =", last)


if __name__ == "__main__":
    verify_moving_record_threshold()
