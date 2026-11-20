"""An exact rational enclosure of the weighted theta fourth moment at T=10^-7.

See FOUR_BUNDLE_REVIEW.zh.md for the Poisson formula and the Fourier tail.
All arithmetic in this script is integer or Fraction arithmetic.
"""
from fractions import Fraction as F
from math import factorial


def atan_bounds(x, pairs):
    """Alternating Taylor bounds for atan(x), 0 <= x < 1."""
    assert 0 <= x < 1
    terms = 2 * pairs
    lo = sum(((-1)**j * x**(2*j + 1) / (2*j + 1)
              for j in range(terms)), F(0))
    return lo, lo + x**(2*terms + 1) / (2*terms + 1)


def round_out(lo, hi, places=55):
    scale = 10**places
    return F((lo * scale).__floor__(), scale), F((hi * scale).__ceil__(), scale)


def verify():
    a_lo, a_hi = atan_bounds(F(1, 5), 30)
    b_lo, b_hi = atan_bounds(F(1, 239), 12)
    pi_lo, pi_hi = round_out(16*a_lo - 4*b_hi, 16*a_hi - 4*b_lo)
    p, q, T = 113, 355, F(1, 10**7)
    d_lo, d_hi = round_out(2 - 2*pi_hi*p/q, 2 - 2*pi_lo*p/q)
    assert 0 < d_lo <= d_hi

    # t/(t^2+d^2) has a single maximum, so check both interval endpoints.
    min_real = min(pi_lo**2 * t / (4*q*q*(t*t + d_hi*d_hi))
                   for t in (T, 2*T))
    assert min_real > 50
    assert sum(F(50**k, factorial(k)) for k in range(101)) > 10**21
    # The constant Fourier coefficient, multiplied by sqrt(q), has modulus one.
    # Twice the uniform bound 3/sqrt(2) is less than 4.5. For m>=1,
    # m^2 >= 1+3(m-1); sum the resulting geometric progression.
    exp_bound = F(1, 10**21)
    tail = F(9, 2)*exp_bound/(1-exp_bound**3)

    x_lo, x_hi = round_out(T*d_lo/(d_hi*d_hi + 2*T*T),
                           T*d_hi/(d_lo*d_lo + 2*T*T))
    angle_lo = atan_bounds(x_lo, 22)[0]
    angle_hi = atan_bounds(x_hi, 22)[1]
    main_lo = pi_lo**2 * angle_lo / (q*q*d_hi)
    main_hi = pi_hi**2 * angle_hi / (q*q*d_lo)
    left = '154.31470577028288820'
    right = '154.31470577028288822'
    assert F(left) < main_lo*(1-tail)**4
    assert main_hi*(1+tail)**4 < F(right)
    print(f'{left} < E(1/10^7) < {right}')


if __name__ == '__main__':
    verify()
