"""Certify an exact zero by sign changes on the four faces of a box.

The phase varies with K. This is a finite control close to x=1/pi.
Requires mpmath, run python pika_exact_residual_zero.py.
"""
from mpmath import iv

iv.dps = 65
p, q = 113, 355
D = 2 / (2 - iv.j)**2
inv_4p = pow(4*p, -1, q)


def scaled_modes(K, s):
    """Return exp(K) times the normalized nonconstant dual sum at s*T."""
    Z = 5*K / (2*(s-iv.j))
    total = iv.mpc(0)
    for m in (1, 2):
        j = (m + (q if m % 2 else 0)) // 2
        phase = -2*iv.pi*((inv_4p*j**2) % q)/q
        total += (1 if m % 2 else 2)*iv.exp(K - m*m*Z + phase*iv.j)
    U = Z.real
    tail = 2*iv.exp(K-9*U)/(1-iv.exp(-7*U))
    error = iv.mpf([-1, 1])*tail
    return total + iv.mpc(error, error)


def scaled_residual(K, v):
    small = iv.exp(-K)
    s = 1 + v*small
    u1 = scaled_modes(K, s)
    u2 = scaled_modes(K, 2*s)
    z0 = (s-iv.j)/(2*s-iv.j)
    den = 1 + small*u1
    assert abs(den).a > 0
    y = z0*(u2-u1)*(2+small*(u1+u2))/den**2
    return v*((s+1)/(D*(2*s-iv.j)**2)) + ((10*z0-6)*y + 5*small*y*y)/D


def parts(a, b, count):
    for m in range(count):
        start, end = a+(b-a)*m/count, a+(b-a)*(m+1)/count
        yield iv.mpf([start.a, end.b])


def window_tail_bound(K, lower_time, upper_time):
    # t/(1+t^2) attains its minimum on an interval at an endpoint.
    values = [5*K*t/(2*(1+t*t)) for t in (lower_time,upper_time)]
    decay = min(x.a for x in values)
    return 2*iv.exp(-decay)/(1-iv.exp(-3*decay))


def exclude_outer_scales(K, lower_s, upper_s, min_distance):
    # Everywhere |z0| <= 1 and |Pprime(z0)| <= 4, since
    # Pprime(z0) = (-2s - 4i)/(2s - i). Bound both Fourier tails.
    tail_one = window_tail_bound(K,lower_s,upper_s)
    tail_two = window_tail_bound(K,2*lower_s,2*upper_s)
    assert tail_one.b < 1
    ratio_error = (tail_one+tail_two)/(1-tail_one)
    z_error = ratio_error*(2+ratio_error)
    polynomial_error = 4*z_error+5*z_error**2
    # P(z0)=(s-1)(s+1)/(2s-i)^2 exactly.
    main_min = min_distance*(lower_s+1)/(4*upper_s**2+1)
    assert main_min.a-polynomial_error.b > iv.mpf("4e-22")


def fixed_phase_minimum():
    """Bound the minimum for the original phase over all |v| <= 50."""
    K = iv.pi/(20*q*(q/iv.pi-p))
    small = iv.exp(-K)
    all_v = iv.mpf([-50,50])
    lo = abs(D)*small*scaled_residual(K,all_v).imag
    assert lo.a > iv.mpf("3.99716660877867e-22")
    v_trial = -scaled_residual(K,iv.mpf(0)).real
    assert -50 < v_trial < 50
    hi = abs(D)*small*abs(scaled_residual(K,v_trial))
    assert hi.b < iv.mpf("3.99716660877871e-22")
    mid = iv.mpf("1e-5")
    exclude_outer_scales(K,1-mid,1+mid,50*small)
    exclude_outer_scales(K,iv.mpf("0.5"),iv.mpf(2),mid)
    print("3.99716660877867e-22 < min_{1/2<=s<=2} |R_1/pi(s*T*)| < 3.99716660877871e-22")


def main():
    fixed_phase_minimum()
    kl, kh = iv.mpf('46.143'), iv.mpf('46.154')
    vl, vh = iv.mpf('-6.40'), iv.mpf('-6.25')
    Kbox, vbox = iv.mpf([kl.a, kh.b]), iv.mpf([vl.a, vh.b])
    # Also check there is no pole over the full rectangle.
    sbox = 1 + vbox*iv.exp(-Kbox)
    assert abs(1 + iv.exp(-Kbox)*scaled_modes(Kbox,sbox)).a > 0
    for cell in parts(kl,kh,16):
        assert scaled_residual(cell,vl).real.b < 0
        assert scaled_residual(cell,vh).real.a > 0
    for cell in parts(vl,vh,16):
        assert scaled_residual(kl,cell).imag.a > 0
        assert scaled_residual(kh,cell).imag.b < 0
    x = iv.mpf(p)/q + iv.pi/(20*q*q*Kbox)
    gap = x-1/iv.pi
    assert iv.mpf('-2.4e-11') < gap < iv.mpf('-1.6e-11')
    print("Opposite signs on each pair of faces; |W|>0 throughout.")
    print("-2.4e-11 < x-1/pi < -1.6e-11 for every point of the box.")


if __name__ == '__main__':
    main()
