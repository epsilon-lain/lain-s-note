"""Bound E and the log residual at convergent indices 307 and 431 of 1/pi."""
from fractions import Fraction as F
from verify_theta_energy import atan_bounds, round_out


def log_bounds(v, terms=40):
    """Positive atanh expansion of log(v), with geometric tail (1 <= v)."""
    u = (v-1)/(v+1)
    assert 0 <= u < 1
    value = 2*sum((u**(2*k+1)/(2*k+1) for k in range(terms)), F(0))
    tail = 2*u**(2*terms+1)/((2*terms+1)*(1-u*u))
    return value, value+tail


def main():
    # At these indices q has at most 217 digits; 480 places resolve q^2 errors.
    a, b = atan_bounds(F(1,5), 180)
    c, d = atan_bounds(F(1,239), 52)
    pi_lo, pi_hi = round_out(16*a-4*d, 16*b-4*c, 480)
    lower, upper = 1/pi_hi, 1/pi_lo
    p0,p,q0,q = 0,1,1,0
    angle_lo,angle_hi = atan_bounds(F(1,3),20)
    # For these specific classes |k0|^4 is 1/4 (q=8 mod16) or 1 (odd q).
    targets = {307: (8, F(1,4), 4, F('55.2055'), F('55.2057')),
               431: (15, F(1), 1, F('10500.70'), F('10500.71'))}
    for j in range(max(targets)+1):
        digit = lower.__floor__()
        assert upper.__floor__() == digit
        p0,p = p,digit*p+p0
        q0,q = q,digit*q+q0
        lower,upper = 1/(upper-digit),1/(lower-digit)
        if j not in targets:
            continue
        res, k0four, L, left, right = targets[j]
        assert q%16 == res
        off_lo, off_hi = 2-2*pi_hi*p/q, 2-2*pi_lo*p/q
        assert off_lo * off_hi > 0
        tl,tu = sorted((abs(off_lo),abs(off_hi)))
        hl,hu = pi_lo**2/(10*q*q*tu), pi_hi**2/(10*q*q*tl)
        # hL^2>100. For all s in [1,2], x=exp(-Re Z_s * L^2)<10^-30:
        # e>2 and 2^100>10^30. All active indices are multiples of L
        # and |k_m/k0|<=3. Relative absolute amplitude error <= 6x/(1-x^3).
        assert L**2*hl > 100 and 2**100 > 10**30
        amplitude_error = F(7, 10**30)
        ml = 10*angle_lo*k0four*hl
        mu = 10*angle_hi*k0four*hu
        assert left < ml*(1-amplitude_error)**4
        assert mu*(1+amplitude_error)**4 < right
        print(f'j={j}, q mod 16={q%16}, {float(left)} < E(T*) < {float(right)}')

        # rho = exp(-L^2 h), y=rho^(1/4), rho<10^-30 and y<10^-7.
        # K=2*k_L/k_0 has |K|=1 for odd q and |K|=6 for q=8 mod 16.
        # At s=2 the normalized first mode is g=K exp(-L^2 Z_2).
        # With F_s/k_0=1+r_s, |r_2-g|<7*rho^4, |r_1|<7*rho*y.
        # d=(r_2-r_1)/(1+r_1) has |d-g|<9*rho*y and |d|<7*rho.
        # z=z0(1+d)^2, so |z-z0-2*z0*g|<19*|z0|*rho*y and
        # |z-z0|<15*|z0|*rho. Using P'(z0) of norm 2 gives
        # ||P(z)| - 4*|z0|*|K|*rho| < 39*|z0|*rho*y.
        # The relative error is <10*y<10^-6, hence log error <1/1000.
        l2,u2 = log_bounds(F(2))
        l54,u54 = log_bounds(F(5,4))
        l10,u10 = 3*l2+l54,3*u2+u54
        # log(4*sqrt(2/5)) = 3*log(2)-log(10)/2.
        lp,up = 3*l2-u10/2,3*u2-l10/2
        if res == 8:
            l3,u3 = log_bounds(F(3))
            lp,up = lp+l2+l3,up+u2+u3
        ln_lo,ln_hi = lp-L**2*hu-F(1,1000),up-L**2*hl+F(1,1000)
        assert ln_hi < 0
        log10_lo,log10_hi = ln_lo/l10,ln_hi/u10
        bounds = {307: ('-475.721','-475.718'),431: ('-1416.968','-1416.966')}
        low,high = bounds[j]
        assert F(low) < log10_lo <= log10_hi < F(high)
        print(f'{low} < log10(abs(R(T*))) < {high}')


if __name__ == '__main__':
    main()
