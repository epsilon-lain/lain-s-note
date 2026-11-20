#!/usr/bin/env python3
"""Exact Gaussian-CM trace filter experiment for E: y^2=x^3-36x.

All elliptic group arithmetic, Gaussian endomorphisms, traces, filters,
cocycles, and rational pi approximants use Fraction arithmetic. The only
floating-point values are decimal summaries of exact rational intervals.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import comb, log10
import json
from pathlib import Path


@dataclass(frozen=True)
class G:
    """Gaussian rational a+bi."""
    r: F = F(0)
    i: F = F(0)

    def __init__(self, r=0, i=0):
        object.__setattr__(self, "r", F(r))
        object.__setattr__(self, "i", F(i))

    def __add__(self, other):
        other = other if isinstance(other, G) else G(other)
        return G(self.r + other.r, self.i + other.i)

    __radd__ = __add__

    def __neg__(self):
        return G(-self.r, -self.i)

    def __sub__(self, other):
        return self + (-other if isinstance(other, G) else G(-other))

    def __rsub__(self, other):
        return G(other) - self

    def __mul__(self, other):
        other = other if isinstance(other, G) else G(other)
        return G(self.r * other.r - self.i * other.i,
                 self.r * other.i + self.i * other.r)

    __rmul__ = __mul__

    def inv(self):
        den = self.r * self.r + self.i * self.i
        return G(self.r / den, -self.i / den)

    def __truediv__(self, other):
        other = other if isinstance(other, G) else G(other)
        return self * other.inv()

    def __rtruediv__(self, other):
        return G(other) * self.inv()


E_A = 36
P_RATIONAL = (F(12), F(-36))


def point_add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and y1 == -y2:
        return None
    if P == Q:
        if y1 == G(0):
            return None
        slope = (G(3) * x1 * x1 - G(E_A)) / (G(2) * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope * slope - x1 - x2
    y3 = slope * (x1 - x3) - y1
    return x3, y3


def point_mul(n, P):
    result, base = None, P
    while n:
        if n & 1:
            result = point_add(result, base)
        base = point_add(base, base)
        n >>= 1
    return result


def rational_point(n):
    return point_mul(n, tuple(G(c) for c in P_RATIONAL))


def xi(P):
    x, y = P
    return -x / y


def tau(P, sign=1):
    x, y = P
    return -x, G(0, sign) * y


def gaussian_endomorphism(P, a, b, sign=1):
    """[a +/- bi]P = [a]P + [b]tau^{+/-1}(P)."""
    return point_add(point_mul(a, P), point_mul(b, tau(P, sign)))


def trace_multiplier(P, a, b):
    """T_(a,b) = (xi((a+bi)P)+xi((a-bi)P))/(2a xi(P))."""
    assert a != 0
    value = (xi(gaussian_endomorphism(P, a, b, +1))
             + xi(gaussian_endomorphism(P, a, b, -1)))
    value /= G(2 * a) * xi(P)
    assert value.i == 0
    return value.r


def rational_cocycles(max_n):
    """R_(n+1)=R_n+(y(nP)-y(P))/(x(nP)-x(P))."""
    P = P_RATIONAL
    x1, y1 = P
    R = {1: F(0), 2: (3*x1*x1 - E_A) / (2*y1)}
    for n in range(2, max_n):
        xn_g, yn_g = rational_point(n)
        assert xn_g.i == 0 and yn_g.i == 0
        xn, yn = xn_g.r, yn_g.r
        R[n + 1] = R[n] + (yn-y1)/(xn-x1)
    return R


def xi_rational(n):
    return xi(rational_point(n)).r


def r1(a, b):
    return a**4 - 10*a*a*b*b + 5*b**4 - 1


def r2(a, b):
    return (a**8 - 36*a**6*b**2 + 126*a**4*b**4
            - 84*a*a*b**6 + 9*b**8 - 6*a**4
            + 60*a*a*b*b - 30*b**4 + 5)


WEIGHTS = {
    (1, 0): F(233, 210),
    (1, 1): F(-13, 105),
    (2, 0): F(1, 70),
}
FIXED_INTEGER_WEIGHTS = {2: 1478890, 3: -241928, 4: 15301}
FIXED_INTEGER_SCALE = 2293200


def cm_filter(n):
    P = rational_point(n)
    z = xi(P).r
    F_value = sum((weight * trace_multiplier(P, a, b)
                   for (a, b), weight in WEIGHTS.items()), F(0))
    return F_value / z


def fixed_integer_filter(n):
    z = xi_rational(n)
    numerator = sum((F(weight) * xi_rational(j*n)
                     for j, weight in FIXED_INTEGER_WEIGHTS.items()), F(0))
    return numerator / (F(FIXED_INTEGER_SCALE) * z * z)


def local_zeta_partial(z, R):
    """G(z)=zeta(L(z))=1/z-sum c_r z^(4r-1)/(2(4r-1))."""
    total = F(1, z)
    for r in range(1, R + 1):
        cr = (-1)**r * comb(2*r, r) * 36**r
        total -= F(cr, 2*(4*r-1)) * z**(4*r-1)
    return total


def local_zeta_interval(z, R=30):
    partial = local_zeta_partial(z, R)
    r = R + 1
    cr = (-1)**r * comb(2*r, r) * 36**r
    next_term = -F(cr, 2*(4*r-1)) * z**(4*r-1)
    assert 144 * abs(z)**4 < 1
    return min(partial, partial + next_term), max(partial, partial + next_term)


def log10_fraction(x):
    x = abs(x)
    return log10(x.numerator) - log10(x.denominator)


def error_log_interval(value, true_interval):
    lo, hi = true_interval
    err_lo, err_hi = value - hi, value - lo
    if err_lo > 0:
        return [log10_fraction(err_lo), log10_fraction(err_hi)]
    if err_hi < 0:
        return [log10_fraction(err_hi), log10_fraction(err_lo)]
    raise AssertionError("error interval unexpectedly contains zero")


def L1(z):
    return z - F(72, 5) * z**5


def atan_interval(base, terms=100):
    partial = sum((F((-1)**r, (2*r+1)*base**(2*r+1))
                   for r in range(terms)), F(0))
    neighbor = partial + F((-1)**terms,
                           (2*terms+1)*base**(2*terms+1))
    return min(partial, neighbor), max(partial, neighbor)


def machin_pi_interval(terms=100):
    lo5, hi5 = atan_interval(5, terms)
    lo239, hi239 = atan_interval(239, terms)
    return 16*lo5 - 4*hi239, 16*hi5 - 4*lo239


def pi_error_log_interval(value, pi_interval):
    lo, hi = pi_interval
    if value < lo:
        return [log10_fraction(lo-value), log10_fraction(hi-value)], "below"
    if value > hi:
        return [log10_fraction(value-hi), log10_fraction(value-lo)], "above"
    raise AssertionError("Machin interval contains rational approximant")


def make_pi_approximant(zeta_filter, R):
    z7, z11 = xi_rational(7), xi_rational(11)
    cocycles = rational_cocycles(11)
    omega_hat = 7*L1(z11) - 11*L1(z7)
    H_hat = (11*cocycles[7] - 7*cocycles[11]
             - 11*zeta_filter(7) + 7*zeta_filter(11))
    return omega_hat * H_hat


def main():
    z7, z11 = xi_rational(7), xi_rational(11)
    assert 144*abs(z7)**4 < 1 and 144*abs(z11)**4 < 1

    # Verify the three exact moment conditions.
    assert sum(WEIGHTS.values()) == 1
    assert sum((w*r1(a, b) for (a, b), w in WEIGHTS.items()), F(0)) == F(5, 6)
    assert sum((w*r2(a, b) for (a, b), w in WEIGHTS.items()), F(0)) == F(-45, 14)

    # Verify the CM automorphism and rationality of all trace values.
    P7 = rational_point(7)
    assert tau(P7, +1) != tau(P7, -1)
    for a, b in WEIGHTS:
        assert isinstance(trace_multiplier(P7, a, b), F)

    pi_interval = machin_pi_interval(100)
    rows = []
    for n, z in ((7, z7), (11, z11)):
        cm_value = cm_filter(n)
        fixed_value = fixed_integer_filter(n)
        true_zeta = local_zeta_interval(z, 30)
        rows.append({
            "n": n,
            "xi_n_numerator_bits": abs(z.numerator).bit_length(),
            "xi_n_denominator_bits": z.denominator.bit_length(),
            "cm_filter_numerator_bits": abs(cm_value.numerator).bit_length(),
            "cm_filter_denominator_bits": cm_value.denominator.bit_length(),
            "cm_filter_log10_error_interval": error_log_interval(cm_value, true_zeta),
            "fixed_integer_filter_numerator_bits": abs(fixed_value.numerator).bit_length(),
            "fixed_integer_filter_denominator_bits": fixed_value.denominator.bit_length(),
            "fixed_integer_filter_log10_error_interval":
                error_log_interval(fixed_value, true_zeta),
        })

    pi_cm = make_pi_approximant(cm_filter, 30)
    pi_fixed = make_pi_approximant(fixed_integer_filter, 30)
    cm_err, cm_side = pi_error_log_interval(pi_cm, pi_interval)
    fixed_err, fixed_side = pi_error_log_interval(pi_fixed, pi_interval)

    R = rational_cocycles(11)
    result = {
        "curve": "E: y^2=x^3-36x",
        "point": "P=(12,-36)",
        "cm_ring": "Z[i]",
        "automorphism": "tau(x,y)=(-x,i*y)",
        "gaussian_trace_filter_weights": {
            "1": "233/210",
            "1+i": "-13/105",
            "2": "1/70"
        },
        "moment_checks": {
            "sum_weights": "1",
            "sum_weight_r1": "5/6",
            "sum_weight_r2": "-45/14"
        },
        "rational_cocycles": {"R7": str(R[7]), "R11": str(R[11])},
        "local_filter_results": rows,
        "pi_approximants": {
            "gaussian_cm": {
                "numerator_bits": abs(pi_cm.numerator).bit_length(),
                "denominator_bits": pi_cm.denominator.bit_length(),
                "log10_error_interval": cm_err,
                "side_of_pi": cm_side
            },
            "fixed_integer_three_point": {
                "numerator_bits": abs(pi_fixed.numerator).bit_length(),
                "denominator_bits": pi_fixed.denominator.bit_length(),
                "log10_error_interval": fixed_err,
                "side_of_pi": fixed_side
            }
        },
        "scope": (
            "All elliptic arithmetic, Gaussian endomorphism traces, filters, "
            "cocycles, and rational approximants are exact. Error intervals "
            "use the exact Machin arctangent enclosure. Local zeta error "
            "intervals use the alternating formal series. This is a finite "
            "Legendre-bridge experiment."
        )
    }
    out = Path(__file__).with_name("gaussian_cm_trace_filter_results.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
