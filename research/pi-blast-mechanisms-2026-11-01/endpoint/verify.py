"""Exact checks of sample parameters and displayed constants.

The universal obstruction is proved in CONSTANT_GAIN_OBSTRUCTION.md.
Finite checks are not a substitute for that proof. Standard library only.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from pathlib import Path
import json


def log_unit_bounds(x, terms=90):
    assert 1 <= x <= 2
    z = (x - 1) / (x + 1)
    lo = 2 * sum((z**(2*k+1) / (2*k+1) for k in range(terms)), F())
    tail = 2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
    return lo, lo+tail


def log_bounds(x):
    assert x > 0
    if x < 1:
        lo, hi = log_bounds(1/x)
        return -hi, -lo
    k = 0
    while x >= 2:
        x /= 2
        k += 1
    lo, hi = log_unit_bounds(x)
    l2, u2 = log_unit_bounds(F(2))
    return lo+k*l2, hi+k*u2


def decimal(x, rounding):
    with localcontext() as ctx:
        ctx.prec = 18
        ctx.rounding = rounding
        return str(Decimal(x.numerator)/Decimal(x.denominator))


def main():
    l2, u2 = log_bounds(F(2))
    clo, chi = l2/4, u2/4
    cases = []
    for eps in (F(1), F(1, 2), F(1, 10), F(1, 100)):
        nu = 2+eps
        b = (F(1, 2) + 1 - 1/nu)/2
        delta = min(F(1, 4), (2*b-1)/(2*b*b))
        theta, A = 1-delta, 1-b*delta
        C = (A/theta+1/A)/2
        B = (A+min(1/C, C*theta))/2
        g0 = nu*A-1-(nu-1)*theta
        eta = g0/(2*nu*A)
        g = g0-nu*A*eta
        assert 0 < theta < A < B < 1
        assert C > 1 and B < 1/C and B < C*theta < 1
        assert 0 < g < g0 < eps*eps/(4*(1+eps))
        assert 0 < eta < eps*eps/(4*(1+eps))
        rlo, rhi = log_bounds(B/A)
        hlo, hhi = log_bounds(C)
        assert rlo > 0 and rhi < eps*hlo/2
        # Bounds for the necessary w_* threshold, with ln(2) also enclosed.
        arglo = 64*(1+eps)**2/(chi*eps**4)
        arghi = 64*(1+eps)**2/(clo*eps**4)
        lower = 16*(1+eps)/eps**3 * log_bounds(arglo)[0]
        upper = 16*(1+eps)/eps**3 * log_bounds(arghi)[1]
        assert lower > 2/eps**3
        cases.append({
            'epsilon': str(eps),
            'theta': str(theta), 'A': str(A), 'B': str(B), 'C': str(C),
            'eta': str(eta),
            'strict_r_less_than_epsilon_h_over_2': True,
            'required_w_star_threshold_interval_display': [decimal(lower, ROUND_FLOOR), decimal(upper, ROUND_CEILING)],
            'threshold_greater_than_2_epsilon_power_minus_3': True,
        })
    report = {
        'status': 'passed', 'arithmetic': 'exact rational intervals',
        'cases': cases,
        'scope': 'Checks admissible sample parameter families and constants; the universal theorem is the analytic proof in the accompanying note.',
        'pi_badly_approximable_decided': False,
    }
    Path(__file__).with_name('verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
