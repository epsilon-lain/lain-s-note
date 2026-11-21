#!/usr/bin/env python3
"""Independent interval checks for the record block and its single pole.

Run with Python 3 and mpmath 1.3.x: python3 verify_record_echo.py
All comparisons use intervals. The JSON display values are rounded midpoints.
The finite tests use the 113/355 convergent of 1/pi and a quadratic control.
"""
import json
from mpmath import iv, mp

iv.dps = mp.dps = 75


def csc2(x):
    s = iv.sin(x)
    assert s.a > 0 or s.b < 0, "precision insufficient to enclose the pole away from zero"
    return 1 / s**2


def nearint(x):
    # Midpoint rounding only proposes an integer. The interval check certifies it.
    h = int(mp.floor(mp.mpf(x.mid) + mp.mpf("0.5")))
    assert abs(x - h) < iv.mpf(1)/2
    return h


def display(x):
    return mp.nstr(mp.mpf(x.mid), 14)


def audit(name, alpha, phase, a, r, s, A, expected_theta=None):
    one = iv.mpf(1)
    d = r*r*alpha - r*a
    e = abs(d)
    assert e > 0 and e < one/8 and r >= 5 and A >= 4
    theta = 1/e - A - one*s/r
    assert theta > 0 and theta < 1
    if expected_theta is not None:
        assert abs(theta - expected_theta) < iv.mpf('1e-65')
    sign = 1 if d > 0 else -1
    assert (a*s + sign) % r == 0
    inverse_a = pow(a, -1, r)
    baseline = (one - one/r**2)/3
    sample = {0, 1, 2, A//4, A-1, A, A+1, 2*A, 5*A+3,
              10**6+3, 10**12+7, 10**24+9, 10**30+17}
    max_error = mp.mpf(0)
    sums = {}
    partial = iv.mpf(0)
    crossover = {}

    def block(j):
        # Direct phase sum, independent of the grid proof and its proxy.
        return sum((csc2(phase*n) for n in range(j*r+1, (j+1)*r+1)),
                   iv.mpf(0)) / r**2

    def check_proxy(j, B):
        nonlocal max_error
        k = ((-nearint(j*d)*inverse_a-1) % r) + 1
        proxy = csc2(iv.pi*d*(j + one*k/r))
        residual = abs(B - proxy)
        assert residual < 40*e
        max_error = max(max_error, mp.mpf(residual.b))
        if j == A:
            assert k == s

    for j in range(A+1):
        B = block(j)
        partial += B
        check_proxy(j, B)
        if j in (0, A):
            sums[j] = B
        K = j+1
        if K in (A, A+1):
            ordinary = sum((one/h**2 for h in range(1, K+1)), iv.mpf(0))
            incoming = sum((1/(h+theta)**2 for h in range(A-K+1, A+1)), iv.mpf(0))
            expected = (ordinary + incoming) / (iv.pi*K*e)**2
            H = partial/K**2
            deficit = H-expected
            assert deficit > 0 and deficit < one*4/K
            crossover[str(K)] = display(H)
    for j in sorted(sample.difference(range(A+1))):
        check_proxy(j, block(j))

    assert abs(sums[0] - baseline - 1/(iv.pi*e)**2) < 5*e + one/r**2
    assert abs(sums[A] - baseline - 1/(iv.pi*e*theta)**2) < 10*e + one/r**2
    return {"case": name, "a": a, "r": r, "A": A,
            "blocks_checked": len(set(range(A+1)) | sample), "largest_j": max(sample),
            "e_approx": display(e), "theta_approx": display(theta),
            "initial_block_approx": display(sums[0]), "translated_block_approx": display(sums[A]),
            "H_at_Kr_approx": crossover,
            "largest_sampled_abs_proxy_error_approx": mp.nstr(max_error, 12)}


def main():
    theta = (iv.sqrt(5)-1)/2
    tail = 11+theta
    control = (5*tail+2)/(17*tail+7)
    cases = [
        audit("reciprocal pi", 1/iv.pi, iv.mpf(1), 113, 355, 333, 292),
        audit("quadratic [0;3,2,2,11,1,1,...]", control, iv.pi*control,
              5, 17, 7, 11, theta),
    ]
    print(json.dumps({"verified": True, "interval_decimal_precision": iv.dps,
                      "cases": cases}, indent=2))


if __name__ == '__main__':
    main()
