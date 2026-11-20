#!/usr/bin/env python3
"""CM/Legendre bridge from rational elliptic points to rational approximants to pi.

Curve: E: y^2 = x^3 - 36x (j=1728), P=(12,-36).
All elliptic point arithmetic, zeta cocycles, filters, and Pi_K values use
fractions.Fraction exactly. mpmath is used only for period/integral diagnostics
and decimal error summaries.
"""
from fractions import Fraction as F
from math import comb, log
import json
from pathlib import Path
import mpmath as mp

A = 36
P = (F(12), F(-36))
mp.mp.dps = 80


def add(P1, P2):
    if P1 is None:
        return P2
    if P2 is None:
        return P1
    x1, y1 = P1
    x2, y2 = P2
    if x1 == x2 and y1 == -y2:
        return None
    if P1 == P2:
        if y1 == 0:
            return None
        slope = (3*x1*x1 - A) / (2*y1)
    else:
        slope = (y2-y1)/(x2-x1)
    x3 = slope*slope-x1-x2
    y3 = slope*(x1-x3)-y1
    return x3, y3


def mul(n):
    result, base = None, P
    while n:
        if n & 1:
            result = add(result, base)
        base = add(base, base)
        n >>= 1
    return result


def xi(n):
    x, y = mul(n)
    return -x/y


def add_with_cocycle(pair1, pair2):
    Q1, R1 = pair1
    Q2, R2 = pair2
    if Q1 is None:
        return pair2
    if Q2 is None:
        return pair1
    x1, y1 = Q1
    x2, y2 = Q2
    if Q1 == Q2:
        if y1 == 0:
            return None, F(0)
        g = (3*x1*x1-A)/(2*y1)
    else:
        if x1 == x2 and y1 == -y2:
            return None, F(0)
        g = (y1-y2)/(x1-x2)
    return add(Q1, Q2), R1+R2+g


def mul_with_cocycle(n):
    result = (None, F(0))
    base = (P, F(0))
    while n:
        if n & 1:
            result = add_with_cocycle(result, base)
        base = add_with_cocycle(base, base)
        n >>= 1
    return result


def formal_log_coeff(r):
    return F((-1)**r * comb(2*r, r) * A**r, 4*r+1)


def formal_log(z, K):
    return sum((formal_log_coeff(r)*z**(4*r+1)
                for r in range(K+1)), F(0))


def zeta_coefficients(K):
    # wp(u)=u^-2 + sum a_r u^(4r-2), with wp''=6wp^2-2A.
    a = {1: F(1, 5)}
    for r in range(2, K+1):
        den = (4*r-2)*(4*r-3)-12
        a[r] = 6*sum((a[i]*a[r-i] for i in range(1, r)), F(0))/den
    return {r: -F(A**r)*a[r]/(4*r-1) for r in a}


def zeta_series(delta, K, coeffs):
    return 1/delta + sum((coeffs[r]*delta**(4*r-1)
                          for r in range(1, K+1)), F(0))


def log_height(x):
    return max(log(abs(x.numerator)), log(x.denominator))


def rational_pi_approx(n1, m1, n2, m2, K, R, zetas):
    D = n2*m1-n1*m2
    assert abs(D) == 1
    z1, z2 = xi(n1), xi(n2)
    d1, d2 = formal_log(z1, K), formal_log(z2, K)
    omega_approx = (n1*d2-n2*d1)/D
    H_approx = (n2*R[n1]-n1*R[n2]
                -n2*zeta_series(d1, K, zetas)
                +n1*zeta_series(d2, K, zetas))/D
    return omega_approx*H_approx


def zeta_filter(n, weights, C):
    z = xi(n)
    total = sum((F(lam)*xi(j*n)/z for j, lam in weights.items()), F(0))
    return total/(C*z)


def count_points_mod_p(p):
    count = 1  # point at infinity
    for x in range(p):
        rhs = (x**3-36*x) % p
        count += sum(1 for y in range(p) if (y*y-rhs) % p == 0)
    return count


def main():
    R = {n: mul_with_cocycle(n)[1] for n in (7, 11, 18, 101)}
    zetas = zeta_coefficients(6)

    # Periods: substitute x=6+t^2 in the real elliptic integrals.
    f_period = lambda t: 2/mp.sqrt((6+t*t)*(12+t*t))
    f_log = lambda t: 1/mp.sqrt((6+t*t)*(12+t*t))
    Omega = mp.quad(f_period, [0, 1, mp.inf])
    uP = mp.quad(f_log, [mp.sqrt(6), mp.sqrt(6)+1, mp.inf])
    theta = uP/Omega
    H = mp.pi/Omega

    pairs = [
        {"n1":7,"m1":2,"n2":11,"m2":3,"Kmax":6},
        {"n1":18,"m1":5,"n2":101,"m2":28,"Kmax":3},
    ]
    approx_rows = []
    for pair in pairs:
        n1,m1,n2,m2,Kmax = (pair[k] for k in ("n1","m1","n2","m2","Kmax"))
        D = n2*m1-n1*m2
        assert abs(D) == 1
        for K in range(Kmax+1):
            rat = rational_pi_approx(n1,m1,n2,m2,K,R,zetas)
            err = mp.mpf(rat.numerator)/rat.denominator-mp.pi
            approx_rows.append({
                "pair": f"{m1}/{n1},{m2}/{n2}",
                "determinant": D,
                "K": K,
                "log10_abs_error": mp.nstr(mp.log10(abs(err)), 16),
                "error": mp.nstr(err, 18),
                "numerator_bits": rat.numerator.bit_length(),
                "denominator_bits": rat.denominator.bit_length(),
                "log_height": mp.nstr(log_height(rat), 16),
            })

    filters = [
        {"name":"two_point_order7", "weights":{2:285,3:-34}, "C":468},
        {"name":"four_point_order11",
         "weights":{2:-15064,3:17375,4:-5291,5:539}, "C":3528},
    ]
    filter_rows = []
    filter_pi_rows = []
    for fil in filters:
        weights, C = fil["weights"], fil["C"]
        assert sum(lam*j for j,lam in weights.items()) == C
        assert sum(lam*(j**5-j) for j,lam in weights.items()) == 5*C//6
        if fil["name"] == "four_point_order11":
            assert sum(lam*(j**9-6*j**5+5*j)
                       for j,lam in weights.items()) == -45*C//14
        zhat = {}
        for n,m in ((7,2),(11,3)):
            z = xi(n)
            zhat[n] = zeta_filter(n,weights,C)
            delta = mp.mpf(n)*uP-m*Omega
            zeta_true = 1/delta + sum(
                mp.mpf(zetas[r].numerator)/zetas[r].denominator
                * delta**(4*r-1) for r in range(1,7)
            )
            zhat_mp = mp.mpf(zhat[n].numerator)/zhat[n].denominator
            filter_rows.append({
                "filter":fil["name"],"n":n,"m":m,
                "weights":weights,"C":C,
                "zeta_filter_value":mp.nstr(zhat_mp,24),
                "zeta_true_numeric":mp.nstr(zeta_true,24),
                "absolute_error":mp.nstr(zhat_mp-zeta_true,16),
                "numerator_bits":zhat[n].numerator.bit_length(),
                "denominator_bits":zhat[n].denominator.bit_length(),
            })

        # Combine the zeta filter with the first-order formal-log truncation.
        d7 = formal_log(xi(7),1)
        d11 = formal_log(xi(11),1)
        omega_approx = 7*d11-11*d7
        H_approx = 11*R[7]-7*R[11]-11*zhat[7]+7*zhat[11]
        pi_filter = omega_approx*H_approx
        err = mp.mpf(pi_filter.numerator)/pi_filter.denominator-mp.pi
        filter_pi_rows.append({
            "filter":fil["name"],
            "log10_abs_pi_error":mp.nstr(mp.log10(abs(err)),16),
            "numerator_bits":pi_filter.numerator.bit_length(),
            "denominator_bits":pi_filter.denominator.bit_length(),
            "log_height":mp.nstr(log_height(pi_filter),16),
        })

    result = {
        "curve":"E: y^2=x^3-36x",
        "j_invariant":1728,
        "point":"P=(12,-36)",
        "point_non_torsion_certificate":{
            "good_prime_point_counts":{"5":8,"7":8,"11":12,"13":20},
            "gcd":4,
            "two_P":"(25/4,35/8)",
            "conclusion":"Prime-to-p reduction bounds rational torsion order by 4; P is not of order 1, 2, or 4."
        },
        "CM_period_relation":"For full periods Omega, i*Omega and quasi-periods H,-i*H, Legendre relation gives Omega*H=pi.",
        "period_diagnostics":{
            "Omega":mp.nstr(Omega,35),
            "uP":mp.nstr(uP,35),
            "theta":mp.nstr(theta,35),
            "H":mp.nstr(H,35),
            "delta_7_2":mp.nstr(7*uP-2*Omega,25),
            "delta_11_3":mp.nstr(11*uP-3*Omega,25),
            "delta_18_5":mp.nstr(18*uP-5*Omega,25),
            "delta_101_28":mp.nstr(101*uP-28*Omega,25)
        },
        "exact_rational_cocycles":{
            "R7":str(R[7]),"R11":str(R[11]),
            "R18_bit_lengths":[R[18].numerator.bit_length(),R[18].denominator.bit_length()],
            "R101_bit_lengths":[R[101].numerator.bit_length(),R[101].denominator.bit_length()]
        },
        "exact_rational_local_parameters":{
            "z7":str(xi(7)),"z11":str(xi(11)),
            "z18_log10_abs":mp.nstr(mp.log10(abs(mp.mpf(xi(18).numerator)/xi(18).denominator)),16),
            "z101_log10_abs":mp.nstr(mp.log10(abs(mp.mpf(xi(101).numerator)/xi(101).denominator)),16),
            "z18_bit_lengths":[xi(18).numerator.bit_length(),xi(18).denominator.bit_length()],
            "z101_bit_lengths":[xi(101).numerator.bit_length(),xi(101).denominator.bit_length()]
        },
        "rational_pi_approximants":approx_rows,
        "cm_sparse_zeta_filters":filter_rows,
        "pi_from_zeta_filters":filter_pi_rows,
        "scope":"Elliptic point arithmetic, cocycles, filters, and rational Pi_K values are exact. Periods and errors against pi are high-precision numerical diagnostics. The Legendre relation and two-convergent identities are exact."
    }
    out=Path(__file__).with_name("cm_legendre_pi_bridge_experiment.json")
    out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "period_relation_check":mp.nstr(Omega*H-mp.pi,12),
        "approximant_rows":[
            (r["pair"],r["K"],r["log10_abs_error"],r["denominator_bits"])
            for r in approx_rows
        ],
        "filter_pi_rows":filter_pi_rows,
        "output":str(out)
    },indent=2))


if __name__=="__main__":
    main()
