"""Exact arithmetic laboratories; standard library only. See THEORY.md.

Run from any directory: python path/to/verify.py
The output is deterministic and is written only alongside this file.
"""
from fractions import Fraction as F
from math import factorial, gcd, isqrt
from pathlib import Path
import json


def floor(x):
    return x.numerator // x.denominator


def ceil(x):
    return -floor(-x)


def ceil_div(a, b):
    return -((-a) // b)


def atan_bounds(k, terms):
    s = sum((F((-1) ** j, (2 * j + 1) * k ** (2 * j + 1))
             for j in range(terms)), F(0))
    next_s = s + F((-1) ** terms, (2 * terms + 1) * k ** (2 * terms + 1))
    return min(s, next_s), max(s, next_s)


def pi_bounds(terms=128):
    # Rational tangent addition verifies the Machin identity and its branch.
    tan2 = F(2, 5) / (1 - F(1, 25))
    # alpha=atan(1/5)<pi/4, hence 2alpha<pi/2; tan(2alpha)<1
    # gives 2alpha<pi/4. beta<alpha then puts 4alpha-beta in (0,pi/2).
    assert 0 < tan2 < 1 and 0 < F(1, 239) < F(1, 5)
    tan4 = 2 * tan2 / (1 - tan2 ** 2)
    assert tan4 == F(120, 119)
    assert (tan4 - F(1, 239)) / (1 + tan4 / 239) == 1
    l5, u5 = atan_bounds(5, terms)
    l239, u239 = atan_bounds(239, terms)
    return 16 * l5 - 4 * u239, 16 * u5 - 4 * l239


def dyadic_bounds(lo, hi, bits=256):
    den = 1 << bits
    lower, upper = floor(lo * den), ceil(hi * den)
    assert upper - lower == 1
    return lower, upper, den


def decimal_enclosure(lo, hi, digits=16):
    scale = 10 ** digits
    return [str(F(floor(lo * scale), scale)), str(F(ceil(hi * scale), scale))]


def abs_bounds(lo, hi):
    assert lo <= hi
    return (F(0) if lo <= 0 <= hi else min(abs(lo), abs(hi)), max(abs(lo), abs(hi)))


def enumerate_jet_box(lower, upper, den, p=355, q=113, C=20):
    """Enumerate every possible coefficient, without root/square pruning.

    D=theta=1/4; a>=0 suffices by central symmetry, and a=0 is
    included independently. Every interval is an outward enclosure.
    """
    assert den % 4 == 0
    den2, qq = den * den, q * q
    value_den = 4 * qq * den2
    found, candidates, boundary = [], 0, 0
    for a in range(C * qq + 1):
        bmin = ceil_div(-den - 8 * a * upper, 4 * den)
        bmax = (den - 8 * a * lower) // (4 * den)
        for b in range(bmin, bmax + 1):
            candidates += 1
            base_lo = a * lower * lower + b * (lower if b >= 0 else upper) * den
            base_hi = a * upper * upper + b * (upper if b >= 0 else lower) * den
            cmin = ceil_div(-den2 - 4 * qq * base_hi, value_den)
            cmax = (den2 - 4 * qq * base_lo) // value_den
            for c in range(cmin, cmax + 1):
                if a == b == c == 0:
                    continue
                deriv_lo, deriv_hi = 2 * a * lower + b * den, 2 * a * upper + b * den
                val_lo, val_hi = base_lo + c * den2, base_hi + c * den2
                inside = (4 * max(abs(deriv_lo), abs(deriv_hi)) <= den
                          and 4 * qq * max(abs(val_lo), abs(val_hi)) <= den2)
                if inside:
                    found.append([a, b, c])
                else:
                    boundary += 1
    assert boundary == 0, "Increase interval precision to resolve a boundary"
    # This comparison occurs AFTER the entire coefficient enumeration.
    expected = [[m * q * q, -2 * m * p * q, m * p * p] for m in range(1, C + 1)]
    assert found == expected
    return {"passed": True, "a_nonnegative_values": C * qq + 1,
            "derivative_candidates": candidates, "unresolved_boundaries": boundary,
            "positive_a_polynomials": found, "nonzero_full_box_count": 2 * len(found),
            "primitive_polynomial": [q*q, -2*p*q, p*p]}


def repair_jet_nonzero(lower, upper, den, escape_lower):
    """Try the costed repair M=1 in a unimodular square-product basis.

    This is a bounded search for a certificate, not a claim of global
    minimization. The proven escape lower bound applies globally.
    """
    p, q, u, v = 355, 113, 7, 22
    assert p*u-q*v == -1
    endpoints = [F(lower, den), F(upper, den)]
    assert all(q*q*x-p*q < 0 for x in endpoints)
    den2 = den*den
    tested, feasible, best = 0, 0, None
    for beta in range(-1000, 1):
        # h=alpha*(qX-p)^2+beta*(qX-p)*(uX-v)+(uX-v)^2
        ends = []
        for x in endpoints:
            s, t = q*q*x-p*q, q*(u*x-v)
            base_d = beta*(t+u*s/q)+2*u*t/q
            ends.extend([(F(-1,4)-base_d)/(2*s), (F(1,4)-base_d)/(2*s)])
        # Each endpoint function is fractional linear in x with no pole
        # in the pi interval; its extrema are attained at the endpoints.
        for alpha in range(ceil(min(ends)), floor(max(ends))+1):
            tested += 1
            a = alpha*q*q+beta*q*u+u*u
            b = -2*alpha*q*p-beta*(q*v+p*u)-2*u*v
            c = alpha*p*p+beta*p*v+v*v
            if a <= 0:
                continue
            dl, dh = 2*a*lower+b*den, 2*a*upper+b*den
            vl = a*lower*lower+b*(lower if b >= 0 else upper)*den+c*den2
            vh = a*upper*upper+b*(upper if b >= 0 else lower)*den+c*den2
            if 4*max(abs(dl), abs(dh)) <= den and 4*q*q*max(abs(vl), abs(vh)) <= den2:
                feasible += 1
                assert a*p*p+b*p*q+c*q*q == 1
                assert gcd(gcd(a,b),c) == 1
                disc = b*b-4*a*c
                if best is None or a < best[0]:
                    best = (a,b,c,disc,alpha,beta)
    assert best is not None
    a,b,c,disc,alpha,beta = best
    assert disc < 0  # irreducible over Q, no rational root
    assert escape_lower < F(a,q*q) < F(101,100)*escape_lower
    return {"passed": True, "search_basis": {"z": [q,-p], "w": [u,-v]},
            "beta_range": [-1000,0], "tested": tested, "certified_feasible": feasible,
            "candidate_coefficients": [a,b,c], "basis_coefficients": [alpha,beta,1],
            "discriminant": disc, "M_at_355_over_113": 1,
            "C_required_by_candidate": str(F(a,q*q)),
            "C_decimal": decimal_enclosure(F(a,q*q),F(a,q*q),10),
            "within_factor_of_global_necessary_cost": "101/100",
            "global_optimality_claimed": False}


def solve(matrix, rhs):
    n = len(rhs)
    a = [list(map(F, row)) + [F(b)] for row, b in zip(matrix, rhs)]
    for k in range(n):
        pivot = next(i for i in range(k, n) if a[i][k])
        a[k], a[pivot] = a[pivot], a[k]
        d = a[k][k]
        a[k] = [v / d for v in a[k]]
        for i in range(n):
            if i != k:
                d = a[i][k]
                a[i] = [u - d * v for u, v in zip(a[i], a[k])]
    return [row[-1] for row in a]


def det(matrix):
    a = [list(map(F, row)) for row in matrix]
    result = F(1)
    for k in range(len(a)):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return F(0)
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            result = -result
        d = a[k][k]
        result *= d
        for i in range(k + 1, len(a)):
            t = a[i][k] / d
            for j in range(k + 1, len(a)):
                a[i][j] -= t * a[k][j]
    return result


def recurrences(nmax):
    A, B = [1, 1], [0, 1]
    for n in range(1, nmax):
        A.append((2*n+1)*A[-1] + n*n*A[-2])
        B.append((2*n+1)*B[-1] + n*n*B[-2])
    return A, B


def moments(z, n):
    values = [F(z)]
    for k in range(n):
        values.append(F(1, 2*k+1) - values[-1])
    return values


def moment_roots(A, B):
    rows = []
    for m in range(11):
        roots = []
        for shift in (0, 1):
            # Congruence basis [1,(1+t),t(1+t),...] yields Schur root.
            v = [F(1, 2*i+2*shift+1) for i in range(m)]
            H = [[F(1, 2*i+2*j+2*shift+1) + F(1, 2*i+2*j+2*shift+3)
                  for j in range(m)] for i in range(m)]
            correction = sum((x*y for x, y in zip(v, solve(H, v))), F(0)) if m else F(0)
            root = correction if shift == 0 else 1 - correction
            n = 2*m + shift
            assert root == F(B[n], A[n])
            roots.append(str(4*root))
            # Independent raw Hankel determinants: affine identity at z=0,1,2.
            if m <= 6:
                ds = []
                for z in (0, 1, 2):
                    mu = moments(z, 2*m+shift)
                    ds.append(det([[mu[i+j+shift] for j in range(m+1)] for i in range(m+1)]))
                assert ds[2] - ds[1] == ds[1] - ds[0] != 0
                assert -ds[0] / (ds[1] - ds[0]) == root
        rows.append({"m": m, "pi_lower_root": roots[0], "pi_upper_root": roots[1]})
    return {"passed": True, "schur_m_range": [0, 10], "raw_hankel_m_range": [0, 6], "rows": rows}


def coefficient_costs(A, B, lo, hi):
    rows = []
    for n in range(1, 65):
        g = gcd(A[n], 4*B[n])
        q, p = A[n] // g, 4*B[n] // g
        elo, ehi = abs_bounds(lo - F(p, q), hi - F(p, q))
        assert elo > 0
        # Independent standard Legendre integral enclosure; strict bounds.
        low = F(2*factorial(n)**2, (2*n+1)*A[n]**2)
        high = 2*low
        assert low < elo < ehi < high
        # q interval for K=100, theta=1/2; use upper error (safe window).
        # |L|=q_n*|pi-p_n/q_n|. Empty intervals are recorded, not omitted.
        start = max(1, ceil(F(q, 100)))
        stop = floor(1 / (2*q*ehi))
        rows.append({"n": n, "raw_A_bits": A[n].bit_length(),
                     "raw_4B_bits": (4*B[n]).bit_length(), "gcd": str(g),
                     "p": str(p), "q": str(q), "q_bits": q.bit_length(),
                     "quality": decimal_enclosure(q*q*elo, q*q*ehi),
                     "K100_theta_half_safe_window": [str(start), str(stop)],
                     "window_nonempty": start <= stop})
    return {"passed": True, "n_range": [1, 64], "K": 100, "theta": "1/2", "rows": rows}


def quadratic_sources(lo, hi):
    u, v = [1, 3], [F(0), F(5)]
    for n in range(1, 25):
        factor = 11*n*n + 11*n + 3
        un = F(factor*u[-1] + n*n*u[-2], (n+1)**2)
        assert un.denominator == 1
        u.append(un.numerator)
        v.append((factor*v[-1] + n*n*v[-2]) / (n+1)**2)
    rows = []
    for n in range(1, 25):
        r = 6*v[n]/u[n]
        delta = F(30, (n+1)**2*u[n]*u[n+1])
        rn = 6*v[n+1]/u[n+1]
        assert abs(r-rn) == delta
        assert min(r, rn) < lo*lo < hi*hi < max(r, rn)
        den = r.denominator
        square = isqrt(r.numerator)**2 == r.numerator and isqrt(den)**2 == den
        assert not square, "All-p lower certificate needs an independent square exception"
        budget = 1/(2*den*delta)
        horizon = isqrt(budget.numerator // budget.denominator)
        assert 2*den*horizon*horizon*delta <= 1
        assert 2*den*(horizon+1)**2*delta > 1
        rows.append({"n": n, "A": str(r.numerator), "B": str(den),
                     "source_is_rational_square": square, "delta": str(delta),
                     "all_p_safe_denominator_horizon": str(horizon),
                     "quality_lower_constant": str(F(1, 16*den))})
    return {"passed": True, "source_range": [1, 24], "rows": rows}


def main():
    lo, hi = pi_bounds()
    assert F(3) < lo < hi < F(4)
    l, u, den = dyadic_bounds(lo, hi)
    p, q = 355, 113
    kl, ku = abs_bounds(q*q*lo-p*q, q*q*hi-p*q)
    assert F(3406311930, 10**12) < kl < ku < F(3406311931, 10**12)
    assert F(1,4) + ku/4 + 20*ku*ku < 1
    assert F(1,4) + 40*ku < 1
    assert 40*ku < F(1,4) and 20*ku*ku < F(1,4)
    # f(k)=(3/4-k/4)/k^2 is strictly decreasing on the tiny interval.
    escape_lo = (F(3,4)-ku/4)/(ku*ku)
    escape_hi = (F(3,4)-kl/4)/(kl*kl)
    assert escape_lo > 64564
    A, B = recurrences(64)
    result = {
        "status": "PASS; no assertion of pi badly approximable or its negation",
        "arithmetic": "Python integers/Fraction; alternating-tail bounds; no floating point decisions",
        "pi_interval": {"machin_terms_each": 128, "lo": str(lo), "hi": str(hi),
                        "dyadic_bits": 256, "dyadic_lower_numerator": str(l),
                        "dyadic_upper_numerator": str(u)},
        "jet_theorem_instance": {"p": p, "q": q, "C": 20, "D": "1/4", "theta": "1/4",
                                 "volume": 10, "kappa": decimal_enclosure(kl, ku, 20),
                                 "necessary_C_for_nonzero_M": decimal_enclosure(escape_lo, escape_hi, 10)},
        "jet_exhaustion": enumerate_jet_box(l, u, den),
        "jet_costed_repair": repair_jet_nonzero(l, u, den, escape_lo),
        "moments": moment_roots(A, B),
        "linear_costs": coefficient_costs(A, B, lo, hi),
        "quadratic_sources": quadratic_sources(lo, hi),
    }
    path = Path(__file__).with_name("verification.json")
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print("PASS: complete quadratic jet box, 40 nonzero polynomials, one primitive square direction")
    print("PASS: costed irreducible repair M=1; coefficient budget within 1% of necessary bound")
    print("PASS: 22 exact Schur roots and 14 raw Hankel determinants; same Lambert family")
    print("PASS: primitive coefficient/error/window costs n=1..64; pi^2 sources n=1..24")
    print("Saved verification.json; no infinite pi lower-bound conclusion")


if __name__ == "__main__":
    main()
