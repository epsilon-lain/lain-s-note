"""Exact finite spectral observation, using only rational target intervals.

The search is a diagnostic, not a new pi approximation algorithm or theorem.
No target decimal, continued-fraction digits, or floating-point test is used.
"""
from fractions import Fraction as F
from math import gcd, isqrt
import json
from pathlib import Path


def floor(x):
    return x.numerator // x.denominator


def ceil(x):
    return -floor(-x)


def sqrt_floor(x):
    assert x >= 0
    return isqrt(x.numerator // x.denominator)


def atan_interval(inv, terms=90):
    value = sum((F((-1)**k, (2*k+1)*inv**(2*k+1))
                 for k in range(terms)), F(0))
    nxt = F((-1)**terms, (2*terms+1)*inv**(2*terms+1))
    return min(value, value+nxt), max(value, value+nxt)


def pi_interval():
    a, b = atan_interval(5)
    c, d = atan_interval(239)
    return 16*a-4*d, 16*b-4*c


def residual_range(p, q, target):
    a, b = target
    left, right = p-q*b, p-q*a
    lo = F(0) if left <= 0 <= right else min(abs(left), abs(right))
    return lo, max(abs(left), abs(right))


def finite_spectrum(target, scale, cutoff, primitive=True):
    """Enclose all q>0 modes with q²/scale+scale*(p-alpha*q)²<=cutoff.

    Each returned row is certain or unresolved. Omitted modes are certified
    above cutoff for every alpha in target. Interval inputs may be negative.
    """
    a, b = target
    assert a <= b and scale > 0 and cutoff > 0
    qmax = sqrt_floor(scale*cutoff)
    radius = sqrt_floor(cutoff/scale)+1
    rows = []
    for q in range(1, qmax+1):
        for p in range(ceil(q*a-radius), floor(q*b+radius)+1):
            if primitive and gcd(p, q) != 1:
                continue
            lo, hi = residual_range(p, q, target)
            emin, emax = F(q*q)/scale+scale*lo*lo, F(q*q)/scale+scale*hi*hi
            if emin <= cutoff:
                rows.append(dict(p=p, q=q, emin=emin, emax=emax,
                                 quality_lo=q*lo, quality_hi=q*hi,
                                 certain=emax <= cutoff))
    return qmax, rows


def encode(x):
    if isinstance(x, F):
        return dict(exact=str(x), display=float(x))
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def run():
    checks = 0
    # Exhaust every p in a larger integer box for rational targets; test the
    # completeness boundary and the primitive/harmonic distinction.
    for target in [F(-7, 3), F(0), F(1, 3), F(355, 113), F(7, 5)]:
        for scale in [F(1, 3), F(1), F(7), F(100)]:
            for cutoff in [F(1, 10), F(1), F(3)]:
                for primitive in [False, True]:
                    cap, rows = finite_spectrum((target, target), scale, cutoff, primitive)
                    observed = {(r['p'], r['q']) for r in rows}
                    expected = set()
                    for q in range(1, cap+5):
                        for p in range(-100, 101):
                            if primitive and gcd(p, q) != 1:
                                continue
                            energy = F(q*q)/scale+scale*(p-target*q)**2
                            if energy <= cutoff:
                                expected.add((p, q))
                    assert observed == expected
                    assert all(r['emin'] == r['emax'] and r['certain'] for r in rows)
                    checks += 1
    # Interval validity: all exact sample targets within interval obey bounds.
    box = (F(31, 10), F(32, 10))
    cap, rows = finite_spectrum(box, F(100), F(2))
    by_pair = {(r['p'], r['q']): r for r in rows}
    for j in range(21):
        target = box[0]+(box[1]-box[0])*F(j, 20)
        _, point_rows = finite_spectrum((target, target), F(100), F(2))
        for r in point_rows:
            parent = by_pair[r['p'], r['q']]
            assert parent['emin'] <= r['emin'] <= parent['emax']
            checks += 1
    target = pi_interval()
    scale, cutoff = F(4_000_000), F(1, 10)
    cap, rows = finite_spectrum(target, scale, cutoff)
    assert len(rows) == 1 and (rows[0]['p'], rows[0]['q']) == (355, 113)
    assert rows[0]['certain']
    assert rows[0]['quality_hi'] < rows[0]['emax']/2
    assert rows[0]['quality_hi'] < F(341, 100000)
    # Harmonics are modes of the same reduced fraction, not new discoveries.
    _, harmonic_rows = finite_spectrum(target, scale, cutoff, primitive=False)
    assert {(r['p'], r['q']) for r in harmonic_rows} == {(355*k, 113*k) for k in (1, 2, 3)}
    # Spectral sensitivity to interval width; no invented precision.
    sensitivity = []
    for digits in [4, 6, 8, 10]:
        m = 10**digits
        rough = (F(floor(target[0]*m), m), F(floor(target[1]*m)+1, m))
        _, rr = finite_spectrum(rough, scale, cutoff)
        sensitivity.append(dict(decimal_places=digits,
                                possible_modes=len(rr),
                                certain_modes=sum(r['certain'] for r in rr)))
    result = dict(status='PASS', exact_case_checks=checks,
                  model='q^2/s+s*(p-alpha*q)^2; q>0',
                  source='Machin alternating rational interval; not pi decimals',
                  scale=scale, cutoff=cutoff, exhaustive_q_cap=cap,
                  primitive_modes=rows, harmonic_modes=harmonic_rows,
                  interval_sensitivity=sensitivity,
                  scope='Known 355/113 is a control, not a new approximation record.')
    out = Path(__file__).with_name('spectral_checks.json')
    out.write_text(json.dumps(encode(result), indent=2)+'\n')
    print(json.dumps(dict(status='PASS', exact_checks=checks, q_cap=cap,
                          modes=[(r['p'], r['q']) for r in rows],
                          quality_upper=float(rows[0]['quality_hi']),
                          interval_sensitivity=sensitivity), indent=2))


if __name__ == '__main__':
    run()
