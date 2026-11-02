"""Exact root-family audit. Python standard library; optional mpmath cross-check.

Run: python3 verify.py
All acceptance decisions use integers/Fraction, never floating point.
"""
from fractions import Fraction as F
from math import gcd, factorial
from pathlib import Path
from decimal import Decimal, localcontext
import json

HERE = Path(__file__).resolve().parent


def atan_bounds(base, terms):
    # Integrating the finite geometric sum gives the alternating remainder bound.
    v = F(0)
    for j in range(terms):
        v += F((-1) ** j, (2*j + 1) * base ** (2*j + 1))
    nxt = F((-1) ** terms, (2*terms + 1) * base ** (2*terms + 1))
    return min(v, v+nxt), max(v, v+nxt)


def pi_bounds():
    # Deliberately different term counts from the independent agent certificate.
    a, b = atan_bounds(5, 220)
    c, d = atan_bounds(239, 75)
    return 16*a-4*d, 16*b-4*c


def pair(v):
    return [str(v.numerator), str(v.denominator)]


def decimal(v, digits=45):
    with localcontext() as ctx:
        ctx.prec = digits
        return str(Decimal(v.numerator)/Decimal(v.denominator))


def quality_bounds(p, q, lo, hi):
    a, b = q*q*lo-p*q, q*q*hi-p*q
    lower = F(0) if a <= 0 <= b else min(abs(a), abs(b))
    return lower, max(abs(a), abs(b))


def common_cf_prefix(lo, hi):
    out = []
    while True:
        a, b = lo.numerator//lo.denominator, hi.numerator//hi.denominator
        if a != b:
            return out
        out.append(a)
        l, h = lo-a, hi-a
        if l == 0:
            return out
        lo, hi = 1/h, 1/l


def convergents(prefix):
    p0, p1, q0, q1 = 0, 1, 1, 0
    for a in prefix:
        p0, p1 = p1, a*p1+p0
        q0, q1 = q1, a*q1+q0
        yield p1, q1


def finite_certificate(lo, hi, bound):
    prefix = common_cf_prefix(lo, hi)
    rows, records = [], []
    best = None
    for index, (p, q) in enumerate(convergents(prefix)):
        assert gcd(p, q) == 1
        if q > bound:
            assert rows and rows[-1]['q_int'] <= bound
            for r in rows:
                r.pop('q_int')
            return {'bound': str(bound), 'common_prefix_length': len(prefix),
                    'checked_convergents': len(rows), 'rows': rows,
                    'records': records, 'first_excluded_convergent':
                    {'index': index, 'p': str(p), 'q': str(q)}}
        l, u = quality_bounds(p, q, lo, hi)
        row = {'index': index, 'p': str(p), 'q': str(q), 'q_int': q,
               'quality_lower': pair(l), 'quality_upper': pair(u),
               'quality_display': decimal((l+u)/2)}
        rows.append(row)
        if best is None or u < best[0]:
            best = (l, u, p, q)
            records.append({'p': str(p), 'q': str(q), 'quality_display': decimal((l+u)/2)})
        else:
            assert l > best[1], 'Intervals do not separate candidate from running record'
    raise AssertionError('Insufficient pi interval to certify all convergents up to bound')


def independent_scan(lo, hi, limit):
    # Convert to an outward-rounded decimal grid, then scan integer residuals.
    scale = 10**90
    lower = lo.numerator*scale//lo.denominator
    upper = -((-hi.numerator*scale)//hi.denominator)
    assert F(lower, scale) <= lo < hi <= F(upper, scale)
    assert upper-lower <= 2
    best_lower = best_upper = None
    records, primitive_count = [], 0
    for q in range(1, limit+1):
        p = (2*q*lower+scale)//(2*scale)
        assert p == (2*q*upper+scale)//(2*scale), 'Nearest numerator unresolved'
        assert 3*q <= p <= 4*q
        if gcd(p, q) != 1:
            continue
        primitive_count += 1
        a, b = q*lower-p*scale, q*upper-p*scale
        l = 0 if a <= 0 <= b else q*min(abs(a), abs(b))
        u = q*max(abs(a), abs(b))
        if best_lower is None or u < best_lower:
            best_lower, best_upper = l, u
            records.append({'p': str(p), 'q': str(q),
                            'quality_lower': pair(F(l, scale)),
                            'quality_upper': pair(F(u, scale))})
        else:
            assert l > best_upper
    assert [(r['p'],r['q']) for r in records] == [('3','1'),('22','7'),('355','113')]
    return {'limit': limit, 'primitive_nearest_candidates': primitive_count,
            'records': records,
            'coverage_argument': 'Any primitive fraction with quality < 1/2 has nearest numerator; nonprimitive nearest pairs recur at a smaller denominator.'}


def check_root_algebra():
    checks = 0
    for q in range(1, 31):
        for k in range(q+1):
            r = F(3*q+k, q)
            assert all(q*r-3*q-n != 0 for n in range(q+1) if n != k)
            derivative = F(q)
            for n in range(q+1):
                if n != k:
                    derivative *= q*r-3*q-n
            assert derivative == q*(-1)**(q-k)*factorial(k)*factorial(q-k)
            checks += 1
    # Exact divisibility/old-root retention under denominator dilation.
    dilation_checks = 0
    for q in range(1, 11):
        for m in range(2, 6):
            x = F(22, 7)
            def poly(d):
                v = F(1)
                for n in range(d+1):
                    v *= d*(x-3)-n
                return v
            rhs = m**(q+1)*poly(q)
            for n in range(m*q+1):
                if n % m:
                    rhs *= m*q*(x-3)-n
            assert poly(m*q) == rhs
            dilation_checks += 1
    return {'distinct_root_derivatives': checks, 'dilation_identities': dilation_checks}


def optional_numerical_checks():
    try:
        import mpmath as mp
    except ImportError:
        return {'skipped': 'mpmath is optional and not required by exact certificate'}
    mp.mp.dps = 100
    checks = 0
    for x in (mp.pi, 2+mp.sqrt(2), mp.mpf(10)/3):
        for q in (1, 2, 7, 20, 31, 113):
            t = q*(x-3)
            product = mp.fprod(t-n for n in range(q+1))
            rhs = (-1)**q*mp.gamma(t+1)*mp.gamma(q+1-t)*mp.sin(mp.pi*t)/mp.pi
            assert mp.almosteq(product, rhs, rel_eps=mp.mpf('1e-80'))
            delta = abs(t-mp.nint(t))
            signal = abs(mp.sin(mp.pi*t))
            assert 2*delta <= signal+mp.mpf('1e-90')
            assert signal <= mp.pi*delta+mp.mpf('1e-90')
            a=x-3
            rho=a**a*(1-a)**(1-a)
            assert abs(product)/mp.factorial(q) <= (q+1)*rho**q/mp.pi
            checks += 1
    return {'nonrigorous_cross_checks': checks, 'note': 'Proofs are analytic; numerical checks are not acceptance certificates.'}


def main():
    lo, hi = pi_bounds()
    assert lo < hi and hi-lo < F(1, 10**300)
    finite = finite_certificate(lo, hi, 10**100)
    assert [(r['p'],r['q']) for r in finite['records']] == [('3','1'),('22','7'),('355','113')]
    best_lo, best_hi = quality_bounds(355, 113, lo, hi)
    assert best_hi < F(1, 2)
    # Stronger than running-record checks: every other convergent is above final best.
    for row in finite['rows']:
        if (row['p'], row['q']) != ('355', '113'):
            assert F(*map(int, row['quality_lower'])) > best_hi
    result = {'pi_lower': pair(lo), 'pi_upper': pair(hi),
              'pi_interval_width_display': decimal(hi-lo),
              'machin_terms': [220,75], 'finite_certificate': finite,
              'brute_force': independent_scan(lo, hi, 10**6),
              'root_algebra': check_root_algebra(),
              'analytic_identity_cross_checks': optional_numerical_checks(),
              'limits': 'Finite certificate plus exact identities only. No new pi record and no infinite quality-zero construction.'}
    (HERE/'verification.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'passed': True, 'bound': finite['bound'],
                      'checked_convergents': finite['checked_convergents'],
                      'records': finite['records'], 'root_algebra': result['root_algebra'],
                      'brute_force_candidates': result['brute_force']['primitive_nearest_candidates'],
                      'pi_interval_width': result['pi_interval_width_display']}, ensure_ascii=False))


if __name__ == '__main__':
    main()
