"""Exact finite certificates + high-precision checks; not a proof of P2848.

Run in this directory: python verify.py
Optional numerical checks require mpmath. Exact certificates use the standard library.
"""
from fractions import Fraction as Q
from math import gcd
import json


def atan_interval(k, terms=400):
    s = sum((Q((-1)**j, (2*j+1)*k**(2*j+1)) for j in range(terms)), Q(0))
    r = Q((-1)**terms, (2*terms+1)*k**(2*terms+1))
    return min(s, s+r), max(s, s+r)


def pi_interval():
    a, b = atan_interval(5)
    c, d = atan_interval(239)
    return 16*a-4*d, 16*b-4*c


def depth_interval(p, q, lo, hi):
    # h(b)=p^2/b-b*q^2 is strictly decreasing for b>0.
    ends = [Q(p*p)/hi-hi*q*q, Q(p*p)/lo-lo*q*q]
    assert ends[0]*ends[1] > 0
    return min(map(abs, ends)), max(map(abs, ends))


def bezout(a, b):
    r0, r1, s0, s1, t0, t1 = a, b, 1, 0, 0, 1
    while r1:
        k = r0//r1
        r0, r1, s0, s1, t0, t1 = r1, r0-k*r1, s1, s0-k*s1, t1, t0-k*t1
    assert r0 == 1
    return s0, t0


def decimal_enclosure(lo, hi, places=45):
    """Outward-rounded decimal strings, certified using integers only."""
    scale = 10**places
    lower = lo.numerator*scale//lo.denominator
    upper = -((-hi.numerator*scale)//hi.denominator)
    def render(n):
        return f'{n//scale}.{n%scale:0{places}d}'
    return [render(lower), render(upper)]


def exact_checks():
    lo, hi = pi_interval()
    # Existing repository witness, CF index 430; not a new approximation record.
    p = int('1901870728566923076090143944714770339621590768313546337192526115562704339680963564320007808107929370299752345187688835741387003036853361285671158059867702399073227994426905220194699766118756059055619036488502928002591')
    q = int('605384255146420326102361023215940531716391478150345020739231253172134740688232476946000058713774549796561447468267746412874022717544100946587144148739626803435133473281606663121381125761746030151344353855924025288111')
    assert gcd(p, q) == 1
    ml, mh = depth_interval(p, q, lo, hi)
    assert mh < lo/Q(180**2)
    _, small_hi = depth_interval(355, 113, lo, hi)
    assert small_hi < lo/Q(21**2)
    # Exact norm preservation and unimodularity of the two Pell transformations.
    checks = 0
    for D, a, b in [(2, 1, 1), (3, 2, 1)]:
        det = a*a-D*b*b
        assert abs(det) == 1
        for x in range(-30, 31):
            for y in range(-30, 31):
                xx, yy = a*x+D*b*y, b*x+a*y
                assert xx*xx-D*yy*yy == det*(x*x-D*y*y)
                checks += 1
    return lo, hi, {'pell_norm_checks': checks,
                    'old_index_430_theta_at_balance_strictly_above': 180,
                    '355_over_113_theta_at_balance_strictly_above': 21,
                    'index_430_depth_interval_outward_decimal': decimal_enclosure(ml, mh)}


def numerical_checks(lo, hi):
    import mpmath as mp
    mp.mp.dps = 180
    pi = (mp.mpf(lo.numerator)/lo.denominator + mp.mpf(hi.numerator)/hi.denominator)/2

    def vector(beta, t, p, q):
        return (beta**t*(p-beta*q)/mp.sqrt(2*beta),
                beta**(-t)*(p+beta*q)/mp.sqrt(2*beta))

    def theta_reduced(beta, t, p, q):
        s, k = bezout(p, q)
        # p*s + q*k=1, hence det((p,q),(-k,s))=1.
        v, w = vector(beta, t, p, q), vector(beta, t, -k, s)
        f = sum(x*x for x in v)
        u = sum(x*y for x, y in zip(v, w))/f
        u -= mp.nint(u)
        assert abs(v[0]*w[1]-v[1]*w[0]-1) < mp.mpf('1e-140')
        B = int(mp.ceil(mp.sqrt(450*f)))+1
        K = int(mp.ceil(mp.sqrt(450*f)/pi))+1
        # Poisson transform in the coordinate parallel to v; cutoff checks only.
        val = mp.sqrt(pi/f)*mp.fsum(
            mp.exp(-j*j/f)*(1+2*mp.fsum(
                mp.exp(-pi*pi*k*k/f)*mp.cos(2*pi*k*j*u) for k in range(1, K+1)))
            for j in range(-B, B+1))
        return f, val

    rows = []
    for beta, p, q, pp, qp in [(mp.sqrt(2), 7, 5, 3, 2),
                              (mp.sqrt(3), 7, 4, 5, 3),
                              (pi, 22, 7, 3, 1), (pi, 355, 113, 333, 106)]:
        t0 = mp.log(abs((p+beta*q)/(p-beta*q)))/(2*mp.log(beta))
        m = abs(p*p-beta*beta*q*q)/beta
        alpha = (pp-beta*qp)/(beta*q-p)
        cf_mass = alpha+mp.mpf(qp)/q
        assert abs(1/m-cf_mass/(2+(mp.mpf(p)/q-beta)/beta)) < mp.mpf('1e-130')
        for h in [mp.mpf('0'), mp.mpf('0.25'), mp.mpf('-0.25')]:
            f, theta = theta_reduced(beta, t0+h, p, q)
            assert abs(f-m*mp.cosh(2*h*mp.log(beta))) < mp.mpf('1e-140')
            lead = mp.sqrt(pi/f)
            error = lead*2*mp.exp(-pi*pi/f)/(1-mp.exp(-3*pi*pi/f))
            error += 2*(1+lead)*mp.exp(-1/f)/(1-mp.exp(-3/f))
            assert -mp.mpf('1e-140') < theta-lead < error+mp.mpf('1e-140')
            rows.append({'base': mp.nstr(beta, 20), 'p': p, 'q': q, 'offset': str(h),
                         'depth': mp.nstr(f, 25), 'theta_numeric': mp.nstr(theta, 25),
                         'sqrt_pi_over_depth': mp.nstr(lead, 25),
                         'analytic_error_bound_numeric': mp.nstr(error, 12)})
    # Full theta periodicity checked at several points using reindexed bases.
    for D, a, b in [(2, 1, 1), (3, 2, 1)]:
        beta = mp.sqrt(D)
        period = mp.log(a+b*beta)/mp.log(beta)
        for t in map(mp.mpf, ['0', '.3', '1.1']):
            _, first = theta_reduced(beta, t, 1, 0)
            _, second = theta_reduced(beta, t+period, a, b)
            assert abs(first-second) < mp.mpf('1e-130')
    return rows


if __name__ == '__main__':
    lo, hi, exact = exact_checks()
    try:
        numeric = numerical_checks(lo, hi)
    except ImportError:
        numeric = 'Skipped: install mpmath for optional numerical checks.'
    print(json.dumps({'exact_certificates': exact, 'numerical_checks_not_interval_proofs': numeric,
                      'status': 'No proof or disproof of P2848; no Lean verification.'}, indent=2))
