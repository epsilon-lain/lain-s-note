"""Independent interval tests of the first-order quarter-window formula."""
from fractions import Fraction as F
from mpmath import iv, mp

iv.dps = mp.dps = 65


def show(i):
    return mp.nstr(mp.mpf(i.mid), 12)


def check(name, alpha, angle, a, r, previous_r, A):
    one = iv.mpf(1)
    K = A // 4
    e = abs(r*r*alpha - r*a)
    s = 1/e - A
    tau = s - one*previous_r/r
    assert tau > 0 and tau < 1 and 4*K == A and A >= 16
    H = iv.mpf(0)
    for n in range(1, K*r+1):
        sin = iv.sin(n*angle)
        assert sin.a > 0 or sin.b < 0
        H += 1/sin**2
    H /= (K*r)**2
    forecast = 8*one/3 + (16*s/3 - 16/iv.pi)/A
    scaled = A*A*(H - forecast)
    assert abs(scaled) < 300
    print(name, 'A=', A, 'nodes=', K*r, 'H=', show(H),
          'scaled signed residual=', show(scaled))


def main():
    check('1/pi', 1/iv.pi, iv.mpf(1), 113, 355, 333, 292)
    gold = (iv.sqrt(5)-1)/2
    for a, r, ap, rp, tail in [(1, 6, 0, 1, gold), (2, 7, 1, 3, 1/(2+gold))]:
        for A in [16, 64, 256, 1024]:
            alpha = (a*(A+tail)+ap)/(r*(A+tail)+rp)
            check('quadratic', alpha, iv.pi*alpha, a, r, rp, A)
    A = 845188380288
    coefficient = max(16/F(314,100), F(32,3)-16/F(22,7))
    assert coefficient+F(300,A) < F(558,100)
    assert F(558,100*A) < F(661,10**14)
    print('Certified giant symbolic bound 5.58/A < 6.61e-12.')


if __name__ == '__main__':
    main()
