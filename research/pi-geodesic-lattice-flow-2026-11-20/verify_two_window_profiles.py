"""Rational certificates for the two window audits. Requires Python 3."""
from fractions import Fraction as F
from math import factorial, gcd


def atan_interval(b, terms):
    s = sum((F((-1)**k, (2*k+1)*b**(2*k+1)) for k in range(terms)), F(0))
    next_term = F((-1)**terms, (2*terms+1)*b**(2*terms+1))
    return min(s,s+next_term), max(s,s+next_term)


def pi_interval(terms):
    a,b = atan_interval(5,terms)
    c,d = atan_interval(239,terms)
    return 16*a-4*d, 16*b-4*c


def abs_interval(lo,hi):
    if lo <= 0 <= hi:
        return F(0),max(-lo,hi)
    return min(abs(lo),abs(hi)),max(abs(lo),abs(hi))


def sinc_polynomial(x, last):
    return sum(((-1)**j*x**(2*j)/factorial(2*j+1) for j in range(last+1)),F(0))


def certify_last_half_height():
    lo,hi = pi_interval(14)
    a,q,Q = 113,355,256
    w0,w1 = abs_interval(1-F(a,q)*hi, 1-F(a,q)*lo)
    assert 0 < w0 <= w1 < F(1,100)
    assert w1/lo < F(1,8*Q**2)  # uniqueness of nearest shell node
    u250 = 250*Q**2*w1
    u251 = 251*Q**2*w0
    assert 0 < u250 < u251 < 2
    # sinc decreases on [0,pi]. Alternating Taylor bounds through x^10, x^8.
    lower = sinc_polynomial(u250,5)
    denom_lower = 1-w1*w1/6
    assert lower>0 and lower**2 > F(5001,10000)
    upper = sinc_polynomial(u251,4) / denom_lower
    assert upper**2 < F(498,1000)
    tail = F(27,251**2)
    assert 8*hi**2/3 < 27 and F(498,1000)+tail < F(1,2)
    assert 1/(lo*denom_lower)**2 + tail < F(1,2)
    print("Q=256: half-height holds exactly at integer H=1,...,250")


def certify_truncation_gap():
    p0,p1 = pi_interval(10)
    scale = 10**9
    lo = F(3141592653,scale)
    hi = F(3141592654,scale)
    assert lo < p0 < p1 < hi
    Q = 512
    count = total = 0
    for q in range(Q,2*Q):
        a = (F(q,1)/lo + F(1,2)).__floor__()
        assert a == (F(q,1)/hi + F(1,2)).__floor__()
        if gcd(q,a) != 1:
            continue
        t0,t1 = abs_interval(q-a*hi,q-a*lo)
        sin_lower = t0-t1**3/6
        if sin_lower>0 and sin_lower**8*Q**3 >= 1:
            contribution = F(q*q,Q**4)/t1**2
            total += (scale*contribution).__floor__()
            count += 1
    discarded_floor = F(total,scale)
    old_upper = 8*hi**2/(9*Q) + F(1,Q**3)
    assert discarded_floor > F(18,1000) > old_upper
    print("Q=512: discarded normalized terms alone exceed 0.018;",
          count,"certified terms, lower sum in billionths:",total)


def certify_secondary_nodes():
    lo,hi = pi_interval(120)
    nodes = [
        (67, 85226204065545021247, 267746016585660808088, 644296),
        (230, 1092952505136224121241866812364478845496740650511249959836121401165445,
         3433611560858522460944982514029752488389302644609030681931573361922021, 1017963),
    ]
    for j,a,e,target in nodes:
        Q=2**j
        assert Q <= e < 2*Q and gcd(a,e)==1
        t0,t1 = abs_interval(e-a*hi, e-a*lo)
        assert 0 < t0 <= t1 < 1 and t1**8*Q**3 < 1
        h0,h1 = F(e*e,Q**4)/t1**2, F(e*e,Q**4)/t0**2
        assert F(target,10**6)<h0<=h1<F(target+1,10**6)
        d = pow(a,-1,e)
        b = (a*d-1)//e
        left, right = F(b,d),F(a-b,e-d)
        assert left<F(a,e)<right and (right < 1/hi or 1/lo < left)
        print(f"Q=2^{j}: non-intermediate node height in ({target},{target+1})/10^6")


if __name__ == "__main__":
    certify_last_half_height()
    certify_truncation_gap()
    certify_secondary_nodes()
