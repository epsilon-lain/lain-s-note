"""Three crossing marks, decoded using outward-rounded rational intervals only."""
from fractions import Fraction as F
from math import ceil, floor


def prefix_value(digits, tail):
    value = tail
    for digit in reversed(digits):
        value = digit + 1 / value
    return value


def convergents(digits):
    before, previous = (0, 1), (1, digits[0])  # (denominator, numerator)
    yield previous
    for digit in digits[1:]:
        current = tuple(x + digit*y for x, y in zip(before, previous))
        yield current
        before, previous = previous, current


def crossing(left, right, target):
    q, p = left
    Q, P = right
    h, H = p - target*q, P - target*Q
    assert h*H < 0
    # Linear interpolation of the actual segment, with no inverse formulas.
    return (q*H - Q*h)/(H-h)


def outward(lo, hi, scale):
    return F(floor(lo*scale), scale), F(ceil(hi*scale), scale)


def mark_intervals(digits, target_bounds, scale):
    points = list(convergents(digits))
    marks = []
    for v, w in zip(points, points[1:]):
        # This fractional-linear function is monotone on the enclosure.
        for node in (v,w):
            q,p = node
            assert (p-q*target_bounds[0])*(p-q*target_bounds[1]) > 0
        bounds = sorted(crossing(v, w, target) for target in target_bounds)
        marks.append(outward(*bounds, scale))
    return marks


def candidates(left, right):
    low = right[0]/left[1] - 1
    high = right[1]/left[0] - 1
    return range(max(1, floor(low)+1), ceil(2*high))


def decode(marks):
    left, center, right = marks
    u, v, w = [(1/hi, 1/lo) for lo, hi in marks]
    D = u[0]-v[1], u[1]-v[0]
    E = v[0]-w[1], v[1]-w[0]
    pairs = [(d,e) for d in candidates(left, center) for e in candidates(center, right)]
    certified = []
    for d, e in pairs:
        lower = D[0]/d + E[0]/e - v[1]
        upper = D[1]/d + E[1]/e - v[0]
        if lower <= 0 <= upper:
            certified.append((d,e))
    return len(pairs), certified


def finite_geometry_check(digits, target):
    # Rational test ray with the prescribed continued-fraction prefix and
    # a nonintegral tail beyond it; all tested entries precede termination.
    vertices = list(convergents(digits))
    marks = [crossing(a,b,target) for a,b in zip(vertices,vertices[1:])]
    b = [abs(p-target*q) for q,p in vertices]
    lobes = [b[n]*(marks[n]-marks[n-1])/2 for n in range(1,len(marks))]
    tau = [b[n]/(b[n]+b[n+1])-F(1,2) for n in range(len(marks))]
    for n, area in enumerate(lobes,1):
        gap = marks[n]-marks[n-1]
        assert F(1,12) < area < F(1,2)
        assert gap**2/(2*area*marks[n]*marks[n-1]) == digits[n+1]
        assert 2*area == tau[n-1]+tau[n]
    for start in range(len(lobes)):
        signed = F(0)
        for stop in range(start,len(lobes)):
            signed += (-1)**(stop-start)*lobes[stop]
            assert 2*signed == tau[start] + (-1)**(stop-start)*tau[stop+1]
    return len(lobes)


if __name__ == '__main__':
    phi_lo, phi_hi = F(1618033,10**6), F(1618034,10**6)
    assert phi_lo**2-phi_lo-1 < 0 < phi_hi**2-phi_hi-1
    digits = [0, 2, 17, 9]
    bounds = sorted(prefix_value(digits, t) for t in (phi_lo,phi_hi))
    scale = 10**5
    marks = mark_intervals(digits, bounds, scale)
    count, recovered = decode(marks)
    assert count > 1 and recovered == [(17,9)]
    print('Crossing enclosures: all endpoints divided by', scale,
          [(int(lo*scale), int(hi*scale)) for lo,hi in marks])
    print('Candidate pairs before/after the reciprocal-mark test:', count, recovered)
    patterns = [[0]+[1]*15, [3,7,15,1,292,1,1,1,2,1,3,1,14],
                [0]+[x for n in range(1,11) for x in (1,n*n)]]
    tested = sum(finite_geometry_check(prefix, prefix_value(prefix,F(19,8))) for prefix in patterns)
    print('Exact finite-ray lobe and signed interval checks:', tested, 'lobes')
