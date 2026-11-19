"""Certify consecutive sine blocks with rational enclosures."""
from fractions import Fraction
from math import ceil, floor, gcd

SCALE = 10**25


def outward(lower, upper, scale=SCALE):
    return Fraction(floor(lower * scale), scale), Fraction(ceil(upper * scale), scale)


def arctan_bounds(reciprocal, terms):
    """An even alternating partial sum and the next odd partial sum."""
    assert terms % 2 == 0
    lower = sum(Fraction((-1)**k, (2*k+1)*reciprocal**(2*k+1)) for k in range(terms))
    return lower, lower + Fraction(1, (2*terms+1)*reciprocal**(2*terms+1))


def pi_bounds():
    lo5, hi5 = arctan_bounds(5, 24)
    lo239, hi239 = arctan_bounds(239, 12)
    # Machin's identity follows by tangent addition, with the angle in (0, pi/2).
    return outward(16*lo5-4*hi239, 16*hi5-4*lo239)


def sine_bounds(lower, upper):
    lower, upper = outward(lower, upper)
    center, radius = (lower + upper)/2, (upper - lower)/2
    assert abs(center) < 2
    terms = 14
    partial, term = center, center
    for k in range(1, terms):
        term *= -center**2 / ((2*k)*(2*k+1))
        partial += term
    remainder = abs(term*center**2/((2*terms)*(2*terms+1)))
    return outward(partial-remainder-radius, partial+remainder+radius)


def reciprocal_sine_square(k, target, pi_interval):
    center = k*(target[0]+target[1])/2
    integer = floor(center+Fraction(1, 2))
    residuals = [k*a-integer for a in target]
    angles = [p*r for p in pi_interval for r in residuals]
    lower, upper = sine_bounds(min(angles), max(angles))
    assert lower*upper > 0
    magnitudes = sorted((abs(lower), abs(upper)))
    return outward(1/magnitudes[1]**2, 1/magnitudes[0]**2)


def block(q, shift, target, pi_interval):
    values = [reciprocal_sine_square(k+shift*q, target, pi_interval) for k in range(1, q)]
    return tuple(sum(v[index] for v in values)/q**2 for index in (0, 1))


def quantify(name, p, q, target, pi_interval, leading_tolerance):
    assert gcd(p, q) == 1
    eta = sorted(abs(q*q*a-p*q) for a in target)
    assert 0 < eta[0] <= eta[1] <= Fraction(1, 4)
    base = block(q, 0, target, pi_interval)
    next_block = block(q, 1, target, pi_interval)
    diff = next_block[0]-base[1], next_block[1]-base[0]
    ratio = diff[0]/eta[1]**2, diff[1]/eta[0]**2
    assert Fraction(2, 81) < ratio[0] <= ratio[1] < 256
    leading = tuple(2*x*x*Fraction(q**4-1, 15*q**4) for x in pi_interval)
    assert max(abs(a-b) for a in ratio for b in leading) < leading_tolerance
    readable_scale = 10**9
    intervals = (eta, base, diff, ratio)
    print(name, [(floor(lo*readable_scale), ceil(hi*readable_scale)) for lo, hi in intervals])


if __name__ == '__main__':
    pi_interval = pi_bounds()
    phi = tuple(1/x for x in reversed(pi_interval))
    print('Enclosures for |eta|, B_0, B_1-B_0, (B_1-B_0)/eta², in billionths:')
    quantify('reciprocal of pi / 113/355', 113, 355, phi, pi_interval, Fraction(1, 25))
    p, q = 42, 127
    for sign in (1, -1):
        alpha = Fraction(p, q)+Fraction(sign, 10**5*q**2)
        quantify(f'rational target, signed displacement {sign}', p, q, (alpha, alpha),
                 pi_interval, Fraction(1, 5000))
