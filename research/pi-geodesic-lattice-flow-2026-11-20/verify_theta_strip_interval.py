"""Independently enclose the signed near-axis union for the fixed theta angle."""
from mpmath import iv

iv.dps = 30
I = iv.mpf
pi = iv.pi
sqrt = iv.sqrt


def kernel(h, t, shift, alternating):
    frequency = I(2*h)/pi + shift
    rate = pi*pi/(2*t)
    lower_index = int(float(frequency.mid))
    if not (I(lower_index) <= frequency and frequency < I(lower_index+1)):
        raise ArithmeticError("uncertified floor", frequency, lower_index)
    two_terms = I(0)
    for j in (lower_index, lower_index+1):
        sign = -1 if alternating and j % 2 else 1
        two_terms += sign * iv.exp(-rate*(frequency-j)**2)
    tail_first = iv.exp(-rate)
    return two_terms + I([-1,1]) * (2*tail_first/(1-tail_first))


def row(h, t):
    unshifted = kernel(h, t, 0, h % 2 == 1)
    if h % 2:
        weighted = I(3)/4 * unshifted
    else:
        shifted = kernel(h, t, I(1)/2, False)
        weighted = I(5)/4 * unshifted + (-1)**(h//2) * shifted
    return iv.exp(-t*h*h/2) * sqrt(pi/(2*t)) * weighted


def weight(n):
    return I(1)/2 if n % 2 else I(3)/2


def summand(h, k, t):
    n, m = (h+k)//2, (h-k)//2
    u = t*(h*h+k*k)/4
    return weight(n)*weight(m)*iv.cos(I(2*h*k))*(16*iv.exp(-2*u)-9*iv.exp(-u))


def signed_union(t, radius):
    strip = I(0)
    for h in range(1, radius+1):
        strip += 2*(16*row(h,t)-9*row(h,t/2))
        if h % 2 == 0:
            strip -= 2*summand(h,0,t)
    overlap = I(0)
    for h in range(1, radius+1):
        for k in range(1, radius+1):
            if h % 2 == k % 2:
                overlap += 4*summand(h,k,t)
    return sqrt(t)*(2*strip-overlap)


if __name__ == '__main__':
    for denominator in (16,10):
        for exponent in (12,16,20,24,28):
            radius = 1
            while (radius+1)**denominator <= 2**exponent:
                radius += 1
            print(denominator, exponent, radius, signed_union(I(1)/2**exponent, radius))
    target = signed_union(I(1)/2**20, 4)
    assert I('-0.142026071210') < target < I('-0.142026071209')
