"""Finite checks for the theta-wave note. Run with Python and mpmath."""
import argparse
import hashlib
import sys
from pathlib import Path

from mpmath import iv

def wave_certificates(first=17, last=21, cutoff=14000):
    """Enclose |W(t/2)|/|W(t)|, t=2^-m, by rigorous Gaussian tails."""
    iv.dps = 32
    one = iv.mpf(1)
    times = [one / 2**m for m in range(first, last + 2)]
    real = [3 * one / 2 for t in times]
    imag = [iv.mpf(0) for t in times]
    for n in range(1, cutoff + 1):
        square = n * n
        weight = 1 if n % 2 else 3
        cosine = weight * iv.cos(2 * square)
        sine = weight * iv.sin(2 * square)
        gaussian = iv.exp(-times[-1] * square)
        for i in reversed(range(len(times))):
            real[i] += gaussian * cosine
            imag[i] += gaussian * sine
            gaussian *= gaussian
    norms = []
    for t, x, y in zip(times, real, imag):
        # Paired tail has modulus at most 3 int_N^infty exp(-t x^2) dx.
        tail = 3 * iv.exp(-t * cutoff**2) / (2 * t * cutoff)
        norms.append(iv.sqrt(x**2 + y**2) + iv.mpf([-1, 1]) * tail)
    ratios = [upper / lower for lower, upper in zip(norms, norms[1:])]
    normalized = [iv.sqrt(iv.sqrt(t)) * norm for t, norm in zip(times, norms)]
    return normalized, ratios

def check_large_seed(path):
    """Use the interval certified by the companion GMP program, not a theta sum."""
    data = path.read_bytes()
    expected = "61971ced12e4fb6b02d0cdd55945e93c8b0cfd7d25b5340405ad35360c6cba42"
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError("seed file differs from the certified input")
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    p, q = map(int, data.split())
    scale = 10**26
    low = 558229268763056064281
    high = low + 1
    first = (10 * p * p - 1).bit_length()
    last = (p * q * scale // (16 * high)).bit_length() - 1
    assert (first, last) == (535137, 535144)
    assert p < 4 * q and 8 * high * 10000 < scale
    # Each pair brackets 1000*t/tau. All divisions below are positive.
    for m, lo, hi in [(535145, 7352, 7353), (535146, 3676, 3677),
                      (535147, 1838, 1839)]:
        numerator = 1000 * p * q * scale
        two_power = 1 << (m + 1)
        assert lo * two_power * high < numerator < hi * two_power * low
    print("Ten growth steps at m=535137..535146, followed by a break at 535147.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--witness", type=Path)
    args = parser.parse_args()
    if args.witness:
        check_large_seed(args.witness)
    else:
        heights, ratios = wave_certificates(first=16, last=21)
        for m, height, ratio in zip(range(16, 22), heights, ratios):
            assert (ratio > iv.mpf(4) / 3) if 17 <= m <= 20 else (ratio < iv.mpf(4) / 3)
            assert height > 1
            print(m, height, ratio)
