"""Exact rational checks of modular pulse heights, thresholds and peak spacing.

A rational probe of a rational-function identity suffices for these finite
checks. Targets coinciding with a chart's base point have their own formula.
"""
from fractions import Fraction
from itertools import combinations

def pulse(alpha, x, y):
    q = x.denominator
    return y / (q*q * ((alpha - x)**2 + y*y))

def run():
    alphas = [Fraction(13, 31), Fraction(41, 29), Fraction(-17, 43)]
    slopes = sorted(set(Fraction(p, q) for q in range(1, 13) for p in range(-2, 2*q + 1)))
    heights = [Fraction(1), Fraction(3, 2), Fraction(2), Fraction(5)]
    pulse_count = pair_count = 0
    for alpha in alphas:
        for x in slopes:
            delta = abs(alpha - x)
            if not delta:
                continue
            q = x.denominator
            peak = 1 / (2 * q*q * delta)
            for ratio in (Fraction(1, 5), Fraction(2, 3), Fraction(1), Fraction(3, 2), Fraction(5)):
                y = delta * ratio
                h = pulse(alpha, x, y)
                assert h == peak * 2 * ratio / (1 + ratio**2)
                assert h == pulse(alpha, x, delta / ratio)
                for threshold in heights:
                    radius = 1 / (2 * threshold * q*q)
                    assert (h > threshold) == ((alpha-x)**2 + (y-radius)**2 < radius**2)
                pulse_count += 1
        for x, z in combinations(slopes, 2):
            ex, ez = x - alpha, z - alpha
            if not ex or not ez:
                continue
            hx, hz = abs(ex), abs(ez)
            amplitude_x = 1 / (2 * x.denominator**2 * hx)
            amplitude_z = 1 / (2 * z.denominator**2 * hz)
            area = abs(x.numerator*z.denominator-z.numerator*x.denominator)
            squared_gap_function = (hx + hz if ex*ez < 0 else hx - hz)**2 / (4*hx*hz)
            assert squared_gap_function == area**2 * amplitude_x * amplitude_z
            pair_count += 1
    print(f"Checked {pulse_count} rational pulse samples and {pair_count} pairs of signed errors.")
if __name__ == "__main__":
    run()
