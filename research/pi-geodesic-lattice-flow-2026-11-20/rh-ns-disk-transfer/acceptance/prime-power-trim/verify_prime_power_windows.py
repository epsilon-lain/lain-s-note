"""Rebuild hull areas and check the single-residue criterion against a full search."""

import argparse
import json
from math import gcd, isqrt
from pathlib import Path


def twice_hull_area(radius):
    chain = []
    for x in range(radius + 1):
        y = isqrt(radius**2 - x**2)
        while len(chain) > 1:
            ax, ay = chain[-2]
            bx, by = chain[-1]
            if (bx - ax) * (y - ay) < (by - ay) * (x - ax):
                break
            chain.pop()
        chain.append((x, y))
    return 4 * sum(ay * bx - ax * by for (ax, ay), (bx, by) in zip(chain, chain[1:]))


def common_modulus(prime, exponent):
    largest_power = prime ** (2 * exponent // 3)
    scale = prime**exponent
    common_denominator = largest_power if prime == 2 else 2 * largest_power
    return 2 * scale**2 // common_denominator


def small_enough_denominator(numerator, denominator, radius):
    reduced_denominator = denominator // gcd(numerator, denominator)
    return reduced_denominator**3 <= radius**2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()
    data = json.loads(args.certificate.read_text(encoding="utf-8"))
    count = 0
    unique_count = 0
    for group in ("windows_H1", "windows_H40"):
        for row in data[group]:
            radius = row["N"]
            area = twice_hull_area(radius)
            assert area == row["M"]
            denominator = 2 * radius**2
            limit = row["actual_d_upper"]
            step = common_modulus(row["ell"], row["n"])
            compressed = [
                d for d in range(area % step, limit + 1, step)
                if small_enough_denominator(area - d, denominator, radius)
            ]
            exhaustive = [
                d for d in range(limit + 1)
                if small_enough_denominator(area - d, denominator, radius)
            ]
            assert compressed == exhaustive
            if step > limit:
                unique_count += 1
                assert len(compressed) <= 1
            assert bool(compressed) == bool(row["hit_count"])
            count += 1

    radius = 81
    area = twice_hull_area(radius)
    modulus = common_modulus(3, 4)
    trim = area % modulus
    local = 8 * 404 * 2 % modulus
    assert (area, modulus, trim, local) == (41064, 729, 240, 632)
    assert trim**3 < 20**3 * radius**2 < local**3
    assert gcd(area - trim, 2 * radius**2) == 1458
    print(f"Verified all {count} prime-power windows, including {unique_count} single-residue windows.")
    print("N=81, H=20: global trim 240 hits; the local residue 632 would miss.")


if __name__ == "__main__":
    main()
