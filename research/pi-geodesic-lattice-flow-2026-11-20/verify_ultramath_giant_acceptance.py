"""Independently check the delivered ZIP, the finite CRT slice, and the new bound."""

import argparse
from fractions import Fraction
from hashlib import file_digest, sha256
from math import gcd
from zipfile import ZipFile


ZIP_DIGEST = "5797e0f7a8f2737b65d8ea1125595cd7a2d8f4671f5b0e4bd681050c65530e44"
PREFIX = "giant_record_round_bundle/"


def check_archive(path):
    with open(path, "rb") as source:
        assert file_digest(source, "sha256").hexdigest() == ZIP_DIGEST
    with ZipFile(path) as archive:
        assert archive.testzip() is None
        manifest = archive.read(PREFIX + "MANIFEST.sha256")
        lines = manifest.splitlines(keepends=True)
        entries = [line.decode().split(maxsplit=1) for line in lines]
        names = {name.strip() for _, name in entries}
        expected = {PREFIX + name for name in names | {"MANIFEST.sha256"}}
        assert len(lines) == len(names) == 416
        assert len(archive.namelist()) == len(expected) == 417
        assert set(archive.namelist()) == expected
        for digest, name in entries:
            assert sha256(archive.read(PREFIX + name.strip())).hexdigest() == digest
        selfhash = archive.read(PREFIX + "MANIFEST.selfhash.txt").split()[0].decode()
        assert lines[-1].endswith(b"  MANIFEST.selfhash.txt\n")
        assert sha256(b"".join(lines[:-1])).hexdigest() == selfhash
        assert sha256(manifest).hexdigest() != selfhash
    print("Outer ZIP, CRC and 416 payload hashes verified. Self-hash covers 415 lines.")


def atan_interval(base):
    terms = 48
    lower = sum(Fraction((-1)**j, (2*j + 1) * base**(2*j + 1)) for j in range(terms))
    upper = lower + Fraction(1, (2*terms + 1) * base**(2*terms + 1))
    return lower, upper


def verify_crt_slice():
    a_lo, a_hi = atan_interval(5)
    b_lo, b_hi = atan_interval(239)
    pi_lo, pi_hi = 16*a_lo - 4*b_hi, 16*a_hi - 4*b_lo
    scale = 1 << 192
    alpha_lo, alpha_hi = scale // pi_hi, -(-scale // pi_lo)
    numerator, r = 113, 355
    seed_upper = r * max(abs(r*alpha_lo - numerator*scale),
                         abs(r*alpha_hi - numerator*scale))
    survivors = []
    for m in range(1, 100001):
        n = r*m
        denom = n*scale
        lower = -((seed_upper - n*n*alpha_lo) // denom)
        upper = (n*n*alpha_hi + seed_upper) // denom
        survivors.extend((m, v) for v in range(lower, upper + 1) if gcd(v, n) == 1)
    assert survivors == [(1, 113)]
    print("100,000 denominators divisible by 355 checked; only the seed survives.")


def verify_dalzell_range():
    index = 383_059_989_912
    digit = 845_188_380_288
    # The rational 1.618033 is less than the golden ratio.
    assert 0 < 2*1618033 - 10**6 and (2*1618033 - 10**6)**2 < 5*10**12
    assert 1618033**1000 > (1 << 694)*10**6000
    r_exponent = 694*((index - 1)//1000)
    assert digit > 1 << 39
    exponent = 2*r_exponent + 39
    # 20736^s < 2^(717s/50) and 24*s < 2^41 on the stated range.
    assert 20736**50 < 1 << 717
    max_sum = ((exponent - 41)*50)//717
    assert max_sum == 37_077_215_113
    assert 24*max_sum < 1 << 41
    assert 717*max_sum + 50*41 <= 50*exponent
    print(f"Dalzell–Lucas enlarged range certified: k+l <= {max_sum:,}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("zip_path", nargs="?", help="Optional original ZIP for hash verification")
    arguments = parser.parse_args()
    if arguments.zip_path:
        check_archive(arguments.zip_path)
    verify_crt_slice()
    verify_dalzell_range()
