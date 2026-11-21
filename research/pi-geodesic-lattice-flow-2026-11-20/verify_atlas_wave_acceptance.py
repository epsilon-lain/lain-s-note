"""Verify the delivered ZIP, the exact counterexamples and archived check groups.

Usage: python3 verify_atlas_wave_acceptance.py /path/to/giant_atlas_wave_bundle.zip
This checks the recorded logs. It does not run the programs inside the ZIP.
"""
import argparse
import ast
from collections import Counter
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path
from zipfile import ZipFile


def verify_bundle(path):
    expected = "d79409971a8a0f229566655a4af7c69c4e3c169a4c57c89b3c5045abdaa0be0e"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
    prefix = "giant_atlas_wave_bundle/"
    with ZipFile(path) as bundle:
        names = bundle.namelist()
        assert len(names) == len(set(names)) == 602 and bundle.testzip() is None
        manifest_name = prefix + "MANIFEST.sha256"
        manifest = bundle.read(manifest_name).decode("utf-8").splitlines()
        paths = []
        for line in manifest:
            digest, relative_path = line.split("  ", 1)
            member = prefix + relative_path
            assert hashlib.sha256(bundle.read(member)).hexdigest() == digest
            paths.append(member)
        assert len(paths) == len(set(paths)) == 601
        assert set(paths) == set(names) - {manifest_name}
        scripts = [n for n in names if n.endswith(".py")]
        for name in scripts:
            ast.parse(bundle.read(name), filename=name)
        assert len(scripts) == 78
        outcomes = json.loads(bundle.read(prefix + "VERDICTS.json"))
        a, b = (outcomes["swarms"][group]["directions"]
                for group in ("smwq2qriq", "smwq2qrl4"))
        assert Counter(x["verdict"] for x in a) == {"accepted": 14, "accepted_with_errata": 2}
        assert Counter(x["verdict"] for x in b) == {"accepted": 16}
        log = bundle.read(prefix + "coordinator/checks/packaged/smoke.final.driver.log")
        groups = [line.removeprefix(b"PASS: ") for line in log.splitlines() if line.startswith(b"PASS: ")]
        assert len(groups) == len(set(groups)) == 20
        assert sum(group.endswith(b"_lean") for group in groups) == 2
    print("Verified ZIP, 601 payload hashes, 78 Python sources, 32 verdicts and 20 named test groups.")


def verify_counterexamples():
    r, c, d, p, q = 2, 1, 1, 1, 2
    common = gcd(r*p+c, r*q+d)
    output = Fraction(r*p+c, r*q+d)
    s = Fraction(common, output.denominator)
    assert r*s == Fraction(r*common**2, r*q+d) == Fraction(2, 5)
    assert Fraction(common**2, q+1) == Fraction(1, 3)

    # For a=22, r=69, |sin(a)| <= |22-7*pi|.
    lo, hi = Fraction(314159, 100000), Fraction(314160, 100000)
    assert max(abs(22-7*lo), abs(22-7*hi)) < Fraction(9, 1000)
    assert 22-69/lo > Fraction(36, 1000)
    assert max(abs(22*lo-69), abs(22*hi-69)) < lo/2
    print("Verified the wrong affine gain and wrong sine argument by explicit counterexamples.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    verify_bundle(parser.parse_args().bundle)
    verify_counterexamples()
