#!/usr/bin/env python3
"""Verify fixed snapshot files and record errors using an independent Machin formula.

Usage: python verify_snapshot.py SNAPSHOT_DIRECTORY [OUTPUT_JSON]
Requires gmpy2. No network, no original scan, no invocation of uploaded code.
The index verification is supplied separately by the exact CF reconstruction.
"""
import hashlib
import json
from pathlib import Path
import sys
import time
import gmpy2 as g


def rat(pair):
    return g.mpq(g.mpz(pair[0], 16), g.mpz(pair[1], 16))


def main():
    root = Path(sys.argv[1])
    manifest = json.loads((root / 'manifest.json').read_text())
    for name, metadata in manifest.items():
        data = (root / name).read_bytes()
        assert len(data) == metadata['bytes'], name
        assert hashlib.sha256(data).hexdigest() == metadata['sha256'], name
    records = [json.loads(p.read_text()) for p in sorted((root / 'records').glob('*.json'))]
    bits = max(1024, 2 * max(x['q_bits'] for x in records) + 512)
    start = time.monotonic()
    # Include rounding error of the atan INPUTS as well as the atan outputs.
    with g.context(precision=bits, round=g.RoundDown):
        a5lo = g.atan(g.mpfr(g.mpq(1, 5)))
        a239lo = g.atan(g.mpfr(g.mpq(1, 239)))
    with g.context(precision=bits, round=g.RoundUp):
        a5hi = g.atan(g.mpfr(g.mpq(1, 5)))
        a239hi = g.atan(g.mpfr(g.mpq(1, 239)))
    with g.context(precision=bits, round=g.RoundDown):
        lower = g.mpq(16 * a5lo - 4 * a239hi)
    with g.context(precision=bits, round=g.RoundUp):
        upper = g.mpq(16 * a5hi - 4 * a239lo)
    assert 3 < lower < upper < 4
    rows = []
    for rec in records:
        p, q = g.mpz(rec['p_hex'], 16), g.mpz(rec['q_hex'], 16)
        assert q > 0 and g.gcd(p, q) == 1
        left, right = q * lower - p, q * upper - p
        assert not left <= 0 <= right
        qlo = q * min(abs(left), abs(right))
        qhi = q * max(abs(left), abs(right))
        assert rat(rec['quality_lower_hex']) <= qlo <= qhi <= rat(rec['quality_upper_hex'])
        idx = rec['convergent_index_zero_based']
        determinant = -1 if idx % 2 == 0 else 1
        previous_q = g.mpz(0) if q == 1 else (determinant * g.invert(p, q)) % q
        previous_p = (p * previous_q - determinant) // q
        assert p * previous_q - previous_p * q == determinant
        endpoints = [(previous_p - x * previous_q) / (x * q - p) for x in (lower, upper)]
        tlo, thi = min(endpoints), max(endpoints)
        a = g.mpz(rec['following_partial_quotient'])
        assert g.floor(tlo) == g.floor(thi) == a
        with g.context(precision=100):
            rows.append({
                'index': idx, 'q_bits': q.bit_length(),
                'q_decimal_digits': len(q.digits(10)),
                'quality': str(g.mpfr((qlo + qhi) / 2)),
                'a': str(a), 'beta': str(g.mpfr(g.mpq(previous_q, q))),
                'theta': str(g.mpfr((tlo + thi) / 2 - a)),
                'quality_certificate_verified': True,
                'following_partial_quotient_verified': True,
            })
    result = {'snapshot_manifest_entries_verified': len(manifest),
              'pi_formula': '16 atan(1/5) - 4 atan(1/239); directed rounding including inputs',
              'precision_bits': bits, 'seconds': time.monotonic() - start,
              'gmpy2_version': g.version(), 'records': rows,
              'scope': 'Verifies the supplied candidates, not full reported 10.18m coverage.'}
    if len(sys.argv) > 2:
        Path(sys.argv[2]).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
