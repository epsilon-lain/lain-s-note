"""Summarize the existing exact prime table; do not extend the prime scan.

Integer checks audit the table's reflection/count identities. Floating-point
logarithms are descriptive only and certify no asymptotic inequality.
"""
from pathlib import Path
import csv
import hashlib
import json
import math


def main():
    base = Path(__file__).resolve().parent
    source = base.parent / 'pi-initial-value-sieve-2026-11-02/active-primes/active_prime_roots.csv'
    raw = source.read_bytes()
    blob = hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest()
    assert blob == 'f9c06a0b04955a0328d5c5a25616b07135d04629'
    rows = []
    for entry in csv.DictReader(raw.decode().splitlines()):
        p, z = int(entry['p']), int(entry['root_count'])
        roots = [int(s) for s in entry['roots'].split()]
        midpoint = int(entry['middle_is_root'] == 'True')
        assert roots == sorted(set(roots)) and z == len(roots)
        assert all(0 <= r < p for r in roots)
        assert {p-1-r for r in roots} == set(roots)
        assert midpoint == int((p-1)//2 in roots)
        assert (z-midpoint) % 2 == 0
        pairs = (z-midpoint)//2
        active = int(z > 0)
        excess = pairs+midpoint-active
        assert excess >= 0
        assert z+midpoint == 2*(active+excess)
        rows.append((p,z,midpoint,excess))
    snapshots = []
    for x in (1000,2000,5000,10000,20000):
        rr = [r for r in rows if r[0] <= x]
        total = math.fsum(z*math.log(p) for p,z,m,e in rr)
        mid = math.fsum(m*math.log(p) for p,z,m,e in rr)
        active = math.fsum(math.log(p) for p,z,m,e in rr if z)
        loss = math.fsum(e*math.log(p) for p,z,m,e in rr)
        snapshots.append(dict(x=x,zero_first_moment_over_x=total/x,
                              midpoint_weight_over_x=mid/x,
                              active_weight_over_x=active/x,
                              paired_upper_envelope_over_x=(total+mid)/(2*x),
                              envelope_overcount_over_x=loss/x))
    out = dict(status='passed exact count identities',
               source_commit='9ac9c9065a038a87a6ab587b8d87c65f689f393b',
               source_blob_sha1=blob,source_sha256=hashlib.sha256(raw).hexdigest(),
               rows_checked=len(rows),new_prime_classifications=0,
               snapshots=snapshots,
               scope='Existing finite table only. Logarithms are descriptive; no first-moment or capacity limit is proved by these numbers.')
    (base/'existing_pool_moments.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__ == '__main__':
    main()
