#!/usr/bin/env python3
"""Exact root scans for U_r mod p, with p odd and 0<=r<p.

Holt reflection reduces the production scan to 0<=r<=(p-1)//2.
All recurrence operations are integer arithmetic modulo p. Logarithms are
used only for descriptive mass summaries, not for classifying primes.
"""

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
from time import perf_counter


def primes_to(n):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, math.isqrt(n) + 1):
        if sieve[p]:
            sieve[p*p::p] = b"\x00" * (((n-p*p)//p)+1)
    return [p for p in range(3, n + 1, 2) if sieve[p]]


def values_mod_prime(p, last):
    inv = [0, 1] + [0] * max(0, last-1)
    for k in range(2, last+1):
        inv[k] = -(p//k) * inv[p % k] % p
    values = [1, 2]
    for k in range(2, last+1):
        values.append((2*(2*k-1)*values[-1]+4*(k-1)*values[-2])*inv[k] % p)
    return values[:last+1]


def roots_half(p):
    half = (p-1)//2
    vals = values_mod_prime(p, half)
    low = [r for r, val in enumerate(vals) if val == 0]
    roots = sorted(set(low + [p-1-r for r in low]))
    return roots


def weighted_stats(records):
    theta = sum(math.log(row["p"]) for row in records)
    active = [row for row in records if row["root_count"]]
    theta_z = sum(math.log(row["p"]) for row in active)
    return {
        "prime_count": len(records),
        "active_count": len(active),
        "active_count_ratio": len(active)/len(records) if records else None,
        "theta_odd": theta,
        "theta_active": theta_z,
        "theta_active_over_theta": theta_z/theta if theta else None,
    }


def run(limit, output_dir):
    started = perf_counter()
    records = []
    reflection_checks = 0
    for p in primes_to(limit):
        roots = roots_half(p)
        if p <= 503:
            full = values_mod_prime(p, p-1)
            assert roots == [r for r, val in enumerate(full) if val == 0]
            # Check the full Holt congruence as well as just zero symmetry.
            for r in range((p-1)//2+1):
                assert full[p-1-r] == pow(-4, (p-1)//2-r, p)*full[r] % p
                reflection_checks += 1
        records.append({
            "p": p,
            "root_count": len(roots),
            "first_root": roots[0] if roots else None,
            "middle_is_root": (p-1)//2 in roots,
            "roots": roots,
        })
    snapshots = []
    for x in [50,100,200,500,1000,2000,5000,10000,20000]:
        if x <= limit:
            stat = weighted_stats([r for r in records if r["p"] <= x])
            stat.update({"x": x, "theta_active_over_x": stat["theta_active"]/x})
            snapshots.append(stat)
    groups = {}
    for mod in [8,16]:
        groups[str(mod)] = {
            str(a): weighted_stats([r for r in records if r["p"] % mod == a])
            for a in range(1,mod,2)
        }
    running_theta_z = 0.0
    prefix_maxima = {lo: {"p": None, "theta_active_over_p": 0.0}
                     for lo in [0,1000,5000,10000] if lo <= limit}
    for row in records:
        if row["root_count"]:
            running_theta_z += math.log(row["p"])
            for lo, maximum in prefix_maxima.items():
                ratio = running_theta_z/row["p"]
                if row["p"] >= lo and ratio > maximum["theta_active_over_p"]:
                    maximum.update({"p": row["p"], "theta_active_over_p": ratio})
    summary = {
        "limit": limit,
        "status": "PASS",
        "elapsed_seconds": perf_counter()-started,
        "full_holt_congruence_checks": reflection_checks,
        "scope": "Finite exact classification, not an asymptotic density theorem.",
        "source": "https://cs.uwaterloo.ca/journals/JIS/VOL9/Noe/noe35.pdf, equations (13), (14), Theorems 8.5 and 8.6",
        "sequence": "U_n=[x^n](1+2x+2x^2)^n; discriminant -4",
        "snapshots": snapshots,
        "groups": groups,
        "prefix_ratio_maxima_by_lower_cutoff": prefix_maxima,
        "root_count_distribution": dict(sorted(Counter(r["root_count"] for r in records).items())),
        "first_active_primes": [r["p"] for r in records if r["root_count"]][:30],
        "first_inactive_primes": [r["p"] for r in records if not r["root_count"]][:30],
        "middle_root_primes": [r["p"] for r in records if r["middle_is_root"]],
    }
    output_dir.mkdir(exist_ok=True,parents=True)
    (output_dir/"active_prime_scan.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    (output_dir/"active_prime_roots.json").write_text(json.dumps(records,separators=(",",":"))+"\n",encoding="utf-8")
    with (output_dir/"active_prime_roots.csv").open("w",newline="",encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["p","root_count","first_root","middle_is_root","roots"])
        for r in records:
            writer.writerow([r["p"],r["root_count"],r["first_root"],r["middle_is_root"]," ".join(map(str,r["roots"]))])
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit",type=int,default=20000)
    parser.add_argument("--output-dir",type=Path,default=Path(__file__).parent)
    args = parser.parse_args()
    run(args.limit,args.output_dir)
