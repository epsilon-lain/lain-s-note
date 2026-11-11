#!/usr/bin/env python3
"""Independent exact audit of observation-gcd properties.

Uses a different identity for pi and never imports the generator. Older
continued-fraction windows are reused only after their archived hashes and
independent verification have been checked.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path
import json

BITS = 2304
SCALE = 1 << BITS
DIRECT_LIMIT = 100000
HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "finite-observation-hypothesis-2026-11-11"


def atan_bounds(inv):
    # A deliberately separate implementation: stop on an exact remainder
    # bound, direct each term outwards, then append the alternating remainder.
    low = high = 0
    power = inv
    index = 0
    while True:
        den = (2*index+1)*power
        signed = SCALE if index % 2 == 0 else -SCALE
        if den > (SCALE << 32):
            low += min(0, signed // den)
            high += max(0, -((-signed) // den))
            break
        low += signed // den
        high += -((-signed) // den)
        index += 1
        power *= inv*inv
    return low, high, index


def independent_pi():
    a, b, terms2 = atan_bounds(2)
    c, d, terms3 = atan_bounds(3)
    return 4*(a+c), 4*(b+d), [terms2, terms3]


def integer_part(low, high, denominator=SCALE):
    a, b = low // denominator, high // denominator
    assert a == b, "Independent precision failed to decide an integer"
    return a


def readout(B, low, high):
    return integer_part(2*B*low+SCALE, 2*B*high+SCALE, 2*SCALE)


def source_certificate():
    manifest = json.loads((SOURCE / "manifest.json").read_text())
    hashes = {}
    for name in ("results.json", "verification.json", "experiment.py", "verify.py"):
        data = (SOURCE / name).read_bytes()
        hashes[name] = sha256(data).hexdigest()
        assert hashes[name] == manifest["files"][name]["sha256"]
    original = json.loads((SOURCE / "results.json").read_text())
    audit = json.loads((SOURCE / "verification.json").read_text())
    assert original["status"] == audit["status"] == "PASS"
    assert audit["inputs_sha256"]["results.json"] == hashes["results.json"]
    assert audit["inputs_sha256"]["experiment.py"] == hashes["experiment.py"]
    assert audit["continued_fractions"]["centers_verified"] == len(original["centers"])
    centers = [(int(c["p"]), int(c["q"]), int(c["full_rounding_multiplier_K"]))
               for c in original["centers"]]
    return original, centers, hashes


def exact_count(H, threshold, centers):
    total = 0
    for _, q, last in centers:
        # q^2 < H / threshold is necessary. Not using this pruning keeps
        # the independent implementation different and makes no difference
        # to the exact sum.
        upper = min(last, H // q)
        first = threshold*q + 1
        if upper >= first:
            total += upper-first+1
    return total


def direct_scan(low, high):
    histogram = Counter()
    divisible = {d:0 for d in range(1, 21)}
    residue = {str(d):[0]*4 for d in range(1, 6)}
    residue["other"] = [0]*4
    events = {t:[] for t in (1, 10, 100, 1000)}
    stream = sha256()
    adjacent_stream = sha256()
    pair_checks = 0
    g_sum = 0
    reciprocal_lower = reciprocal_upper = 0
    reciprocal_bits = 128
    previous = None
    for B in range(1, DIRECT_LIMIT+1):
        A = readout(B, low, high)
        g = gcd(B, A)
        m = A-3*B
        histogram[g] += 1
        g_sum += g
        for d in divisible:
            divisible[d] += (g % d == 0)
        lower_residual = B*low-A*SCALE
        upper_residual = B*high-A*SCALE
        quartile = integer_part(4*lower_residual+2*SCALE,
                                4*upper_residual+2*SCALE)
        assert 0 <= quartile < 4
        residue[str(g) if g <= 5 else "other"][quartile] += 1
        reciprocal_lower += (1 << reciprocal_bits)//g
        reciprocal_upper += -((-(1 << reciprocal_bits))//g)
        for threshold in events:
            if g*g > threshold*B:
                events[threshold].append(B)
        stream.update(f"{B},{A},{g},{quartile}\n".encode())
        if previous is not None:
            old_B, old_A, old_g, old_m = previous
            if old_B >= 4:
                carry = m-old_m
                assert carry in (0, 1)
                determinant = old_m if carry == 0 else old_B-old_m
                assert determinant > 0
                assert determinant % (old_g*g) == 0
                assert (old_g*g)**2 <= determinant**2 < old_B*B
                assert gcd(old_g, g) == 1
                adjacent_stream.update(f"{old_B},{old_g},{g},{old_m},{carry},{determinant}\n".encode())
                pair_checks += 1
        previous = (B, A, g, m)
    return {"full_histogram":histogram, "divisibility":divisible,
            "residual":residue, "events":events, "g_sum":g_sum,
            "readout_sha256":stream.hexdigest(), "adjacent_checks":pair_checks,
            "adjacent_sha256":adjacent_stream.hexdigest(),
            "reciprocal_bits":reciprocal_bits,
            "reciprocal_lower":reciprocal_lower,
            "reciprocal_upper":reciprocal_upper}


def main():
    generated_bytes = (HERE / "results.json").read_bytes()
    generated = json.loads(generated_bytes)
    assert generated["status"] == "PASS"
    assert generated["script_sha256"] == sha256((HERE / "experiment.py").read_bytes()).hexdigest()
    original, centers, source_hashes = source_certificate()
    reused = generated["reused_certificate"]
    assert reused["results_sha256"] == source_hashes["results.json"]
    assert reused["manifest_sha256"] == sha256((SOURCE / "manifest.json").read_bytes()).hexdigest()
    assert reused["center_count"] == len(centers) == 969
    assert int(reused["N_max"]) == int(original["resolution_limit"]) == 10**1000

    lo, hi, terms = independent_pi()
    own_den = 1 << generated["pi_interval"]["bits"]
    old_lo = int(generated["pi_interval"]["lower_numerator"])
    old_hi = int(generated["pi_interval"]["upper_numerator"])
    assert old_lo*SCALE < lo*own_den < hi*own_den < old_hi*SCALE
    scanned = direct_scan(lo, hi)
    snap = generated["direct_scan"]["snapshots"][str(DIRECT_LIMIT)]
    assert snap["count"] == DIRECT_LIMIT
    hist = scanned["full_histogram"]
    assert snap["g_histogram_1_to_20"] == {str(d):hist[d] for d in range(1, 21)}
    assert snap["g_above_20_count"] == sum(c for g,c in hist.items() if g > 20)
    assert snap["divisibility_counts_1_to_20"] == {str(d):scanned["divisibility"][d] for d in range(1, 21)}
    assert snap["residual_quartiles_by_g"] == scanned["residual"]
    assert snap["g_sum"] == str(scanned["g_sum"])
    assert snap["readout_stream_sha256"] == scanned["readout_sha256"]
    assert snap["adjacent_stream_sha256"] == scanned["adjacent_sha256"]
    assert snap["adjacent_checks"] == scanned["adjacent_checks"]
    mean = snap["mean_reciprocal_g_dyadic_bound"]
    assert mean["bits"] == scanned["reciprocal_bits"]
    assert int(mean["lower_sum_numerator"]) == scanned["reciprocal_lower"]
    assert int(mean["upper_sum_numerator"]) == scanned["reciprocal_lower"]+DIRECT_LIMIT
    assert int(mean["denominator"]) == DIRECT_LIMIT << scanned["reciprocal_bits"]
    exact_mean = sum((F(c,g) for g,c in hist.items()), F(0))/DIRECT_LIMIT
    assert F(mean["decimal_interval"][0]) <= exact_mean <= F(mean["decimal_interval"][1])

    for threshold, actual in scanned["events"].items():
        claimed = snap["high_G_by_threshold"][str(threshold)]
        assert claimed["count"] == len(actual) == exact_count(DIRECT_LIMIT, threshold, centers)
        digest = sha256("".join(f"{B}\n" for B in actual).encode()).hexdigest()
        assert claimed["sha256"] == digest
        if "bases" in claimed:
            assert claimed["bases"] == actual
        else:
            assert claimed["bases_first_20"] == actual[:20]
            assert claimed["bases_last_20"] == actual[-20:]
            expanded = []
            for segment in claimed["arithmetic_segments"]:
                index = segment["center_index"]
                p, q, K = centers[index]
                assert int(segment["q"]) == q
                first = int(segment["k_first"])
                last = int(segment["k_last"])
                assert first == threshold*q+1
                assert last == min(K, DIRECT_LIMIT//q)
                expanded.extend(k*q for k in range(first,last+1))
            assert sorted(expanded) == actual

    count_checks = 0
    for row in generated["large_counts"]:
        H = int(row["N"])
        assert H == 10**row["decimal_exponent"] <= 10**1000
        for key, claimed in row["counts"].items():
            threshold = int(key)
            actual = exact_count(H, threshold, centers)
            assert int(claimed) == actual
            numerator, denominator = row["squared_normalized_counts"][key]
            assert int(numerator) == actual*actual and int(denominator) == H
            count_checks += 1

    # Independently reconstruct all published cross-center checkpoints.
    checkpoints = []
    for index,(p,q,K) in enumerate(centers):
        upper = min(K,10**1000//q)
        for k in {q+1,upper}:
            if q < k <= upper:
                checkpoints.append((k*q,k,index,p,q))
    checkpoints.sort()
    pair_hash = sha256()
    checks = separated = 0
    for i in range(len(checkpoints)):
        for gap in (1,7,31):
            j = i+gap
            if j >= len(checkpoints):
                continue
            B,g,index,p,q = checkpoints[i]
            C,h,other,r,s = checkpoints[j]
            if index == other:
                continue
            assert B < C
            determinant = abs(g*p*C-h*r*B)
            # Center determinant times g*h is an independent exact check.
            assert determinant == g*h*abs(p*s-r*q)
            assert determinant > 0 and 2*determinant < B+C
            assert 4*g*g*h*h < (B+C)**2
            for threshold in (1,10,100,1000):
                if g*g > threshold*B and h*h > threshold*C:
                    assert C*C-(4*threshold*threshold-2)*B*C+B*B > 0
                    separated += 1
            pair_hash.update(f"{i},{j},{B},{C},{g},{h},{index},{other}\n".encode())
            checks += 1
    pairs = generated["different_center_pair_checks"]
    assert pairs["checkpoint_count"] == len(checkpoints)
    assert pairs["pair_count"] == checks
    assert pairs["threshold_separation_checks"] == separated
    assert pairs["sha256"] == pair_hash.hexdigest()

    for item in generated["neighboring_peak_example"]:
        B = item["B"]
        A = readout(B,lo,hi)
        g = gcd(A,B)
        assert item["A"] == A and item["gcd"] == g
        assert (item["reduced_p"],item["reduced_q"]) == (A//g,B//g)
        assert item["G_ratio"] == [str(g*g),str(B)]
        assert F(item["G_decimal_interval"][0]) <= F(g*g,B) <= F(item["G_decimal_interval"][1])

    result = {
        "status":"PASS",
        "scope":"Independent finite property audit; no infinite pi bound",
        "independent_pi":{"method":"4*(atan(1/2)+atan(1/3)) with directed integer bounds",
                          "bits":BITS,"terms":terms,"lower_numerator":str(lo),
                          "upper_numerator":str(hi),"strictly_inside_generator_interval":True},
        "reused_source":{"directory":str(SOURCE.name),"verified_hashes":source_hashes,
                         "independent_source_audit_passes":True,"centers":len(centers)},
        "direct_scan":{"resolutions":[1,DIRECT_LIMIT],
                       "all_histograms_and_divisibility_frequencies_match":True,
                       "all_residual_quartile_counts_match":True,
                       "readout_stream_sha256":scanned["readout_sha256"],
                       "adjacent_stream_sha256":scanned["adjacent_sha256"],
                       "adjacent_checks":scanned["adjacent_checks"],
                       "all_reciprocal_mean_bounds_match":True,
                       "high_counts":{str(t):len(v) for t,v in scanned["events"].items()}},
        "large_count_checks":count_checks,
        "different_center_pairs_verified":checks,
        "threshold_separations_verified":separated,
        "neighboring_peak_examples_verified":len(generated["neighboring_peak_example"]),
        "inputs_sha256":{"experiment.py":sha256((HERE/"experiment.py").read_bytes()).hexdigest(),
                         "results.json":sha256(generated_bytes).hexdigest()},
        "script_sha256":sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/"verification.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"status":"PASS","direct":DIRECT_LIMIT,"large_count_checks":count_checks,
                      "pairs":checks,"high_counts":result["direct_scan"]["high_counts"]}))


if __name__ == "__main__":
    main()
