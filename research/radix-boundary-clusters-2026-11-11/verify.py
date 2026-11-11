#!/usr/bin/env python3
"""Independent exact audit of radix boundary clusters.

No imports from experiment.py.  Pi is reconstructed from
4*(atan(1/2)+atan(1/3)); ordinary Euclidean algorithms on the two
rational endpoints independently certify the continued fraction.
"""
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
import json
import sys

sys.set_int_max_str_digits(0)
BITS = 11520
SCALE = 1 << BITS
LIMIT = 10**1000
DIRECT_LIMIT = 10000
THRESHOLDS = (10, 100, 1000)


def atan_bounds(m):
    # Each term is independently rounded outward; the alternating-series
    # tail has the sign of the first omitted term and smaller magnitude.
    low = high = 0
    denominator_power = m
    n = 0
    while True:
        denominator = (2*n+1)*denominator_power
        if denominator > SCALE * (1 << 48):
            if n % 2:
                low -= (SCALE+denominator-1)//denominator
            else:
                high += (SCALE+denominator-1)//denominator
            return low, high, n
        quotient, remainder = divmod(SCALE, denominator)
        ceiling = quotient + bool(remainder)
        if n % 2:
            low -= ceiling
            high -= quotient
        else:
            low += quotient
            high += ceiling
        denominator_power *= m*m
        n += 1


def euclidean_cf(numerator, denominator):
    digits = []
    while denominator:
        quotient, remainder = divmod(numerator, denominator)
        digits.append(quotient)
        numerator, denominator = denominator, remainder
    return digits


def same_floor(lower, upper, denominator=SCALE):
    assert lower//denominator == upper//denominator, "Insufficient audit precision"
    return lower//denominator


def quality_bounds(p, q, lower, upper):
    a, b = q*lower-p*SCALE, q*upper-p*SCALE
    assert a > 0 or b < 0
    if a > 0:
        return q*a, q*b, "below"
    return -q*b, -q*a, "above"


def three_digits(base, lower, upper):
    # Deliberately use iterative remainder propagation, rather than the
    # generator's floors of b*pi, b^2*pi and b^3*pi.
    rlow, rhigh = lower-3*SCALE, upper-3*SCALE
    answer = []
    for _ in range(3):
        rlow *= base
        rhigh *= base
        d = same_floor(rlow, rhigh)
        assert 0 <= d < base
        answer.append(d)
        rlow -= d*SCALE
        rhigh -= d*SCALE
    return answer


def digit_numerator(base, digits):
    if digits[1] == 0:
        return digits[2]
    if digits[1] == base-1:
        return base-1-digits[2]
    return base


def digest_bases(bases):
    return sha256(("\n".join(map(str,sorted(bases)))+"\n").encode()).hexdigest()


def main():
    folder = Path(__file__).parent
    source = json.loads((folder/"results.json").read_text())
    assert source["status"] == "PASS"
    assert source["base_limit"] == str(LIMIT)
    assert source["script_sha256"] == sha256((folder/"experiment.py").read_bytes()).hexdigest()
    l2, h2, n2 = atan_bounds(2)
    l3, h3, n3 = atan_bounds(3)
    lower, upper = 4*(l2+l3), 4*(h2+h3)
    assert 3*SCALE < lower < upper < 4*SCALE
    original = source["pi_interval"]
    scale_original = 1 << original["bits"]
    assert int(original["lower_numerator"])*SCALE < lower*scale_original
    assert upper*scale_original < int(original["upper_numerator"])*SCALE

    left = euclidean_cf(lower,SCALE)
    right = euclidean_cf(upper,SCALE)
    common = []
    for a,b in zip(left,right):
        if a != b:
            break
        common.append(a)
    needed = source["certified_cf_digits"]
    assert common[:len(needed)] == needed
    assert len(common) >= len(needed)
    rows = []
    pm2,pm1,qm2,qm1 = 0,1,1,0
    for index,a in enumerate(common):
        p,q = a*pm1+pm2,a*qm1+qm2
        assert gcd(p,q) == 1
        if q > LIMIT:
            assert {"p":str(p),"q":str(q)} == source["first_convergent_above_limit"]
            assert index+1 == len(needed)
            break
        rows.append((p,q))
        pm2,pm1,qm2,qm1 = pm1,p,qm1,q
    else:
        raise AssertionError("Common continued-fraction prefix too short")
    assert len(rows) == source["convergents_within_limit"]

    # Use the larger 1/2 envelope and then independently reject points
    # outside 7/20; this is distinct from the generator's 7/20 bound.
    c_sets = {t:set() for t in THRESHOLDS}
    delta_sets = {t:set() for t in THRESHOLDS}
    zero = []
    envelope_hash = sha256()
    envelope_count = 0
    family_by_index = {f["convergent_index"]:f for f in source["families"]}
    observed_families = set()
    observed_members = 0
    best_positive = None
    cluster_lift_checks = 0
    raw_centers = {}
    for index,(p,q) in enumerate(rows):
        klow,khigh,side = quality_bounds(p,q,lower,upper)
        kmax = min(LIMIT//q,isqrt((SCALE-1)//(2*klow)))
        family = family_by_index.get(index)
        if family:
            assert (family["p"],family["q"],family["side_of_fraction"]) == (str(p),str(q),side)
            display_low,display_high = map(Fraction,family["primitive_quality_decimal_interval"])
            assert display_low < Fraction(klow,SCALE)
            assert Fraction(khigh,SCALE) < display_high
            stored = {m["multiplier"]:m for m in family["members"]}
        else:
            stored = {}
        for k in range(1,kmax+1):
            b = k*q
            if b < 4:
                continue
            dl,du = k*k*klow,k*k*khigh
            if 20*dl > 7*SCALE:
                continue
            assert 20*du < 7*SCALE
            nearest = same_floor(2*b*lower+SCALE,2*b*upper+SCALE,2*SCALE)
            assert nearest == k*p and gcd(nearest,b) == k
            assert b not in raw_centers
            raw_centers[b] = (p,q)
            ds = three_digits(b,lower,upper)
            a = digit_numerator(b,ds)
            assert a == same_floor(b*dl,b*du)
            assert a == same_floor(k**3*q*klow,k**3*q*khigh)
            if q >= 4:
                base_digits = three_digits(q,lower,upper)
                a0 = digit_numerator(q,base_digits)
                assert a0 == same_floor(q*klow,q*khigh)
                fractional_low = q*klow-a0*SCALE
                fractional_high = q*khigh-a0*SCALE
                assert a == k**3*a0+same_floor(k**3*fractional_low,k**3*fractional_high)
                cluster_lift_checks += 1
            envelope_count += 1
            envelope_hash.update(f"{index},{k},{b},{a},{ds[0]},{ds[1]},{ds[2]}\n".encode())
            c_hit = [t for t in THRESHOLDS if t*a < b]
            d_hit = []
            for t in THRESHOLDS:
                if t*du < SCALE:
                    d_hit.append(t)
                else:
                    assert t*dl > SCALE
            for t in c_hit:
                c_sets[t].add(b)
            for t in d_hit:
                delta_sets[t].add(b)
            if a == 0:
                zero.append(b)
            elif best_positive is None or Fraction(a,b) < best_positive[0]:
                best_positive = (Fraction(a,b),a,b,p,q,k)
            if c_hit or d_hit or a == 0:
                m = stored[k]
                assert m == {"multiplier":k,"gcd_of_nearest_numerator_and_base":k,
                             "first_three_digits":[str(z) for z in ds],
                             "C_numerator_over_base":str(a),"C_below_reciprocals":c_hit,
                             "Delta_below_reciprocals":d_hit}
                observed_families.add(index)
                observed_members += 1
            else:
                assert k not in stored
    assert observed_families == set(family_by_index)
    assert observed_members == sum(len(f["members"]) for f in source["families"])
    # Independently certify every untruncated cluster length, not only
    # members that happen to fit below the finite base cutoff.
    length_checks = 0
    maximum_summary = {}
    for t in THRESHOLDS:
        measured = []
        for index,(p,q) in enumerate(rows):
            klow,khigh,_ = quality_bounds(p,q,lower,upper)
            maximum = isqrt((SCALE-1)//(t*klow))
            assert t*maximum*maximum*khigh < SCALE
            assert t*(maximum+1)**2*klow > SCALE
            capped = min(maximum,LIMIT//q)
            measured.append((maximum,capped,index,p,q))
            if index in family_by_index:
                first = max(1,(4+q-1)//q)
                assert family_by_index[index]["Delta_multiplier_limits"][str(t)] == {
                    "untruncated_maximum":maximum,"B_capped_maximum":capped,
                    "minimum_base_four_multiplier":first,
                    "members_in_base_range":max(0,capped-first+1)}
            length_checks += 1
        record = source["Delta_cluster_maxima_over_centers_with_q_at_most_B"][str(t)]
        for position,key,witness_key in ((0,"maximum_untruncated_multiplier","untruncated_maximizers"),
                (1,"maximum_B_capped_multiplier","B_capped_maximizers")):
            best = max(row[position] for row in measured)
            witnesses = [row for row in measured if row[position] == best and best > 0]
            assert record[key] == best
            actual = record[witness_key]
            assert len(actual) == len(witnesses)
            for entry,(maximum,capped,index,p,q) in zip(actual,witnesses):
                assert entry["convergent_index"] == index
                assert (entry["p"],entry["q"]) == (str(p),str(q))
                assert entry["untruncated_maximum"] == maximum
                assert entry["B_capped_maximum"] == capped
                klow,khigh,_ = quality_bounds(p,q,lower,upper)
                display_low,display_high = map(Fraction,entry["primitive_quality_decimal_interval"])
                assert display_low < Fraction(klow,SCALE) < Fraction(khigh,SCALE) < display_high
        maximum_summary[str(t)] = {"untruncated":record["maximum_untruncated_multiplier"],
                                  "B_capped":record["maximum_B_capped_multiplier"]}
    assert envelope_count == source["envelope_candidate_count"]
    assert envelope_hash.hexdigest() == source["envelope_candidate_sha256"]
    assert [str(b) for b in sorted(zero)] == source["all_C_zero_bases"]
    for title,sets in (("C_thresholds",c_sets),("Delta_thresholds",delta_sets)):
        for t,bases in sets.items():
            entry = source[title][str(t)]
            assert entry["count"] == len(bases)
            assert entry["sorted_bases_sha256"] == digest_bases(bases)
            assert entry["independent_reduced_centers"] == len({raw_centers[b] for b in bases})
            assert entry["first_base"] == (str(min(bases)) if bases else None)
            assert entry["last_base"] == (str(max(bases)) if bases else None)
    value,a,b,p,q,k = best_positive
    assert value < Fraction(1,10)
    assert source["global_minimum_positive_C"] == {"numerator":str(a),"denominator":str(b),
        "primitive_p":str(p),"primitive_q":str(q),"multiplier":k}

    # Direct digit/remainder audit at every small base, not only candidates.
    direct_c = {t:set() for t in THRESHOLDS}
    direct_d = {t:set() for t in THRESHOLDS}
    direct_zeros = []
    direct_hash = sha256()
    boundary_count = 0
    for b in range(4,DIRECT_LIMIT+1):
        nearest = same_floor(2*b*lower+SCALE,2*b*upper+SCALE,2*SCALE)
        dl,du,_ = quality_bounds(nearest,b,lower,upper)
        ds = three_digits(b,lower,upper)
        a = digit_numerator(b,ds)
        assert a == min(b,same_floor(b*dl,b*du))
        boundary_count += a < b
        if a == 0:
            direct_zeros.append(b)
        for t in THRESHOLDS:
            if t*a < b:
                direct_c[t].add(b)
            if t*du < SCALE:
                direct_d[t].add(b)
            else:
                assert t*dl > SCALE
        direct_hash.update(f"{b},{nearest},{','.join(map(str,ds))},{a}\n".encode())
    for t in THRESHOLDS:
        assert direct_c[t] == {b for b in c_sets[t] if b <= DIRECT_LIMIT}
        assert direct_d[t] == {b for b in delta_sets[t] if b <= DIRECT_LIMIT}
    assert direct_zeros == sorted(b for b in zero if b <= DIRECT_LIMIT)

    # For different reduced centers, test the exact determinant inequality.
    # Adjacent entries suffice to certify every forbidden ratio interval.
    gap_checks = 0
    for t,bases in delta_sets.items():
        ordered = sorted(bases)
        for b,c in zip(ordered,ordered[1:]):
            if raw_centers[b] != raw_centers[c]:
                assert b*b+c*c > t*b*c
                # In particular these strict 1/t thresholds forbid any
                # two distinct centers within a factor of two.
                assert c > 2*b
                gap_checks += 1

    result = {"status":"PASS","scope":"Independent finite certificate audit; no infinite claim for pi",
        "independent_pi":{"method":"4*(atan(1/2)+atan(1/3)) with directed integer term bounds",
             "bits":BITS,"terms":[n2,n3],"lower_numerator":str(lower),"upper_numerator":str(upper),
             "width_numerator":str(upper-lower),"strictly_inside_generator_interval":True},
        "continued_fractions":{"method":"two independent ordinary Euclidean endpoint expansions",
             "common_prefix_digits":len(common),"required_digits":len(needed),"all_required_digits_match":True},
        "envelope_candidates_verified":envelope_count,"envelope_sha256":envelope_hash.hexdigest(),
        "families_verified":len(observed_families),"published_members_verified":observed_members,
        "third_digit_cubic_lift_checks":cluster_lift_checks,"distinct_center_gap_checks":gap_checks,
        "cluster_length_checks":length_checks,"verified_cluster_length_maxima":maximum_summary,
        "C_counts":{str(t):len(v) for t,v in c_sets.items()},
        "Delta_counts":{str(t):len(v) for t,v in delta_sets.items()},"all_C_zero_bases":sorted(zero),
        "minimum_positive_C_exactly_matches":True,
        "direct_scan":{"base_range":[4,DIRECT_LIMIT],"bases":DIRECT_LIMIT-3,
             "boundary_second_digit_count":boundary_count,"sha256":direct_hash.hexdigest(),
             "all_threshold_sets_and_zero_set_match":True},
        "inputs_sha256":{"experiment.py":sha256((folder/"experiment.py").read_bytes()).hexdigest(),
                         "results.json":sha256((folder/"results.json").read_bytes()).hexdigest()},
        "script_sha256":sha256(Path(__file__).read_bytes()).hexdigest()}
    (folder/"verification.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"status":"PASS","envelope_candidates":envelope_count,
        "families":len(observed_families),"direct_bases":DIRECT_LIMIT-3,
        "C_zero_bases":sorted(zero),"independent_common_cf_digits":len(common)}))


if __name__ == "__main__":
    main()
