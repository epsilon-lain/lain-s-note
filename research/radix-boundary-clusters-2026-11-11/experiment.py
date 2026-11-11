#!/usr/bin/env python3
"""Certified variable-radix boundary clusters for pi; standard library only.

This classifies finite ranges, not all bases. Legendre's theorem supplies
completeness of the convergent list; integer intervals certify all decisions.
"""
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
import json

BITS = 11200
SCALE = 1 << BITS
LIMIT = 10 ** 1000
DIRECT_LIMIT = 100000
THRESHOLDS = (10, 100, 1000)


def atan_interval(inv):
    """Outward-rounded dyadic alternating-series interval for atan(1/inv)."""
    lo = hi = k = 0
    power = inv
    while True:
        denominator = (2*k+1)*power
        if denominator.bit_length() > BITS + 40:
            break
        numerator = SCALE if k % 2 == 0 else -SCALE
        lo += numerator // denominator
        hi += -((-numerator) // denominator)
        power *= inv*inv
        k += 1
    numerator = SCALE if k % 2 == 0 else -SCALE
    lo += min(0, numerator // denominator)
    hi += max(0, -((-numerator) // denominator))
    return lo, hi, k


def pi_interval():
    l5, u5, n5 = atan_interval(5)
    l239, u239, n239 = atan_interval(239)
    return 16*l5-4*u239, 16*u5-4*l239, [n5, n239]


def floor_interval(lo, hi, denominator=SCALE):
    a, b = lo // denominator, hi // denominator
    assert a == b, "Increase precision to certify integer part"
    return a


def decimal_interval(lo, hi, places=30):
    base = 10 ** places
    a = lo*base // SCALE
    b = -((-hi*base) // SCALE)
    def display(n):
        return str(n//base) + "." + str(n % base).zfill(places)
    return [display(a), display(b)]


def quality(p, q, lo, hi):
    assert gcd(p, q) == 1
    a, b = q*lo-p*SCALE, q*hi-p*SCALE
    if a > 0:
        return q*a, q*b, "below"
    assert b < 0, "Increase precision to separate pi from the fraction"
    return -q*b, -q*a, "above"


def convergents(lo, hi):
    lower, upper = Fraction(lo, SCALE), Fraction(hi, SCALE)
    pm2, pm1, qm2, qm1 = 0, 1, 1, 0
    digits, rows = [], []
    while True:
        a = lower.numerator // lower.denominator
        assert a == upper.numerator // upper.denominator
        digits.append(a)
        p, q = a*pm1+pm2, a*qm1+qm2
        if q > LIMIT:
            return digits, rows, (p, q)
        rows.append((p, q))
        pm2, pm1, qm2, qm1 = pm1, p, qm1, q
        lower, upper = 1/(upper-a), 1/(lower-a)


def radix_digits(base, lo, hi):
    power, previous = 1, 3
    digits = []
    for _ in range(3):
        power *= base
        part = floor_interval(power*lo, power*hi)
        digits.append(part-base*previous)
        previous = part
    assert all(0 <= d < base for d in digits)
    return digits


def digital_index(base, digits):
    _, d2, d3 = digits
    if d2 == 0:
        return d3
    if d2 == base-1:
        return base-1-d3
    return base


def main():
    lo, hi, terms = pi_interval()
    assert 3*SCALE < lo < hi < 4*SCALE
    digits, rows, next_row = convergents(lo, hi)
    c_hits = {str(t): set() for t in THRESHOLDS}
    delta_hits = {str(t): set() for t in THRESHOLDS}
    zero = []
    families = []
    stream = sha256()
    count_envelope = 0
    min_positive = None
    delta_maxima = {str(t): {"maximum_untruncated_multiplier":0,
                             "maximum_B_capped_multiplier":0,
                             "untruncated_maximizers":[], "B_capped_maximizers":[]}
                    for t in THRESHOLDS}
    for index, (p, q) in enumerate(rows):
        l, u, side = quality(p, q, lo, hi)
        assert l > 0
        multiplier_limits = {}
        for t in THRESHOLDS:
            maximum = isqrt((SCALE-1)//(t*l))
            assert t*maximum*maximum*u < SCALE
            assert t*(maximum+1)*(maximum+1)*l > SCALE
            capped = min(maximum, LIMIT//q)
            multiplier_limits[str(t)] = {"untruncated_maximum":maximum,
                                         "B_capped_maximum":capped,
                                         "minimum_base_four_multiplier":max(1,(4+q-1)//q),
                                         "members_in_base_range":max(0,capped-max(1,(4+q-1)//q)+1)}
            witness = {"convergent_index":index,"p":str(p),"q":str(q),
                       "primitive_quality_decimal_interval":decimal_interval(l,u),
                       "untruncated_maximum":maximum,"B_capped_maximum":capped}
            for value, name, witness_key in ((maximum,"maximum_untruncated_multiplier","untruncated_maximizers"),
                                             (capped,"maximum_B_capped_multiplier","B_capped_maximizers")):
                entry = delta_maxima[str(t)]
                if value > entry[name]:
                    entry[name] = value
                    entry[witness_key] = [witness]
                elif value == entry[name] and value > 0:
                    entry[witness_key].append(witness)
        # Every C(b)<1/10, b>=4 has Delta(b)<1/10+1/4=7/20.
        # k^2*l/SCALE<7/20 is a safe superset of the true possibilities.
        kmax = min(LIMIT//q, isqrt((7*SCALE-1)//(20*l)))
        members = []
        for k in range(max(1, (4+q-1)//q), kmax+1):
            dl, du, b = k*k*l, k*k*u, k*q
            if 20*dl >= 7*SCALE:
                continue
            assert 20*du < 7*SCALE, "Unresolved envelope endpoint"
            assert gcd(k*p, b) == k
            raw_nearest = floor_interval(2*b*lo+SCALE, 2*b*hi+SCALE, 2*SCALE)
            assert raw_nearest == k*p
            a = floor_interval(b*dl, b*du)
            ds = radix_digits(b, lo, hi)
            assert digital_index(b, ds) == a
            assert ds[1] == (0 if side == "below" else b-1)
            count_envelope += 1
            stream.update(f"{index},{k},{b},{a},{ds[0]},{ds[1]},{ds[2]}\n".encode())
            qualifying_c, qualifying_delta = [], []
            for t in THRESHOLDS:
                if t*a < b:
                    c_hits[str(t)].add(b)
                    qualifying_c.append(t)
                if t*du < SCALE:
                    delta_hits[str(t)].add(b)
                    qualifying_delta.append(t)
                else:
                    assert t*dl > SCALE, "Unresolved quality threshold"
            if a == 0:
                zero.append(b)
            elif min_positive is None or a*min_positive[1] < min_positive[0]*b:
                min_positive = (a, b, p, q, k)
            if qualifying_c or qualifying_delta or a == 0:
                members.append({"multiplier":k, "gcd_of_nearest_numerator_and_base":k,
                                "first_three_digits":[str(d) for d in ds],
                                "C_numerator_over_base":str(a),
                                "C_below_reciprocals":qualifying_c,
                                "Delta_below_reciprocals":qualifying_delta})
        if members:
            families.append({"convergent_index":index,"p":str(p),"q":str(q),
                             "side_of_fraction":side,
                             "primitive_quality_decimal_interval":decimal_interval(l,u),
                             "Delta_multiplier_limits":multiplier_limits,
                             "members":members})

    assert min_positive and 10*min_positive[0] < min_positive[1]
    # Thus this is the global positive minimum of C on the finite base range:
    # every smaller positive value would belong to the complete C<1/10 set.
    direct_hash = sha256()
    direct_c = {str(t): set() for t in THRESHOLDS}
    direct_delta = {str(t): set() for t in THRESHOLDS}
    direct_zero = []
    boundary_count = 0
    for b in range(4, DIRECT_LIMIT+1):
        ds = radix_digits(b, lo, hi)
        a = digital_index(b, ds)
        raw_nearest = floor_interval(2*b*lo+SCALE, 2*b*hi+SCALE, 2*SCALE)
        signed_lo, signed_hi = b*lo-raw_nearest*SCALE, b*hi-raw_nearest*SCALE
        if signed_lo > 0:
            dl, du = b*signed_lo, b*signed_hi
        else:
            assert signed_hi < 0
            dl, du = -b*signed_hi, -b*signed_lo
        floor_scaled_delta = floor_interval(b*dl,b*du)
        assert a == min(b, floor_scaled_delta)
        if a < b:
            boundary_count += 1
        if a == 0:
            direct_zero.append(b)
        for t in THRESHOLDS:
            if t*a < b:
                direct_c[str(t)].add(b)
            if t*du < SCALE:
                direct_delta[str(t)].add(b)
            else:
                assert t*dl > SCALE
        direct_hash.update(f"{b},{raw_nearest},{','.join(map(str,ds))},{a}\n".encode())
    for t in THRESHOLDS:
        assert direct_c[str(t)] == {b for b in c_hits[str(t)] if b <= DIRECT_LIMIT}
        assert direct_delta[str(t)] == {b for b in delta_hits[str(t)] if b <= DIRECT_LIMIT}
    assert direct_zero == sorted(b for b in zero if b <= DIRECT_LIMIT)

    def summarize_sets(sets):
        result = {}
        for t, bases in sets.items():
            ordered = sorted(bases)
            count_centers = sum(any(int(t) in m[("C" if sets is c_hits else "Delta")+"_below_reciprocals"]
                                    for m in f["members"]) for f in families)
            result[t] = {"count":len(ordered),"independent_reduced_centers":count_centers,
                         "first_base":str(ordered[0]) if ordered else None,
                         "last_base":str(ordered[-1]) if ordered else None,
                         "sorted_bases_sha256":sha256(("\n".join(map(str,ordered))+"\n").encode()).hexdigest(),
                         "bases_up_to_direct_limit":[str(b) for b in ordered if b <= DIRECT_LIMIT]}
        return result

    a, b, p, q, k = min_positive
    result = {"status":"PASS", "scope":"Complete finite range 4<=integer base b<=10^1000; not an infinite theorem",
              "definitions":{"Delta":"b*||b*pi||", "C":"min(1,floor(b*Delta)/b)",
                             "member_base":"multiplier*q", "member_nearest_numerator":"multiplier*p",
                             "C_value":"C_numerator_over_base/member_base",
                             "threshold_key":"t means strict threshold 1/t"},
              "pi_interval":{"method":"16 atan(1/5)-4 atan(1/239), directed dyadic rounding and alternating remainder",
                             "bits":BITS,"terms":terms,"lower_numerator":str(lo),"upper_numerator":str(hi),
                             "denominator_power_of_two":BITS,"width_numerator":str(hi-lo)},
              "base_limit":str(LIMIT), "certified_cf_digits":digits,
              "convergents_within_limit":len(rows),
              "first_convergent_above_limit":{"p":str(next_row[0]),"q":str(next_row[1])},
              "safe_envelope":"Delta<7/20", "envelope_candidate_count":count_envelope,
              "envelope_candidate_sha256":stream.hexdigest(),
              "C_thresholds":summarize_sets(c_hits),"Delta_thresholds":summarize_sets(delta_hits),
              "Delta_cluster_maxima_over_centers_with_q_at_most_B":delta_maxima,
              "all_C_zero_bases":[str(b) for b in sorted(zero)],
              "global_minimum_positive_C":{"numerator":str(a),"denominator":str(b),
                                           "primitive_p":str(p),"primitive_q":str(q),"multiplier":k},
              "families":families,
              "direct_crosscheck":{"base_range":[4,DIRECT_LIMIT],"count":DIRECT_LIMIT-3,
                                   "boundary_second_digit_count":boundary_count,
                                   "candidate_stream_sha256":direct_hash.hexdigest(),
                                   "all_three_threshold_sets_and_zeros_match":True},
              "script_sha256":sha256(Path(__file__).read_bytes()).hexdigest()}
    target = Path(__file__).with_name("results.json")
    target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":"PASS","convergents":len(rows),"envelope_candidates":count_envelope,
                      "families":len(families),"C_counts":{t:len(v) for t,v in c_hits.items()},
                      "Delta_counts":{t:len(v) for t,v in delta_hits.items()},
                      "C_zero_bases":sorted(zero),"direct_bases":DIRECT_LIMIT-3,
                      "bytes":target.stat().st_size},ensure_ascii=False))


if __name__ == "__main__":
    main()
