"""Certify both shell objectives for alpha=1/pi using outward-rounded integers.

Optionally pass the supplied windows_semiconvergent.csv to compare the winners.
"""
import argparse
import csv
from bisect import bisect_right
from hashlib import sha256

SCALE = 10**3400
LAST_SHELL = 4992


def atan_bounds(inverse):
    lo = hi = 0
    power = inverse
    k = 0
    sign = 1
    while True:
        term = SCALE // ((2*k + 1) * power)
        if term == 0:
            return lo-1, hi+1
        lo += sign * term - int(sign < 0)
        hi += sign * term + int(sign > 0)
        power *= inverse**2
        k += 1
        sign = -sign


def reciprocal_pi():
    a, b = atan_bounds(5)
    c, d = atan_bounds(239)
    return SCALE**2 // (16*b-4*c), SCALE**2 // (16*a-4*d) + 1


def common_convergents(lo, hi):
    lower_n, lower_d = lo, SCALE
    upper_n, upper_d = hi, SCALE
    prev_p, p, prev_r, r = 0, 1, 1, 0
    records = [(1, 0)]
    digits = []
    while lower_d and upper_d:
        a, lower_rem = divmod(lower_n, lower_d)
        b, upper_rem = divmod(upper_n, upper_d)
        if a != b:
            break
        digits.append(a)
        prev_p, p = p, a*p+prev_p
        prev_r, r = r, a*r+prev_r
        records.append((p, r))
        lower_n, lower_d = lower_d, lower_rem
        upper_n, upper_d = upper_d, upper_rem
    return records, digits


def det_one_candidates(records, denominators, q):
    start = bisect_right(denominators, q) - 1
    result = set()
    for n in range(start, len(records)-1):
        p, r = records[n]
        if r >= 2*q:
            break
        if r >= q:
            result.add((p, r))
        pp, rr = records[n-1]
        first = max(r, q)
        last = min(records[n+1][1], 2*q)-1
        for sign in (-1, 1):
            left = -((sign*rr-first)//r)
            right = (last-sign*rr)//r
            if left <= right:
                for k in (left, right):
                    result.add((k*p+sign*pp, k*r+sign*rr))
    assert result
    return result


def farey_candidates(records, denominators, q):
    B = 2*q-1
    n = bisect_right(denominators, B) - 1
    p, r = records[n]
    pp, rr = records[n-1]
    k = (B-rr)//r
    a, b = pp+k*p, rr+k*r
    result = set()
    for (u,v),(x,y) in (((p,r),(a,b)),((a,b),(p,r))):
        if v < q:
            h = (B+y)//v
            u, v = h*u-x, h*v-y
        assert q <= v <= B
        result.add((u,v))
    return result


def score_interval(a, b, q, lo, hi, objective):
    left = b*lo - a*SCALE
    right = b*hi - a*SCALE
    assert left*right > 0
    small, large = sorted((abs(left), abs(right)))
    if objective == "M":
        return small*b, large*b
    return (q*q*small)//b, -((-q*q*large)//b)


def winner(candidates, q, lo, hi, objective):
    scores = [(score_interval(a,b,q,lo,hi,objective),a,b) for a,b in candidates]
    scores.sort()
    if len(scores) > 1:
        assert scores[0][0][1] < scores[1][0][0]
    return scores[0]


def decimal_outer(bounds, places=10):
    a,b = bounds
    den = 10**places
    x,y = a*den//SCALE, -((-b*den)//SCALE)
    return tuple(f"{v//den}.{v%den:0{places}d}" for v in (x,y))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--compare", help="m5/results/windows_semiconvergent.csv")
    args = ap.parse_args()
    lo, hi = reciprocal_pi()
    records, digits = common_convergents(lo, hi)
    assert len(digits) >= 3003 and digits[432] == 20776
    denominators = [r for _, r in records]
    assert denominators[-1] >= 2**(LAST_SHELL+1)
    original = {}
    if args.compare:
        with open(args.compare) as f:
            for row in csv.DictReader(f):
                original[int(row["j"])] = int(row["a"]), int(row["b"]), int(row["semi_min_certified"])
    winners = {"M": [], "D": []}
    corrected = {}
    differing = 0
    certificate = sha256()
    for j in range(LAST_SHELL+1):
        q = 1 << j
        mass = winner(det_one_candidates(records,denominators,q), q,lo,hi,"M")
        dist = winner(farey_candidates(records,denominators,q), q,lo,hi,"D")
        differing += (mass[1:] != dist[1:])
        if j in original and mass[1:] != original[j][:2]:
            corrected[j] = original[j][2]
        if j in (0,730,731):
            print(j,"M",decimal_outer(mass[0]),"D",decimal_outer(dist[0]))
        for key,row in (("M",mass),("D",dist)):
            interval, a,b = row
            winners[key].append((interval, j))
            certificate.update(f"{key} {j} {a} {b}\n".encode("ascii"))
    for key, rows in winners.items():
        (bounds,j), (second,_) = sorted(rows, reverse=True)[:2]
        assert bounds[0]>second[1]
        print("Maximum",key,"j",j,decimal_outer(bounds))
    print("Different D and M winners:",differing)
    print("Sorted by shell then M/D, winners SHA-256:",certificate.hexdigest())
    if original:
        print("Different from input:",len(corrected),"of which labelled certified:",sum(corrected.values()))
        print("Changed shell indices:", sorted(corrected))


if __name__ == "__main__":
    main()
