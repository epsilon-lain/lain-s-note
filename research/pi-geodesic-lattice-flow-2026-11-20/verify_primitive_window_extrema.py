"""Certify both primitive minima in 4,993 dyadic windows for 1/pi.

Python 3 standard library. An optional dyadic_window_approximation.zip argument
compares the reported m5 numerators/denominators with the certified optima.
"""
import argparse
import csv
from bisect import bisect_left
from functools import cmp_to_key
from io import TextIOWrapper
from math import gcd
from zipfile import ZipFile


def scaled_atan(base, scale):
    """Enclose scale*atan(1/base) using outward-rounded alternating terms."""
    lo = hi = 0
    power = base
    index = 0
    while True:
        term = scale // ((2 * index + 1) * power)
        if term == 0:
            return (lo, hi + 1) if index % 2 == 0 else (lo - 1, hi)
        if index % 2 == 0:
            lo, hi = lo + term, hi + term + 1
        else:
            lo, hi = lo - term - 1, hi - term
        index += 1
        power *= base * base


def reciprocal_pi(scale):
    a, b = scaled_atan(5, scale)
    c, d = scaled_atan(239, scale)
    p_lo, p_hi = 16*a - 4*d, 16*b - 4*c
    return scale**2 // p_hi, -(-scale**2 // p_lo)


def convergents(lo, hi, scale, limit):
    lo_n, lo_d = lo, scale
    hi_n, hi_d = hi, scale
    old_a, a, old_r, r = 0, 1, 1, 0
    points = []
    while r < limit:
        digit, lo_rem = divmod(lo_n, lo_d)
        other, hi_rem = divmod(hi_n, hi_d)
        assert digit == other and lo_rem and hi_rem, "increase interval precision"
        old_a, a = a, digit * a + old_a
        old_r, r = r, digit * r + old_r
        points.append((a, r))
        lo_n, lo_d, hi_n, hi_d = hi_d, hi_rem, lo_d, lo_rem
    return points


def candidates(Q, a, r, ap, rp):
    result = set()
    if Q <= r:
        result.add((a, r))
    for m in (-1, 1):
        first = -((m * rp - Q) // r)
        last = (2*Q - 1 - m*rp) // r
        if first <= last:
            for k in {first, last}:
                result.add((k*a + m*ap, k*r + m*rp))
    return result


def compare(x, y):
    value = x[0]*y[1] - x[1]*y[0]
    return (value > 0) - (value < 0)


rational_key = cmp_to_key(compare)


def certified_minimum(ranges):
    # Each entry is a lower bound, an upper bound, and a candidate.
    lower, upper, winner = min(ranges, key=lambda t: rational_key(t[0]))
    assert all(c == winner or compare(upper, other) < 0 for other, _, c in ranges)
    return lower, upper, winner


def minima(Q, points, lo, hi, scale):
    quality, distance = [], []
    for p, b in points:
        assert Q <= b < 2*Q
        u, v = lo*b - scale*p, hi*b - scale*p
        left = 0 if u <= 0 <= v else min(abs(u), abs(v))
        right = max(abs(u), abs(v))
        # Common factor 1/scale is applied only when formatting.
        quality.append(((b*left, 1), (b*right, 1), (p, b)))
        distance.append(((Q*Q*left, b), (Q*Q*right, b), (p, b)))
    return certified_minimum(distance), certified_minimum(quality)


def fixed_digits(value, scale, digits=12):
    n, d = value
    return n * 10**digits // (d*scale)


def optional_table(path):
    with ZipFile(path) as archive, archive.open('m5/results/windows_semiconvergent.csv') as raw:
        return {int(row['j']): (int(row['a']), int(row['b']), int(row['semi_min_certified']))
                for row in csv.DictReader(TextIOWrapper(raw, encoding='utf-8'))}


def ordinary_or_intermediate(point, convergents, denominators):
    p, b = point
    i = bisect_left(denominators, b)
    if denominators[i] == b and convergents[i][0] == p:
        return True
    if i == 0:
        return False
    ap, rp = convergents[i-2] if i >= 2 else (1, 0)
    a, r = convergents[i-1]
    k, remainder = divmod(b-rp, r)
    return remainder == 0 and k > 0 and p == ap + k*a


def export_rows(path, samples, scale):
    with open(path, 'w', newline='', encoding='utf-8') as output:
        writer = csv.writer(output)
        writer.writerow(('j','D_numerator','D_denominator','D_floor_times_10^12',
                         'M_numerator','M_denominator','M_floor_times_10^12'))
        for j, rows in enumerate(samples):
            values = [j]
            for lb, ub, (p, b) in rows:
                value = fixed_digits(lb, scale)
                assert value == fixed_digits(ub, scale)
                values.extend((p, b, value))
            writer.writerow(values)


def run(table, csv_path=None):
    scale = 1 << 10300
    lo, hi = reciprocal_pi(scale)
    points = convergents(lo, hi, scale, 2**4993)
    denominators = [r for _, r in points]
    samples = []
    exceptions = []
    for j in range(4993):
        Q = 1 << j
        i = bisect_left(denominators, 2*Q) - 1
        a, r = points[i]
        ap, rp = points[i-1] if i else (1, 0)
        rows = minima(Q, candidates(Q, a, r, ap, rp), lo, hi, scale)
        samples.append(rows)
        if j < 8:
            exhaustive = {(p, b) for b in range(Q, 2*Q) for p in range(b+1) if gcd(p,b)==1}
            assert rows == minima(Q, exhaustive, lo, hi, scale)
        if j in table and rows[1][2] != table[j][:2]:
            exceptions.append(j)
    for which, label in enumerate(('D','M')):
        j = max(range(len(samples)), key=lambda i: rational_key(samples[i][which][0]))
        low, high, _ = samples[j][which]
        assert all(i==j or compare(low, row[which][1])>0 for i,row in enumerate(samples))
        for index in sorted({j, 730, 731}):
            lb, ub, _ = samples[index][which]
            t = fixed_digits(lb, scale)
            assert t == fixed_digits(ub, scale)
            print(f'{label} j={index}: ({t}, {t+1}) * 10^-12' + (' [maximum]' if index == j else ''))
    same = sum(d[2] == m[2] for d,m in samples)
    print(f'Same minimizer for D and M: {same}/{len(samples)} windows')
    for i, label in enumerate(('D','M')):
        outside = sum(not ordinary_or_intermediate(rows[i][2],points,denominators) for rows in samples)
        print(f'{label} winners outside ordinary and intermediate convergents: {outside}')
    if csv_path:
        export_rows(csv_path, samples, scale)
    if table:
        assert set(table) == set(range(1,4993))
        print('m5 quality candidates differing from exact optimum:',len(exceptions))
        print('Of these, marked semi_min_certified:', sum(table[j][2] for j in exceptions))
        print('First differing indices:', exceptions[:20])
        if exceptions:
            j = exceptions[0]
            print('First corrected winner:',j,samples[j][1][2], 'table:',table[j])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', nargs='?', help='Optional original dyadic ZIP.')
    parser.add_argument('--csv', help='Export all winning fractions and enclosures.')
    args = parser.parse_args()
    run(optional_table(args.archive) if args.archive else {}, args.csv)
