"""Check the certificate independently with a different pi identity and a Farey walk."""
import json
import sys
from pathlib import Path

SCALE = 10**360


def check(condition, detail):
    if not condition:
        raise ValueError(detail)


def scaled_atan(d):
    # Each signed term is enclosed by two grid integers. Stop before
    # the first term smaller than one grid unit and bound the signed tail.
    lo = hi = k = 0
    power = d
    while True:
        divisor = (2*k+1)*power
        unit, rem = divmod(SCALE, divisor)
        if unit == 0:
            return (lo-1, hi) if k % 2 else (lo, hi+1)
        upper = unit + bool(rem)
        if k % 2:
            lo, hi = lo-upper, hi-unit
        else:
            lo, hi = lo+unit, hi+upper
        power *= d*d
        k += 1


def verify(data):
    a, b = scaled_atan(2)
    c, d = scaled_atan(3)
    lo, hi = 4*(a+c), 4*(b+d)
    grid = 10**data["pi_grid_digits"]
    check(data["pi_lower_integer"]*SCALE <= lo*grid <= hi*grid
          <= data["pi_upper_integer"]*SCALE, "pi enclosure mismatch")
    target = data["records"][-1]["q"]
    best = None
    records = []

    def consider(p, q):
        nonlocal best
        low, high = q*(q*lo-p*SCALE), q*(q*hi-p*SCALE)
        lower, upper = max(low, -high, 0), max(abs(low), abs(high))
        if best is None or upper < best[0]:
            best = lower, upper
            records.append((p, q, lower*10**15//SCALE, -(-upper*10**15//SCALE)))
        else:
            check(lower > best[1], "undecided score")

    lp, lq, rp, rq = 3, 1, 4, 1
    check(lp*SCALE < lo < hi < rp*SCALE, "initial bracket")
    consider(lp, lq)
    consider(rp, rq)
    steps = 0
    while lq+rq <= target:
        p, q = lp+rp, lq+rq
        consider(p, q)
        if p*SCALE < q*lo:
            lp, lq = p, q
        else:
            check(p*SCALE > q*hi, "unresolved mediant")
            rp, rq = p, q
        steps += 1
    expected = [(r["p"], r["q"], r["cost_floor_1e15"], r["cost_ceil_1e15"])
                for r in data["records"]]
    check(records == expected, "record list mismatch")
    # Derive the CF indices from Euclid on each certified fraction.
    # The two finite expansions have lengths differing by one. The
    # sign of x-p/q selects the convergent index parity (even below).
    for r in data["records"]:
        p, q = r["p"], r["q"]
        num, den = p, q
        length = 0
        while den:
            num, den = den, num % den
            length += 1
        n = length-1
        below = p*SCALE < q*lo
        if (n % 2 == 0) != below:
            n += 1
        check(r["n"] == n, "convergent index mismatch")
    check([(r[0], r[1]) for r in records[:3]] == [(3,1), (22,7), (355,113)]
          and len(records) == 4, "unexpected record pattern")
    check(records[-2][2] > 34*10**11 and records[-1][3] < 25*10**11,
          "isolation lemma margins")
    print(f"PASS: {steps} mediants, {len(records)} records, final q = {target}")


if __name__ == "__main__":
    verify(json.loads(Path(sys.argv[1]).read_text()))
