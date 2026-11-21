"""Check moment cancellation, missing zero prefixes and 29-adic samples."""
from fractions import Fraction as F
from math import comb, isqrt

PAIRS = [(21,20,29), (119,120,169), (697,696,985), (4059,4060,5741)]
WEIGHTS = [F(k, 6759169234) for k in
           [212234236827039075037, -252121705507803843142691,
            8564497927715628053203125, -8312588456444637730797003]]


def power_sum(exponent):
    return sum(w * (F(x,d)**exponent + F(y,d)**exponent)
               for w,(x,y,d) in zip(WEIGHTS, PAIRS))


def valuation(x,p):
    count = 0
    while x % p == 0:
        x //= p
        count += 1
    return count


def main():
    assert sum(WEIGHTS) == 2
    assert all(x*x+y*y == d*d for x,y,d in PAIRS)
    assert all(power_sum(j) == 0 for j in (1,3,5))
    radius = max(F(max(x,y),d) for x,y,d in PAIRS)
    assert radius == F(21,29)
    # Each coefficient of arcsin has modulus <=1. Therefore this bounds
    # |pi-S_n| at every n>=100 without approximation to π.
    assert 2 * sum(map(abs,WEIGHTS)) * radius**203 / (1-radius**2) < 1
    total = F(0)
    zeros = []
    for n in range(101):
        total += F(comb(2*n,n), 4**n * (2*n+1)) * power_sum(2*n+1)
        if total == 0:
            zeros.append(n)
        if n==14:
            expected = 2*n+1 + 1 - valuation(abs(WEIGHTS[0].numerator),29)
            assert valuation(total.denominator,29) == expected
    assert zeros == [0,1,2]

    # At p^2 > 2*n+1 compare the exact folded prefix in the finite field.
    checked = 0
    sums = []
    total = F(0)
    for n in range(65):
        total += F(comb(2*n,n), 4**n * (2*n+1)) * power_sum(2*n+1)
        sums.append(total)
        for p in range(3, 2*n+2, 2):
            if any(p % d == 0 for d in range(2, isqrt(p)+1)):
                continue
            if p*p <= 2*n+1 or any(d % p == 0 for _,_,d in PAIRS):
                continue
            if any(w.denominator % p == 0 for w in WEIGHTS):
                continue
            k = ((2*n+1)//p-1)//2
            difference = p*total - (-1)**((p-1)//2)*sums[k]
            assert difference.denominator % p and difference.numerator % p == 0
            checked += 1
    print("Prime-band congruences checked:", checked)
    print("Only zero partial sums: 0, 1, 2. Non-base denominator rate: 2/7.")


if __name__ == "__main__":
    main()
