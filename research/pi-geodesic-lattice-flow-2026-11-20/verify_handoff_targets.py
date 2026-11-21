"""Run from the research directory with Python 3. All assertions use integers or fractions."""
from fractions import Fraction
from math import comb, gcd
from verify_common_hypotenuse import triangle_data


def prime_powers(value):
    prime = 2
    while prime * prime <= value:
        if value % prime == 0:
            power = 1
            while value % prime == 0:
                power *= prime
                value //= prime
            yield prime, power
        prime += 1
    if value > 1:
        yield value, value


def valuation(value, prime):
    assert value
    count = 0
    while value % prime == 0:
        value //= prime
        count += 1
    return count


def rational_valuation(value, prime):
    return valuation(value.numerator, prime) - valuation(value.denominator, prime)


def ceil_log(value, base):
    power, count = 1, 0
    while power < value:
        power *= base
        count += 1
    return count


def first_pair_term(n):
    m = 2*n + 1
    return Fraction(comb(2*n, n) * (20**m + 21**m), 4**n * m * 29**m)


def predicted_difference(h, j):
    return 2*h - 2*j - valuation(comb(2*j, j), 29) - 2*valuation(2*j, 29)


def main():
    _, hypotenuse, points = triangle_data()
    largest = max(y for _, y in points)
    assert len(points) == 221 and largest == 918959251607736
    assert 10*largest < 9*hypotenuse
    # An integer certificate for (120/169)^(19/20) < 21/29.
    assert 120**19 * 29**20 < 169**19 * 21**20
    print("largest endpoint:", Fraction(largest, hypotenuse), "; points:", len(points))

    checks = 0
    for m in range(21, 62, 2):
        for p, power in prime_powers(m):
            assert m//p in range(1, 20, 2)
            for x in range(power):
                assert pow(x, m, power) == pow(x, m//p, power)
                checks += 1
    print("odd-moment congruences through order 30: residue checks", checks)

    for h in (1, 2):
        n = (29**h - 1)//2
        baseline = rational_valuation(first_pair_term(n), 29)
        block_sum = sum((first_pair_term(k) for k in range(n+1)), Fraction())
        assert rational_valuation(block_sum, 29) == baseline
        for j in range(1, h+2):
            term = first_pair_term(n+j)
            assert rational_valuation(term, 29) - baseline == predicted_difference(h, j)
            block_sum += term
        assert predicted_difference(h, h+1) < 0
    for h in range(2, 201):
        for j in range(1, max(0, h-2*ceil_log(2*h, 29))+1):
            assert predicted_difference(h, j) > 0
        assert predicted_difference(h, h+1) < 0
    print("29-adic changes: direct terms at h=1,2; integer inequalities through h=200")

    # Verify the elementary minimum-denominator formula against full residue cycles.
    for delta in range(2, 56):
        for step in range(delta+1):
            period = delta//gcd(delta, step)
            for t in range(delta):
                minimum = 1
                for _, power in prime_powers(delta):
                    local = gcd(power, step)
                    if t % local:
                        minimum *= power//gcd(power, t)
                actual = min(delta//gcd(delta, t+step*z) for z in range(period))
                assert actual == minimum
    print("minimum-denominator formula: all scalar residue cycles at 2 <= Delta <= 55")


if __name__ == "__main__":
    main()
