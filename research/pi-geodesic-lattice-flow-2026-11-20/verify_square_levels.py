"""Exact local sieve for the two square shapes in CRT_GRID_TWO_AUDITS.zh.md.

Run without arguments to verify levels 2 through 2048. Only the standard
library is used. The reduced fraction is V_n * z_(n,0) = a/P.
"""
from fractions import Fraction
from math import comb, gcd, isqrt, lcm


def primes_through(limit):
    sieve = bytearray([1]) * (limit + 1)
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p::p] = bytes((limit-p*p)//p + 1)
    return [p for p in range(2, limit + 1) if sieve[p]]


def local_test(n, residue, modulus, exponent, primes):
    s_mod_3 = (pow(2, exponent, 3) * modulus) % 3
    allow_a = s_mod_3 == 0
    allow_b = s_mod_3 == 1
    for prime in primes:
        if not (allow_a or allow_b) or prime > 4*n+1:
            break
        if modulus % prime:
            continue
        power = (prime-1)//2
        if allow_a:
            allow_a = pow(residue, power, prime) == 1
        if allow_b:
            allow_b = prime != 3 and pow(3*residue, power, prime) == 1
    return allow_a, allow_b


def verify():
    limit = 2048
    primes = primes_through(4*limit+1)
    denominator = lcm(*range(1, 4*limit+2, 2))
    harmonic = [0]
    for m in range(1, 2*limit+1):
        harmonic.append(harmonic[-1] + (-1)**(m-1)*(denominator//(2*m-1)))

    survivors = []
    for n in range(2, limit+1):
        v = ((4*n+1)*comb(4*n, 2*n)) // (2**n.bit_count())
        coeff, total = 1, 0
        for i in range(n+1):
            total += coeff*harmonic[n+i]
            coeff = coeff*(n-i)//(i+1)
        a, modulus = Fraction(4*v*total, denominator << n).as_integer_ratio()
        exponent = 2+n-n.bit_count()
        residue = -((-1)**n)*a % modulus
        assert modulus % 2 == 1 and gcd(residue, modulus) == 1
        allow_a, allow_b = local_test(n, residue, modulus, exponent, primes)
        if allow_a or allow_b:
            survivors.append((n, modulus, exponent, residue, allow_a, allow_b))

    assert survivors == [(5, 3, 5, 1, True, False), (6, 99, 6, 64, True, False)]
    for n, modulus, exponent, residue, _, _ in survivors:
        S = modulus << exponent
        squares = [x*x for x in range(1, isqrt(S-1)+1, 2)
                   if (x*x-residue) % modulus == 0]
        assert all(3*t < 2*S for t in squares)
        print(f'n={n}, surviving t before positivity check: {squares}')
    print(f'Both positive square shapes excluded at every level 2..{limit}.')


if __name__ == '__main__':
    verify()
