#!/usr/bin/env python3
"""Exact checks of short-prefix Padé gcd bounds and local maxima."""

import argparse
from math import comb, gcd, prod


def valuation(value, prime):
    assert value
    order = 0
    while value % prime == 0:
        value //= prime
        order += 1
    return order


def coefficients(n):
    odd_factorial = 1
    factors = [1]
    for k in range(1, n+1):
        odd_factorial *= 2*k-1
        factors.append(odd_factorial)
    return [comb(2*n-k, k)*factors[n-k] for k in range(n+1)]


def primes_up_to(limit):
    primes = []
    for v in range(2, limit+1):
        if all(v % p for p in primes if p*p <= v):
            primes.append(v)
    return primes


def eval_pair(coeffs, residue, modulus=0):
    real = imag = 0
    for a in reversed(coeffs):
        real, imag = a-residue*imag, residue*real
        if modulus:
            real %= modulus
            imag %= modulus
    return real, imag


def local_maximum(coeffs, prime):
    # q is a local unit. Write t=p/q=prime*u. Removing the common
    # coefficient content lets the lifting tree begin at modest precision.
    scaled = [a*prime**k for k,a in enumerate(coeffs)]
    shift = min(valuation(a, prime) for a in scaled)
    stripped = [a//prime**shift for a in scaled]
    classes = [0]
    modulus = 1
    depth = 0
    while classes:
        high = modulus*prime
        reduced = [a % high for a in stripped]
        classes = [
            r for base in classes for r in range(base, high, modulus)
            if eval_pair(reduced, r, high) == (0,0)
        ]
        if classes:
            depth += 1
        modulus = high
        if len(classes) > 100000:
            raise RuntimeError("residue budget exhausted")
    return shift + depth


def recurrence(p, q, end):
    data = [(1,0), (q,p)]
    for n in range(2, end+1):
        a,b = data[-1]
        c,d = data[-2]
        data.append(((2*n-1)*q*a-p*p*c, (2*n-1)*q*b-p*p*d))
    return data[:end+1]


def check_pair(p, q):
    end = 30
    data = recurrence(p, q, end)
    for n in range(1, end+1):
        a,b = data[n]
        old_a,old_b = data[n-1]
        assert b*old_a-a*old_b == p**(2*n-1)
        coeffs = coefficients(n)
        assert (a,b) == eval_pair([c*q**(n-k) for k,c in enumerate(coeffs)], p)
        h = 0
        while 3**(h+1) <= 2*n-1:
            h += 1
        length = min(n,2*h)
        for prime in primes_up_to(p):
            if prime == 2 or p % prime:
                continue
            content_order = valuation(gcd(a,b), prime)
            base = valuation(coeffs[0], prime)
            if content_order <= base:
                continue
            prefix = 0
            nu = valuation(p, prime)
            has_witness = False
            for k in range(1, length+1):
                prefix += valuation(2*n-2*k+1, prime)
                if k >= 2 and k*nu <= prefix:
                    has_witness = True
            assert has_witness

        if n >= 3 and 6*q <= 2*p <= 7*q:
            for m in range(n+1,end+1):
                c,d = data[m]
                det = d*a-c*b
                assert det > 0
                assert det % (gcd(a,b)*gcd(c,d)) == 0
                assert det <= p**(2*n+1)*prod((2*j-1)*q for j in range(n+2,m+1))


def polynomial_value(coeffs, x, modulus):
    value = 0
    for coefficient in reversed(coeffs):
        value = (value*x + coefficient) % modulus
    return value


def hensel_example(prime, multiplier):
    n = (multiplier*prime*prime+1)//2
    coeffs = coefficients(n)
    baseline = valuation(coeffs[0],prime)
    target = ((2*n+1)*(prime-2)+prime-2)//(prime-1)
    scaled = [c*prime**k for k,c in enumerate(coeffs)]
    divisor = prime**baseline
    assert all(c % divisor == 0 for c in scaled)
    reduced = [c//divisor for c in scaled]
    real = [(-1)**(k//2)*c if k%2 == 0 else 0 for k,c in enumerate(reduced)]
    derivative = [k*c for k,c in enumerate(real) if k >= 1]
    root = next(r for r in range(1,prime) if (r*r+2*multiplier) % prime == 0)
    assert polynomial_value(real,root,prime) == 0
    assert polynomial_value(derivative,root,prime) != 0

    precision = 1
    while precision < target-baseline:
        precision = min(2*precision,target-baseline)
        modulus = prime**precision
        value = polynomial_value(real,root,modulus)
        slope = polynomial_value(derivative,root,modulus)
        root = (root-value*pow(slope,-1,modulus)) % modulus
        assert polynomial_value(real,root,modulus) == 0

    modulus = prime**(target-baseline)
    assert eval_pair(reduced,root,modulus) == (0,0)
    return n,baseline,target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-order", type=int, default=40)
    args = parser.parse_args()

    p,q,n = 51,16,5
    a = coefficients(n)
    assert (945*q**5-a[2]*p**2*q**3) == -3483648000
    real,imag = recurrence(p,q,n)[n]
    assert (real,imag) == (-1859999760,-62126109)
    assert valuation(-3483648000,3) == 5
    assert valuation(gcd(real,imag),3) == 7
    print("Lowest-layer equality fails at (n,p,q)=(5,51,16): predicted 5, actual 7.")

    for prime, multiplier in ((5,3),(7,5),(11,1),(11,9),(13,11),(17,1),(19,1)):
        n,base,lower = hensel_example(prime,multiplier)
        print(f"Constructed l={prime}, n={n}, v_l(H_n)={base}, v_l(c_n)>={lower}.")

    pairs = [(p,q) for p in range(1,75) for q in range(1,36) if gcd(p,q) == 1]
    for p,q in pairs:
        check_pair(p,q)
    print(f"{30*len(pairs)} exact coefficient, prefix and Casoratian checks.")

    exceptions = []
    prime_list = primes_up_to(2*args.max_order-1)
    for n in range(1,args.max_order+1):
        coeffs = coefficients(n)
        for prime in prime_list:
            if prime > 2*n-1:
                break
            if prime != 2:
                maximum = local_maximum(coeffs,prime)
                excess = maximum - valuation(coeffs[0],prime)
                if excess:
                    exceptions.append((n,prime,excess))
    assert [row for row in exceptions if row[0] <= 9] == [(5,3,5)]
    c = gcd(*recurrence(2079840,925121,5)[5])
    assert c == 3**5*945
    print(f"Verified all odd local maxima through n={args.max_order}.")
    print("Positive excess certificates (n, prime, excess):", exceptions)
    assert all(row[2] > 0 for row in exceptions)
    # For n=61 and prime 11 alone the excess already exceeds n*log(3).
    if args.max_order >= 61:
        assert (61,11,109) in exceptions and 11**109 > 3**61


if __name__ == "__main__":
    main()
