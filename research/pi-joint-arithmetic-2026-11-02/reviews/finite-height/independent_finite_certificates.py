"""Independent standard-library audit; never modifies the supplied input.

Reflection fixtures are certified from their defining coefficient determinants
at degree-bound-plus-one integer arguments. Modular irreducibility uses Rabin's
criterion implemented here, not a factorization library.
"""
import json
import math
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
SOURCE = BASE / "imported-next-work/arithmetic_route"
OUT = Path(__file__).resolve().parent


def bareiss(a):
    a = [row[:] for row in a]
    n = len(a)
    sign, prev = 1, 1
    for k in range(n - 1):
        if not a[k][k]:
            row = next((r for r in range(k + 1, n) if a[r][k]), None)
            if row is None:
                return 0
            a[k], a[row] = a[row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                v = a[i][j] * pivot - a[i][k] * a[k][j]
                assert v % prev == 0
                a[i][j] = v // prev
            a[i][k] = 0
        prev = pivot
    return sign * a[-1][-1]


def reflection_scaled(n, sign, x):
    # At Z=2x, row s is multiplied by 2^s s!, making every entry integral.
    return [[sum(
        math.comb(a, k) * x ** (a - k) * 2 ** k
        * math.factorial(s) // math.factorial(s-k)
        * (1 + sign * (-1) ** (a + s - k))
        for k in range(min(a, s) + 1))
        for a in range(n)] for s in range(n)]


def trim(a, p):
    a = [x % p for x in a]
    while a and a[-1] == 0:
        a.pop()
    return a


def sub(a, b, p):
    return trim([(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))], p)


def mul(a, b, p):
    if not a or not b:
        return []
    c = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            c[i+j] = (c[i+j] + u*v) % p
    return trim(c, p)


def rem(a, b, p):
    a, b = trim(a, p), trim(b, p)
    assert b
    inv = pow(b[-1], -1, p)
    while len(a) >= len(b):
        j, c = len(a) - len(b), a[-1] * inv % p
        for i in range(len(b)):
            a[i+j] = (a[i+j] - c*b[i]) % p
        a = trim(a, p)
    return a


def gcd(a, b, p):
    while b:
        a, b = b, rem(a, b, p)
    if not a:
        return []
    inv = pow(a[-1], -1, p)
    return trim([x*inv for x in a], p)


def modpow(a, k, f, p):
    ans = [1]
    a = rem(a, f, p)
    while k:
        if k & 1:
            ans = rem(mul(ans, a, p), f, p)
        a = rem(mul(a, a, p), f, p)
        k >>= 1
    return ans


def prime_divisors(n):
    out, d = [], 2
    while d*d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def rabin_irreducible(f, p):
    f = trim(f, p)
    d = len(f)-1
    assert d >= 1 and f[-1] == 1
    x = rem([0, 1], f, p)
    frobenius = [x]
    for _ in range(d):
        frobenius.append(modpow(frobenius[-1], p, f, p))
    if frobenius[d] != x:
        return False
    return all(gcd(f, sub(frobenius[d//ell], x, p), p) == [1]
               for ell in prime_divisors(d))


def main():
    fixtures = json.loads((SOURCE / "factors.json").read_text())
    certs = json.loads((SOURCE / "modular_irreducibility.json").read_text())
    lookup = {}
    determinant_checks = 0
    records = []
    for row in fixtures["records"]:
        n, r = row["n"], row["n"]//2
        A_n = math.prod(math.factorial(a) for a in range(n))
        A_r = math.prod(math.factorial(a) for a in range(r))
        assert A_n % (A_r*A_r) == 0
        B = A_n // (A_r*A_r)
        assert B == int(row["normalization_B"])
        degree_bound = n*(n-1)//2
        scale = 2**degree_bound * A_n
        for factor in row["factors"]:
            label = factor["sign"]
            eps = 1 if label == "+" else -1
            coeffs = {int(k): int(v) for k,v in factor["coefficients"].items()}
            assert max(coeffs) == r*(r-1) and coeffs[max(coeffs)] == 1
            assert all(k % 2 == 0 for k in coeffs)
            normalization_sign = (-1) ** (r*(r-1)//2 if eps == 1 else r*(r+1)//2)
            for x in range(1, degree_bound+2):
                z = 2*x
                fvalue = sum(v*z**k for k,v in coeffs.items())
                target = scale * normalization_sign * z**r * fvalue
                assert target % B == 0
                assert bareiss(reflection_scaled(n, eps, x)) == target // B
                determinant_checks += 1
            expected_raw_constant = (-1)**(r-1)*(n-1) if eps == 1 else (-1)**r
            assert coeffs[0] == normalization_sign*B*expected_raw_constant
            lookup[n,label] = coeffs
            records.append({"n":n,"sign":label,"degree":max(coeffs),
                            "determinant_evaluations":degree_bound+1,
                            "coefficient_identity_certified":True,
                            "signed_constant_term":str(coeffs[0])})

    factor_checks, prime_checks = 0, 0
    modular_records = []
    for record in certs["records"]:
        n,label,d = record["n"],record["sign"],record["degree"]
        coeffs = lookup[n,label]
        assert max(coeffs) == d
        possible = set(range(1,d))
        primes = []
        for certificate in record["modular_certificates"]:
            p = certificate["prime"]
            assert p > 1 and all(p % j for j in range(2,math.isqrt(p)+1))
            target = trim([coeffs.get(j,0) for j in range(d+1)],p)
            deriv = trim([j*target[j] for j in range(1,len(target))],p)
            assert gcd(target,deriv,p) == [1]
            product, degrees = [1], []
            for desc in certificate["monic_factors_coefficients_descending"]:
                f = trim(list(reversed(desc)),p)
                assert rabin_irreducible(f,p)
                factor_checks += 1
                product = mul(product,f,p)
                degrees.append(len(f)-1)
            assert product == target
            assert degrees == certificate["degrees"]
            sums = {0}
            for degree in degrees:
                sums = sums | {v+degree for v in sums}
            possible &= sums
            assert sorted(possible) == certificate["remaining_proper_degrees"]
            prime_checks += 1
            primes.append(p)
        assert not possible
        modular_records.append({"n":n,"sign":label,"certificate_primes":primes,
                                "Q_irreducible_certified":True})
    result={"status":"passed","determinant_evaluations":determinant_checks,
            "reflection_records":records,"prime_factorization_certificates":prime_checks,
            "modular_factors_Rabin_certified":factor_checks,
            "irreducibility_records":modular_records,
            "arithmetic":"standard-library exact integers and finite-field polynomials",
            "scope":"Each fixture equals its defining reflection determinant by a polynomial degree bound and enough exact evaluations. All ten stored factors are Q-irreducible. No general-n irreducibility assertion."}
    (OUT/"independent_finite_certificates.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ["reflection_records","irreducibility_records"]},indent=2))


if __name__ == "__main__":
    main()
