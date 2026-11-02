"""Standard-library exact checks supplementing the general valuation proofs."""
from fractions import Fraction as F
from math import comb, factorial, gcd
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def v2(value):
    value = abs(value)
    assert value
    return (value & -value).bit_length() - 1


def dfact_odd(n):
    result = 1
    for k in range(1, n + 1, 2):
        result *= k
    return result


def add(p, q):
    result = [F(0)] * max(len(p), len(q))
    for j, a in enumerate(p): result[j] += a
    for j, a in enumerate(q): result[j] += a
    return result


def mul(p, q):
    result = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): result[i + j] += a * b
    return result


def evalpoly(p, x):
    result = F(0)
    for a in reversed(p): result = result * x + a
    return result


def divide_by_x_plus_one(p):
    result = [F(0)] * (len(p) - 1)
    result[-1] = p[-1]
    for j in range(len(p) - 2, 0, -1): result[j - 1] = p[j] - result[j]
    assert p[0] == result[0]
    return result


def primes_upto(n):
    sieve = bytearray(b'\x01') * (n + 1)
    sieve[:2] = b'\x00\x00'
    for p in range(2, n + 1):
        if sieve[p]:
            yield p
            sieve[p * p:n + 1:p] = b'\x00' * len(range(p * p, n + 1, p))


def mod_rational(x, p):
    assert x.denominator % p
    return x.numerator * pow(x.denominator, -1, p) % p


def main():
    # Verify the full v2(A_n) formula independently on a long exact recurrence.
    Aold = A = 1
    for n in range(1, 10001):
        if n > 1:
            k = n - 1
            Aold, A = A, (2 * k + 1) * A + k * k * Aold
        k, residue = divmod(n, 4)
        expected = (2*k, 2*k, 2*k+2, 2*k+3+v2(k+1))[residue]
        assert v2(A) == expected

    # Verify the finite divided-difference polynomial, with no symbolic package.
    R = [F(0)]
    for ell in range(1, 6):
        falling = [F(1)]
        for j in range(ell): falling = mul(falling, [F(-j), F(1)])
        P = [a / (factorial(ell)*dfact_odd(2*ell+1)) for a in mul(falling, falling)]
        P[0] -= evalpoly(P, F(-1))
        R = add(R, divide_by_x_plus_one(P))
    substituted = [F(0)]
    power = [F(1)]
    for coefficient in R:
        substituted = add(substituted, [coefficient*a for a in power])
        power = mul(power, [F(-1), F(2)])
    substituted[0] -= 2
    numerator = [-133890, 106091, -112150, 95560, -57235, 25853, -8200, 1800, -240, 16]
    assert substituted == [F(4*a, 155925) for a in numerator]
    assert all(a.denominator % 2 and a.numerator % 4 == 0 for a in substituted)

    # Store exact normalized solutions and check the digit congruence.
    count = 600
    As, Bs = [1, 1], [0, 1]
    for n in range(1, count):
        As.append((2*n+1)*As[-1]+n*n*As[-2])
        Bs.append((2*n+1)*Bs[-1]+n*n*Bs[-2])
    a = [F(A, factorial(n)) for n, A in enumerate(As)]
    b = [F(B, factorial(n)) for n, B in enumerate(Bs)]
    digit_checks = nonzero_checks = 0
    for p in primes_upto(count):
        if p == 2: continue
        chi = (-1) ** ((p-1)//2)
        for n in range(p, count+1):
            digits = []
            value = n
            while value:
                digits.append(value % p)
                value //= p
            e = len(digits)-1
            rhs = pow(chi, e, p) * mod_rational(b[digits[-1]], p) % p
            for digit in digits[:-1]: rhs = rhs * mod_rational(a[digit], p) % p
            lhs = mod_rational(p**e*b[n], p)
            assert lhs == rhs
            digit_checks += 1
            if rhs:
                g = gcd(As[n], 4*Bs[n])
                expected = sum(n // p**j for j in range(1,e+1))-e
                actual = 0
                while g % p == 0:
                    actual += 1
                    g //= p
                assert actual == expected
                nonzero_checks += 1
    counterexamples=[]
    for n,p in [(72,17),(107,73),(166,131),(170,41),(176,23)]:
        g = gcd(As[n],4*Bs[n]); fact=factorial(n)
        x,y=g,fact;vg=vf=0
        while x%p==0:vg+=1;x//=p
        while y%p==0:vf+=1;y//=p
        assert vg>vf
        counterexamples.append({'n':n,'p':p,'v_p(g)':vg,'v_p(n!)':vf})
    # Independent all-prime-at-once finite check of the lifting cap.
    Aold = A = 1
    Bold, B = 0, 1
    lcm = 1
    for n in range(2, 1001):
        k = n-1
        Aold, A = A, (2*k+1)*A+k*k*Aold
        Bold, B = B, (2*k+1)*B+k*k*Bold
        g = gcd(A, 4*B)
        extra = g // gcd(g, factorial(n))
        extra //= extra & -extra
        lcm = lcm*(n-1)//gcd(lcm, n-1)
        cap = lcm // gcd(lcm, n)
        cap //= cap & -cap
        assert cap % extra == 0
    result={'status':'passed','two_adic_recurrence_range':[1,10000],
            'finite_divided_difference_polynomial':True,
            'digit_congruence_n_max':count,'digit_checks':digit_checks,
            'maximal_loss_valuation_checks':nonzero_checks,
            'uniform_odd_surplus_lcm_cap_range':[2,1000],
            'odd_surplus_counterexamples':counterexamples,
            'scope':'Finite exact checks supplement the general proofs; no subexponential odd-surplus bound is claimed.'}
    (HERE/'gcd_arithmetic_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__': main()
