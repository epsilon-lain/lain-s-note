"""Exact finite checks for the endpoint singular/exterior profile.

Standard library only. Checks use Gaussian rational arithmetic. There are no
floating point singular-value, eigenvalue, pi, or matrix-square-root estimates.
The all-size inequalities are proved in EXTERIOR_SINGULAR_PROFILE.md.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial
from pathlib import Path
import json


@dataclass(frozen=True)
class G:
    re: F = F(0)
    im: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "re", F(self.re))
        object.__setattr__(self, "im", F(self.im))

    @staticmethod
    def lift(x):
        return x if isinstance(x, G) else G(x)

    def __add__(self, other):
        other = G.lift(other)
        return G(self.re+other.re, self.im+other.im)

    __radd__ = __add__

    def __neg__(self):
        return G(-self.re, -self.im)

    def __sub__(self, other):
        return self+-G.lift(other)

    def __rsub__(self, other):
        return G.lift(other)+-self

    def __mul__(self, other):
        other = G.lift(other)
        return G(self.re*other.re-self.im*other.im,
                 self.re*other.im+self.im*other.re)

    __rmul__ = __mul__

    def norm2(self):
        return self.re**2+self.im**2

    def inv(self):
        d = self.norm2()
        assert d
        return G(self.re/d, -self.im/d)

    def __truediv__(self, other):
        return self*G.lift(other).inv()

    def __pow__(self, k):
        if k < 0:
            return self.inv()**(-k)
        result, base = G(1), self
        while k:
            if k % 2:
                result = result*base
            base = base*base
            k //= 2
        return result

    def __bool__(self):
        return bool(self.re or self.im)


def eye(n):
    return [[G(int(i == j)) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(col) for col in zip(*a)]


def matmul(a, b):
    return [[sum((x*y for x, y in zip(row, col)), G())
             for col in zip(*b)] for row in a]


def determinant(a):
    a = [[G.lift(x) for x in row] for row in a]
    n = len(a)
    out = G(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return G()
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            out = -out
        v = a[j][j]
        out *= v
        for i in range(j+1, n):
            c = a[i][j]/v
            for k in range(j+1, n):
                a[i][k] -= c*a[j][k]
            a[i][j] = G()
    return out


def inverse(a):
    n = len(a)
    a = [[G.lift(x) for x in row]+unit for row, unit in zip(a, eye(n))]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        v = a[j][j]
        a[j] = [x/v for x in a[j]]
        for i in range(n):
            if i != j:
                c = a[i][j]
                a[i] = [x-c*y for x, y in zip(a[i], a[j])]
    return [row[n:] for row in a]


def compound(a, k):
    subsets = list(combinations(range(len(a)), k))
    return [[determinant([[a[i][j] for j in cols] for i in rows])
             for cols in subsets] for rows in subsets]


def trace(a):
    return sum((a[i][i] for i in range(len(a))), G())


def factorial_product(n):
    out = 1
    for r in range(n):
        out *= factorial(r)
    return out


def endpoint_matrices(n, z):
    def u(k):
        return (-1)**k*factorial(k)*z**(-k-1)

    def v(k):
        return sum((F((-1)**ell*factorial(k), factorial(k-ell))
                    *z**(-ell-1) for ell in range(k+1)), G())

    return ([[u(r+s) for s in range(n)] for r in range(n)],
            [[v(r+s) for s in range(n)] for r in range(n)])


def main():
    counts = {key: 0 for key in (
        'pascal_inverse', 'factorial_hankel_factorization',
        'pascal_frobenius_bound', 'hilbert_inverse',
        'hilbert_inverse_diagonal_bound', 'hilbert_determinant',
        'endpoint_pascal_congruence', 'endpoint_determinant',
        'phase_factorization', 'gram_determinant_modulus',
        'rational_diagonal_bounds', 'sorted_diagonal_bounds',
        'exterior_multiplicativity', 'determinant_exterior_expansion',
        'relative_determinant_identity')}

    for n in range(2, 13):
        p = [[G(comb(r, s) if s <= r else 0) for s in range(n)]
             for r in range(n)]
        pinv = [[G((-1)**(r-s)*comb(r, s) if s <= r else 0)
                 for s in range(n)] for r in range(n)]
        assert matmul(p, pinv) == eye(n) == matmul(pinv, p)
        counts['pascal_inverse'] += 1
        f = [[G(factorial(r+s)) for s in range(n)] for r in range(n)]
        dfp = [[factorial(r)*p[r][s] for s in range(n)] for r in range(n)]
        assert f == matmul(dfp, transpose(dfp))
        counts['factorial_hankel_factorization'] += 1
        frob = sum(x.norm2() for row in p for x in row)
        assert frob == sum(comb(2*r, r) for r in range(n)) < 4**n
        assert frob == sum(x.norm2() for row in pinv for x in row)
        counts['pascal_frobenius_bound'] += 1

        gram = [[G(F(1, r+s+1)) for s in range(n)] for r in range(n)]
        hinv = inverse(gram)
        assert matmul(gram, hinv) == eye(n)
        counts['hilbert_inverse'] += 1
        for r in range(n):
            diagonal = (2*r+1)*comb(n+r, n-r-1)**2*comb(2*r, r)**2
            assert hinv[r][r] == G(diagonal)
            assert diagonal <= 2*n*2**(8*n)
            counts['hilbert_inverse_diagonal_bound'] += 1
        delta = F(factorial_product(n)**4, factorial_product(2*n))
        assert determinant(gram) == G(delta)
        counts['hilbert_determinant'] += 1

        if n <= 8:
            for z in (G(1, 2), G(-F(3, 2), F(1, 3)),
                      G(0, 6), G(0, 7), G(0, F(44, 7))):
                u, v = endpoint_matrices(n, z)
                assert matmul(matmul(p, u), transpose(p)) == v
                counts['endpoint_pascal_congruence'] += 1
                target = factorial_product(n)**2*z**(-n*n)
                assert determinant(v) == determinant(u) == target
                counts['endpoint_determinant'] += 1
                if not z.re:
                    t = z.im
                    q = [F(factorial(r))/t**r for r in range(n)]
                    dqp = [[q[r]*p[r][s] for s in range(n)] for r in range(n)]
                    k = matmul(dqp, transpose(dqp))
                    expected = [[G(0, 1)**(r+s)*k[r][s]/z
                                 for s in range(n)] for r in range(n)]
                    assert u == expected
                    counts['phase_factorization'] += 1
                    # det(S V S)=det(V)/det(J); no square root is needed.
                    detb2 = (determinant(v)/G(delta)).norm2()
                    expected_modulus = F(factorial_product(2*n),
                                         factorial_product(n)**2)/t**(n*n)
                    assert detb2 == expected_modulus**2
                    counts['gram_determinant_modulus'] += 1

    for n in range(2, 65):
        for t in (F(6), F(7), F(44, 7)):
            q = [F(factorial(r))/t**r for r in range(n)]
            for r, value in enumerate(q):
                assert F(n**r)/(3*t)**n <= value <= n**r
                counts['rational_diagonal_bounds'] += 1
            ordered = sorted(q, reverse=True)
            for j, value in enumerate(ordered, start=1):
                assert F(n**(n-j))/(3*t)**n <= value <= n**(n-j)
                counts['sorted_diagonal_bounds'] += 1

    for n in range(2, 6):
        u, v = endpoint_matrices(n, G(0, F(44, 7)))
        h = [[v[r][s]-u[r][s] for s in range(n)] for r in range(n)]
        m = matmul(inverse(h), v)
        # The rational formal-period specialization is only an algebra check;
        # it is not the analytic actual-period Gram inequality.
        formal_period = G(0, F(44, 7))**(n*n)*determinant(h)/G(factorial_product(n)**2)
        assert determinant(m)*formal_period == G(1)
        counts['relative_determinant_identity'] += 1
        coefficients = []
        for k in range(n+1):
            assert compound(matmul(h, m), k) == matmul(compound(h, k), compound(m, k))
            counts['exterior_multiplicativity'] += 1
            coefficients.append(trace(compound(m, k)))
        for eta in (G(0), G(1, 1), G(-F(2, 3), F(1, 2))):
            pencil = [[h[r][s]+eta*v[r][s] for s in range(n)] for r in range(n)]
            scalar = sum((c*eta**k for k, c in enumerate(coefficients)), G())
            assert determinant(pencil) == determinant(h)*scalar
            counts['determinant_exterior_expansion'] += 1

    result = {'status': 'passed', 'arithmetic': 'exact Gaussian rationals',
              'counts': counts,
              'scope': 'Finite identities checked. Uniform spectral estimates and root counts are proved in the companion note, not inferred numerically.'}
    output = Path(__file__).with_name('exterior_profile_verification.json')
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
