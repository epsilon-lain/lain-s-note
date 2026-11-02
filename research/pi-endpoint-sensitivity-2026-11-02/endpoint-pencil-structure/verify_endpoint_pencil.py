"""Exact audits for ENDPOINT_PENCIL_STRUCTURE.md; Python standard library only.

The all-dimension spectral and asymptotic claims are proved in the note.
Finite tests do not establish pi's bad approximability or any Z-root claim.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import cache
from math import comb, factorial, prod
from pathlib import Path
import json


@dataclass(frozen=True)
class QI:
    real: F = F(0)
    imag: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, 'real', F(self.real))
        object.__setattr__(self, 'imag', F(self.imag))

    @staticmethod
    def coerce(value):
        return value if isinstance(value, QI) else QI(value)

    def __add__(self, other):
        other = self.coerce(other)
        return QI(self.real + other.real, self.imag + other.imag)

    __radd__ = __add__

    def __neg__(self):
        return QI(-self.real, -self.imag)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        return QI(self.real * other.real - self.imag * other.imag,
                  self.real * other.imag + self.imag * other.real)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        squared = other.real**2 + other.imag**2
        if not squared:
            raise ZeroDivisionError
        return self * QI(other.real / squared, -other.imag / squared)

    def __rtruediv__(self, other):
        return self.coerce(other) / self

    def __pow__(self, exponent):
        if exponent < 0:
            return (QI(1) / self) ** (-exponent)
        value, base = QI(1), self
        while exponent:
            if exponent % 2:
                value *= base
            base *= base
            exponent //= 2
        return value

    def __bool__(self):
        return bool(self.real or self.imag)

    def conjugate(self):
        return QI(self.real, -self.imag)


ZERO, ONE = QI(), QI(1)


def identity(n):
    return [[QI(i == j) for j in range(n)] for i in range(n)]


def matmul(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), ZERO)
             for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a):
    return list(map(list, zip(*a)))


def subtract(a, b):
    return [[x-y for x,y in zip(ar,br)] for ar,br in zip(a,b)]


def pencil(a, b, y):
    return [[y*x-z for x,z in zip(ar,br)] for ar,br in zip(a,b)]


def determinant(matrix):
    a = [[QI.coerce(value) for value in row] for row in matrix]
    result = ONE
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return ZERO
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            result = -result
        diagonal = a[i][i]
        result *= diagonal
        for j in range(i+1, len(a)):
            ratio = a[j][i]/diagonal
            for k in range(i+1, len(a)):
                a[j][k] -= ratio*a[i][k]
    return result


def inverse(matrix):
    n = len(matrix)
    right = identity(n)
    a = [[QI.coerce(value) for value in row] + right[i]
         for i,row in enumerate(matrix)]
    for i in range(n):
        pivot = next((j for j in range(i,n) if a[j][i]), None)
        if pivot is None:
            raise ZeroDivisionError('singular matrix')
        a[i],a[pivot] = a[pivot],a[i]
        diagonal = a[i][i]
        a[i] = [entry/diagonal for entry in a[i]]
        for j in range(n):
            if j != i:
                ratio = a[j][i]
                a[j] = [x-ratio*y for x,y in zip(a[j],a[i])]
    return [row[n:] for row in a]


@cache
def A(n):
    return prod(factorial(j) for j in range(n))


def E(n):
    return [[QI(F(1,factorial(s-a))) if s >= a else ZERO
             for a in range(n)] for s in range(n)]


def T(n, z):
    return [[comb(a,s)*z**(a-s) if a >= s else ZERO
             for a in range(n)] for s in range(n)]


def moments(n, z):
    u = [(-1)**k*factorial(k)/z**(k+1) for k in range(2*n-1)]
    v = [sum((F((-1)**ell*factorial(k),factorial(k-ell))/z**(ell+1)
              for ell in range(k+1)),ZERO) for k in range(2*n-1)]
    return ([[u[i+j] for j in range(n)] for i in range(n)],
            [[v[i+j] for j in range(n)] for i in range(n)])


def raw_g(n, z, y):
    e, t = E(n), T(n,z)
    return determinant(pencil(matmul(e,t),matmul(t,e),y))


def exact_matrix_checks():
    counts = dict(moment_translation=0, commutator_inverse_similarity=0,
                  raw_moment_determinants=0, reciprocal_identity=0,
                  self_inversive_identity=0, relative_pencil_identity=0,
                  relative_determinants_and_trace=0)
    z_values = (QI(F(1,2)), QI(-2), QI(0,1), QI(0,5), QI(0,7))
    y_values = (QI(1), QI(-1), QI(F(3,2)), QI(F(1,2),F(2,3)))
    for n in range(2,8):
        one = identity(n)
        b = [[QI(comb(i,k)) if i >= k else ZERO
              for k in range(n)] for i in range(n)]
        e, ei = E(n), inverse(E(n))
        for z in z_values:
            t, ti = T(n,z), T(n,-z)
            assert matmul(t,ti) == one
            c = matmul(matmul(matmul(ti,ei),t),e)
            cm = matmul(matmul(matmul(t,ei),ti),e)
            assert cm == matmul(matmul(t,inverse(c)),ti)
            assert determinant(c) == ONE
            counts['commutator_inverse_similarity'] += 1
            u,v = moments(n,z)
            assert v == matmul(matmul(b,u),transpose(b))
            assert determinant(u) == determinant(v) == A(n)**2/z**(n*n)
            counts['moment_translation'] += 1
            for y in y_values:
                g = raw_g(n,z,y)
                assert g == determinant(pencil(one,c,y))
                assert g == z**(n*n)*determinant(pencil(v,u,y))/A(n)**2
                counts['raw_moment_determinants'] += 1
                assert raw_g(n,-z,y) == (-1)**n*y**n*raw_g(n,z,1/y)
                counts['reciprocal_identity'] += 1
                if not z.real:
                    assert g == ((-1)**n*y**n*
                                 raw_g(n,z,1/y.conjugate()).conjugate())
                    counts['self_inversive_identity'] += 1
            if not z.real:
                h = subtract(v,u)
                p = raw_g(n,z,ONE)
                assert p
                m = matmul(inverse(h),v)
                assert determinant(m) == determinant(subtract(m,one)) == 1/p
                trace = sum((m[i][i] for i in range(n)),ZERO)
                assert trace.real == F(n,2)
                counts['relative_determinants_and_trace'] += 1
                for eta in y_values:
                    assert raw_g(n,z,1+eta)/p == determinant(pencil(m,
                                      [[-x for x in row] for row in one],eta))
                    counts['relative_pencil_identity'] += 1
        for z in z_values:
            assert raw_g(3,z,ONE) == 2*z**3
    return counts


def L(n):
    return 2*3**8*n**3*2**(8*n)*factorial(2*n-2)


def finite_count_certificates():
    records = []
    for n in (32,64,96,128):
        # alpha=1/2; n is even so the threshold is an exact integer.
        threshold = n**(n//2)
        upper = F(44,7)**(n*n)*F(A(n)**2,A(2*n))
        bound = 6*L(n)
        assert bound > threshold
        # If at most k-1 eigenvalues exceed threshold, the determinant
        # would be at most bound**(k-1)*threshold**(n-k+1).
        # det M has modulus > 1/upper, so ratio < 1 excludes that case.
        ratio = upper*threshold**n
        certified_count = 0
        last_ratio = None
        for k in range(1,n+1):
            if ratio >= 1:
                break
            certified_count = k
            last_ratio = ratio
            ratio *= F(bound,threshold)
        records.append({
            'n':n,
            'alpha':'1/2',
            'threshold_formula':'n^(n/2)',
            'endpoint_radius_formula':'n^(-n/2)',
            'certified_root_count_with_algebraic_multiplicity':certified_count,
            'last_comparison_numerator_bits':last_ratio.numerator.bit_length()
                if last_ratio else None,
            'last_comparison_denominator_bits':last_ratio.denominator.bit_length()
                if last_ratio else None,
            'comparison_is_exact_rational':True,
        })
    return records


def main():
    result = {
        'status':'passed',
        'n_matrix_range':[2,7],
        **exact_matrix_checks(),
        'finite_endpoint_root_count_certificates':finite_count_certificates(),
        'scope':('Exact finite algebraic checks and sufficient endpoint-root '
                 'counts. The general identities and asymptotic count are '
                 'proved in the note. Roots are in independent Y at fixed '
                 'Z=2*pi*i, not zeros in Z and not rational approximants.'),
    }
    target = Path(__file__).with_name('endpoint_pencil_verification.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
