"""Exact finite audits for NEAREST_ROOT_SCALE.md.

No floating-point pi, numerical root finder, or external package is used.
The all-n analytic and asymptotic statements are proved in the note.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import permutations
from math import comb, factorial
from pathlib import Path
import json


@dataclass(frozen=True)
class QI:
    real: F = F(0)
    imag: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "real", F(self.real))
        object.__setattr__(self, "imag", F(self.imag))

    @staticmethod
    def cast(x):
        return x if isinstance(x, QI) else QI(x)

    def __add__(self, x):
        x = self.cast(x)
        return QI(self.real + x.real, self.imag + x.imag)

    __radd__ = __add__

    def __neg__(self):
        return QI(-self.real, -self.imag)

    def __sub__(self, x):
        return self + (-self.cast(x))

    def __rsub__(self, x):
        return self.cast(x) - self

    def __mul__(self, x):
        x = self.cast(x)
        return QI(self.real*x.real - self.imag*x.imag,
                  self.real*x.imag + self.imag*x.real)

    __rmul__ = __mul__

    def norm2(self):
        return self.real**2 + self.imag**2

    def conjugate(self):
        return QI(self.real, -self.imag)

    def __truediv__(self, x):
        x = self.cast(x)
        d = x.norm2()
        return self * QI(x.real/d, -x.imag/d)

    def __rtruediv__(self, x):
        return self.cast(x) / self

    def __pow__(self, k):
        if k < 0:
            return (QI(1)/self)**(-k)
        result = QI(1)
        for _ in range(k):
            result *= self
        return result

    def __bool__(self):
        return bool(self.real or self.imag)


def trim(p):
    p = list(p)
    while len(p) > 1 and not p[-1]:
        p.pop()
    return p


def add(p, q):
    out = [F(0)] * max(len(p), len(q))
    for j, x in enumerate(p):
        out[j] += x
    for j, x in enumerate(q):
        out[j] += x
    return trim(out)


def mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            if x and y:
                out[i+j] += x*y
    return trim(out)


def scale(p, a):
    return trim([a*x for x in p])


def derivative(p):
    return [i*p[i] for i in range(1, len(p))] or [F(0)]


def evaluate(p, z):
    out = QI()
    for a in reversed(p):
        out = out*z + a
    return out


def determinant(a):
    a = [[QI.cast(x) for x in row] for row in a]
    out = QI(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return QI()
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            out = -out
        x = a[i][i]
        out *= x
        for j in range(i+1, len(a)):
            ratio = a[j][i]/x
            for k in range(i+1, len(a)):
                a[j][k] -= ratio*a[i][k]
    return out


def polynomial(n):
    # P_n = det(E T - T E); polynomial arithmetic retains all cancellation.
    e = [[F(1, factorial(i-j)) if i >= j else F(0)
          for j in range(n)] for i in range(n)]
    trans = [[[F(0)]*(j-i) + [F(comb(j, i))] if j >= i else [F(0)]
              for j in range(n)] for i in range(n)]
    matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            entry = [F(0)]
            for k in range(n):
                entry = add(entry, scale(trans[k][j], e[i][k]))
                entry = add(entry, scale(trans[i][k], -e[k][j]))
            row.append(entry)
        matrix.append(row)
    result = [F(0)]
    for perm in permutations(range(n)):
        sign = (-1)**sum(perm[i] > perm[j] for i in range(n) for j in range(i+1, n))
        term = [F(sign)]
        for i, j in enumerate(perm):
            term = mul(term, matrix[i][j])
            if not any(term):
                break
        result = add(result, term)
    return result


def raw_value(n, z):
    matrix = []
    for j in (0, 1):
        for s in range(n):
            row = []
            for h in (0, 1):
                for a in range(n):
                    value = QI()
                    for k in range(min(a, s)+1):
                        # [z^s] exp(hz)(jZ+z)^a, with 0^0=1.
                        value += F(comb(a, k)*h**(s-k), factorial(s-k))*(j*z)**(a-k)
                    row.append(value)
            matrix.append(row)
    return determinant(matrix)


def main():
    out = {"arithmetic": "exact rational/Gaussian rational; Python standard library",
           "scope": "Finite identity audits only; asymptotic and disk statements require the note's proofs."}
    data = []
    raw_checks = 0
    for n in range(2, 7):
        p = polynomial(n)
        assert len(p)-1 <= n*(n-1)
        assert not p[0] and any(p)
        assert all(not a or k % 2 == n % 2 for k, a in enumerate(p))
        if n == 2:
            assert p == [F(0), F(0), F(-1)]
        if n == 3:
            assert p == [F(0), F(0), F(0), F(2)]
        for z in (QI(1), QI(F(2, 3), F(3, 2)), QI(0, 6)):
            assert evaluate(p, z) == raw_value(n, z)
            raw_checks += 1
        data.append({"n": n, "degree": len(p)-1,
                     "coefficients_ascending": [str(x) for x in p]})
    out["raw_vs_compressed_checks"] = raw_checks
    out["polynomials"] = data

    # Independent factored examples, including repeated roots.
    derivative_checks = 0
    root_sets = [
        [QI(0), QI(1), QI(-2)],
        [QI(1, 2), QI(1, 2), QI(F(-1, 3), F(5, 2))],
        [QI(j, (-1)**j * F(j+1, 3)) for j in range(-3, 4)],
    ]
    for roots in root_sets:
        p = [F(1)]
        for root in roots:
            p = mul(p, [-root, F(1)])
        for center in (QI(F(2, 7), F(5, 11)), QI(0, 6), QI(F(7, 3), F(-2, 9))):
            assert all(center != root for root in roots)
            ratio = evaluate(derivative(p), center)/evaluate(p, center)
            assert ratio == sum((1/(center-root) for root in roots), QI())
            delta2 = min((center-root).norm2() for root in roots)
            assert delta2*ratio.norm2() <= len(roots)**2
            derivative_checks += 1
    out["factored_log_derivative_and_distance_checks"] = derivative_checks

    endpoint_checks = 0
    for mu in ([QI(1, 2), QI(3, -1)],
               [QI(F(1, 2), 7), QI(F(1, 2), -3), QI(9)],
               [QI(2, 1)]*4):
        assert all(mu)
        roots = [1 - 1/x for x in mu]
        rho2 = max(x.norm2() for x in mu)
        delta2 = min((x-1).norm2() for x in roots)
        trace = sum(mu, QI())
        assert delta2*rho2 == 1
        assert trace.norm2() <= len(mu)**2*rho2
        assert delta2*trace.norm2() <= len(mu)**2
        endpoint_checks += 1
    out["endpoint_trace_root_checks"] = endpoint_checks

    blaschke_checks = 0
    for a in (QI(F(1, 3), F(1, 4)), QI(F(-2, 5), F(3, 7)), QI(0, F(7, 8))):
        assert a.norm2() < 1
        for slope in (F(-3), F(-1, 2), F(0), F(2, 3), F(7)):
            w = QI((1-slope**2)/(1+slope**2), 2*slope/(1+slope**2))
            assert w.norm2() == 1
            assert (1-a.conjugate()*w).norm2() == (w-a).norm2()
            blaschke_checks += 1
    out["blaschke_boundary_identity_checks"] = blaschke_checks
    out["status"] = "PASS"
    target = Path(__file__).with_name("nearest_root_scale_verification.json")
    target.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "polynomials"}, indent=2))


if __name__ == "__main__":
    main()
