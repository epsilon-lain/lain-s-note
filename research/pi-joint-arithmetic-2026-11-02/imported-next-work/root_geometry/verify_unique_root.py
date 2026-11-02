"""Finite audits for UNIQUE_NEAR_PERIOD_ROOT.md; no asymptotic extrapolation.

Exact audits use the Python standard library. Optional mpmath computations
illustrate finite nearest-period roots; they are explicitly not proofs.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import cache
from math import factorial, prod
from pathlib import Path
import json


@cache
def A(n):
    return prod(factorial(j) for j in range(n))


def delta(n):
    return F(A(n)**4, A(2*n))


def K(n, t):
    return F(factorial(2*n-2)**2*factorial(2*n-1),
             factorial(n-1)**4)/t**(2*n-1)


def B(n, k, t):
    m = n-k
    return delta(m)/delta(n)*A(k)*prod(factorial(2*m+j) for j in range(k))/t**(k*(2*n-k))


@dataclass(frozen=True)
class QI:
    re: F = F(0)
    im: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, 're', F(self.re))
        object.__setattr__(self, 'im', F(self.im))

    @staticmethod
    def cast(v):
        return v if isinstance(v, QI) else QI(v)

    def __add__(self, v):
        v = self.cast(v)
        return QI(self.re+v.re, self.im+v.im)

    __radd__ = __add__

    def __neg__(self):
        return QI(-self.re, -self.im)

    def __sub__(self, v):
        return self + (-self.cast(v))

    def __rsub__(self, v):
        return self.cast(v) - self

    def __mul__(self, v):
        v = self.cast(v)
        return QI(self.re*v.re-self.im*v.im, self.re*v.im+self.im*v.re)

    __rmul__ = __mul__

    def norm2(self):
        return self.re*self.re+self.im*self.im

    def __truediv__(self, v):
        v = self.cast(v)
        return self * QI(v.re/v.norm2(), -v.im/v.norm2())

    def __rtruediv__(self, v):
        return self.cast(v)/self

    def __pow__(self, k):
        if k < 0:
            return (1/self)**(-k)
        out = QI(1)
        for _ in range(k):
            out *= self
        return out

    def __bool__(self):
        return bool(self.re or self.im)


def determinant(a):
    a = [[QI.cast(x) for x in row] for row in a]
    out = QI(1)
    for i in range(len(a)):
        j = next((j for j in range(i, len(a)) if a[j][i]), None)
        if j is None:
            return QI()
        if j != i:
            a[i], a[j] = a[j], a[i]
            out = -out
        pivot = a[i][i]
        out *= pivot
        for j in range(i+1, len(a)):
            ratio = a[j][i]/pivot
            for k in range(i+1, len(a)):
                a[j][k] -= ratio*a[i][k]
    return out


def interpolate(values):
    out = [QI()] * len(values)
    basis = [F(1)]
    diff = list(values)
    for k in range(len(values)):
        for j, c in enumerate(basis):
            out[j] += c*diff[0]
        diff = [b-a for a,b in zip(diff,diff[1:])]
        nxt = [F(0)]*(len(basis)+1)
        for j,c in enumerate(basis):
            nxt[j] -= c*F(k,k+1)
            nxt[j+1] += c/F(k+1)
        basis = nxt
    return out


def exact_audits():
    counts = {'B_recurrence':0, 'B_endpoint':0, 'geometric_majorant':0,
              'Laguerre_determinants':0, 'mixed_coefficients':0}
    for t in (F(6), F(44,7), F(7)):
        for n in range(4,31):
            assert B(n,0,t) == 1 and B(n,1,t) == K(n,t)
            assert B(n,n,t) == F(A(n)**2,1)/delta(n)/t**(n*n)
            counts['B_endpoint'] += 1
            for k in range(n):
                m = n-k
                rising = prod(range(2*m,2*m+k))
                assert B(n,k+1,t)/B(n,k,t) == F(factorial(k),rising)*K(m,t)
                counts['B_recurrence'] += 1
            for k in range(1,n+1):
                assert B(n,k,t) <= K(n,t)*K(n-1,t)**(k-1)
                counts['geometric_majorant'] += 1
            assert K(n-1,t)/K(n,t) == F(n-1,8*(2*n-3)**2*(2*n-1))*t*t
    for m in range(7):
        for k in range(1,8):
            value = determinant([[factorial(2*m+i+j) for j in range(k)] for i in range(k)])
            assert value == QI(A(k)*prod(factorial(2*m+j) for j in range(k)))
            counts['Laguerre_determinants'] += 1
    # Synthetic mixed model: physical measure du (unit modulus weight),
    # endpoint parameter Z=6i. It retains the exact mixed-insertion algebra
    # and the nonnegative majorant used in the general proof. It is not the
    # actual e^(2*pi*i*u) period determinant.
    t = F(6)
    z = QI(0,t)
    for n in range(1,7):
        J = [[F(1,i+j+1) for j in range(n)] for i in range(n)]
        V = [[sum((F((-1)**ell*factorial(i+j),factorial(i+j-ell))*z**(-ell-1)
                   for ell in range(i+j+1)), QI()) for j in range(n)] for i in range(n)]
        values = [determinant([[J[i][j]+eta*V[i][j] for j in range(n)]
                              for i in range(n)])/delta(n) for eta in range(n+1)]
        coeffs = interpolate(values)
        assert coeffs[0] == QI(1)
        assert coeffs[n] == determinant(V)/delta(n)
        for k,c in enumerate(coeffs):
            # e^(6k)<3^(6k), an exact rational sufficient majorant.
            majorant = 3**(6*k)*B(n,k,t)
            assert c.norm2() <= majorant**2
            counts['mixed_coefficients'] += 1
    return counts


def numerical_examples():
    import mpmath as mp
    mp.mp.dps = 180
    examples = []
    for n in (4,6,8,10,12):
        E = mp.matrix(n)
        for i in range(n):
            for j in range(n):
                E[i,j] = mp.mpf(1)/factorial(i-j) if i>=j else 0

        def raw(z):
            trans = mp.matrix(n)
            for i in range(n):
                for j in range(n):
                    trans[i,j] = mp.binomial(j,i)*z**(j-i) if j>=i else 0
            return mp.det(E*trans-trans*E)

        t = 2*mp.pi
        omega = 1j*t
        val = raw(omega)
        deriv = mp.diff(raw,omega)
        newton = -val/deriv
        q = lambda x: mp.re(raw(1j*x)/(1j)**n)
        predicted = t+mp.im(newton)
        root = mp.findroot(q, (predicted, predicted+mp.im(newton)/100),
                           solver='secant', tol=mp.mpf('1e-140'), maxsteps=50)
        offset = 1j*(root-t)
        scale = mp.mpf(factorial(2*n-2))**2*factorial(2*n-1)/(mp.mpf(factorial(n-1))**4*t**(2*n-1))
        examples.append({'n':n,'imaginary_root':mp.nstr(root,45),
                         'root_offset':mp.nstr(root-t,16),
                         'offset_over_exact_Newton':mp.nstr(offset/newton,16),
                         'K_times_distance':mp.nstr(scale*abs(offset),16),
                         'absolute_residual':mp.nstr(abs(raw(1j*root)),5),
                         'scope':'high-precision numerical illustration; not a root isolation proof'})
    return examples


if __name__ == '__main__':
    result = {'status':'PASS','exact':exact_audits(),
              'numeric':numerical_examples(),
              'scope':'Finite identity audits and numerical illustrations; all asymptotic root assertions require the proof.'}
    Path(__file__).with_name('unique_root_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
