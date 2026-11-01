"""Exact checks of transverse-displacement Gram growth.

Gaussian-rational arithmetic helpers are reused from the earlier grouped
Gram verifier. This file is standalone and uses only the standard library.
The rational imaginary base node tests a kernel, not periodicity of exp.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from math import comb, factorial
from pathlib import Path
import json


@dataclass(frozen=True)
class QI:
    re: F = F(0)
    im: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, 're', F(self.re))
        object.__setattr__(self, 'im', F(self.im))

    @staticmethod
    def cast(x):
        return x if isinstance(x, QI) else QI(x)

    def __add__(self, other):
        other = self.cast(other)
        return QI(self.re + other.re, self.im + other.im)
    __radd__ = __add__

    def __neg__(self):
        return QI(-self.re, -self.im)

    def __sub__(self, other):
        return self + -self.cast(other)

    def __rsub__(self, other):
        return self.cast(other) - self

    def __mul__(self, other):
        other = self.cast(other)
        return QI(self.re*other.re-self.im*other.im,
                  self.re*other.im+self.im*other.re)
    __rmul__ = __mul__

    def conjugate(self):
        return QI(self.re, -self.im)

    def norm2(self):
        return self.re*self.re+self.im*self.im

    def __truediv__(self, other):
        other = self.cast(other)
        n = other.norm2()
        p = self*other.conjugate()
        return QI(p.re/n, p.im/n)

    def __rtruediv__(self, other):
        return self.cast(other)/self

    def __pow__(self, n):
        if n < 0:
            return (QI(1)/self)**(-n)
        out = QI(1)
        for _ in range(n):
            out = out*self
        return out

    def __bool__(self):
        return bool(self.re or self.im)


def determinant(a):
    a = [[QI.cast(x) for x in row] for row in a]
    out = QI(1)
    for k in range(len(a)):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return QI()
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            out = -out
        p = a[k][k]
        out = out*p
        for i in range(k+1, len(a)):
            q = a[i][k]/p
            for j in range(k+1, len(a)):
                a[i][j] = a[i][j]-q*a[k][j]
    return out


def times(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def adjoint(a):
    return [[QI.cast(a[j][i]).conjugate() for j in range(len(a))]
            for i in range(len(a[0]))]


def log_powers(n):
    out = [[F(int(k == 0)) for k in range(n+1)]]
    log = [F(0)]+[F((-1)**(k+1), k) for k in range(1, n+1)]
    for d in range(1, n+1):
        prev = out[-1]
        out.append([sum(prev[k]*log[s-k] for k in range(s+1))
                    for s in range(n+1)])
    return out


def derivative_kernel(z, wbar, s, ell, radius):
    """Coefficient [t^s u^ell] 1/(1-(z+t)(wbar+u)/R^2)."""
    r2 = F(radius)**2
    base = 1-z*wbar/r2
    out = QI()
    for c in range(min(s, ell)+1):
        degree = s+ell-c
        mult = factorial(degree)//(factorial(s-c)*factorial(ell-c)*factorial(c))
        out += mult*(wbar/r2)**(s-c)*(z/r2)**(ell-c)*(1/r2)**c/base**(degree+1)
    return out


def log_kernel(z, wbar, s, ell, radius):
    lp = log_powers(max(s, ell))
    return sum(lp[d][s]*lp[e][ell]*derivative_kernel(z, wbar, d, e, radius)
               for d in range(s+1) for e in range(ell+1))


def direct_log_kernel(z, wbar, smax, ellmax, radius):
    """Independent bivariate reciprocal recurrence, not derivative conversion."""
    l = [QI(0)]+[QI(F((-1)**(k+1), k)) for k in range(1, max(smax,ellmax)+1)]
    zcoef = [z]+l[1:smax+1]
    wcoef = [wbar]+l[1:ellmax+1]
    den = {(s,e):-zcoef[s]*wcoef[e]/F(radius)**2
           for s in range(smax+1) for e in range(ellmax+1)}
    den[0,0] += 1
    inv = {}
    for s in range(smax+1):
        for e in range(ellmax+1):
            if s==e==0:
                inv[s,e] = 1/den[0,0]
            else:
                inv[s,e] = -sum(den[a,b]*inv[s-a,e-b]
                                for a in range(s+1) for b in range(e+1)
                                if a or b)/den[0,0]
    return inv


def gram(rows, types, nodes, shifts, radius):
    out = []
    for j,s,beta in rows:
        line = []
        for k,ell,gamma in rows:
            transverse = QI()
            for a in types:
                if a >= beta and a >= gamma:
                    transverse += (comb(a,beta)*comb(a,gamma)
                                   *shifts[j]**(a-beta)
                                   *shifts[k].conjugate()**(a-gamma))
            line.append(transverse*log_kernel(nodes[j],nodes[k].conjugate(),s,ell,radius))
        out.append(line)
    return out


def confluent_formula(nodes, n, radius):
    k = len(nodes)
    result = QI(F(radius)**(-k*n*(n-1)))
    scaled = [z/radius for z in nodes]
    for i in range(k):
        for j in range(i+1,k):
            result *= (scaled[j]-scaled[i]).norm2()**(n*n)
    for i in range(k):
        for j in range(k):
            result /= (1-scaled[i]*scaled[j].conjugate())**(n*n)
    return result



def run():
    radius = F(20)
    omega = QI(0,6)
    delta = QI(0,F(1,100))
    records=[]
    for n in range(1,5):
        rows=[(j,s,0) for j in range(2) for s in range(n)]
        nodes=[QI(),omega]
        shifts=[QI(),delta]
        g_log=gram(rows,[0,1],nodes,shifts,radius)
        g_ord=[[sum(shifts[j]**a*shifts[k].conjugate()**a for a in (0,1))
                *derivative_kernel(nodes[j],nodes[k].conjugate(),s,ell,radius)
                for k,ell,_ in rows] for j,s,_ in rows]
        d_log=determinant(g_log)
        d_ord=determinant(g_ord)
        assert d_log == d_ord and d_ord.im == 0 and d_ord.re>0
        zero=determinant(gram(rows,[0,1],nodes,[QI(),QI()],radius))
        expected_zero=QI(radius**(-2*n*(n-1))
                        *(omega.norm2()/(radius**2-omega.norm2()))**(n*n))
        assert zero==expected_zero and zero==confluent_formula(nodes,n,radius)

        labels=[(a,d) for a in (0,1) for d in range(n)]
        z_ord=[[shifts[j]**a*comb(d,s)*nodes[j]**(d-s)/radius**d
                if s<=d else QI() for a,d in labels] for j,s,_ in rows]
        lp=log_powers(n-1)
        z_log=[[shifts[j]**a
                *sum(comb(d,e)*nodes[j]**(d-e)*lp[e][s]
                     for e in range(s+1) if e<=d)/radius**d
                for a,d in labels] for j,s,_ in rows]
        selected=determinant(z_ord)
        assert selected==determinant(z_log)
        assert selected==delta**n*radius**(-n*(n-1))
        selected_sq=selected.norm2()
        assert selected_sq<=d_ord.re
        ratio=d_ord.re/zero.re
        ratio_lb=delta.norm2()**n*((radius**2-omega.norm2())/omega.norm2())**(n*n)
        assert ratio_lb == selected_sq/zero.re and ratio>=ratio_lb
        records.append({
            'jets_per_node':n, 'rows':2*n,
            'gram_determinant':str(d_ord.re),
            'zero_gram_determinant':str(zero.re),
            'selected_minor_modulus_squared':str(selected_sq),
            'gram_ratio_lower_bound':str(ratio_lb),
            'gram_ratio_lower_bound_approx':float(ratio_lb),
            'actual_gram_ratio_approx':float(ratio),
            'ordinary_and_log_jet_checks':True,
        })
    assert records[-1]['gram_ratio_lower_bound_approx']>1
    report={
        'status':'passed',
        'scope':'kernel-only rational imaginary base; no exp(6i)=1 assertion',
        'base_omega':'6i', 'shift_delta':'i/100','radius':str(radius),
        'transverse_types':[0,1], 'cases':records,
    }
    out=Path(__file__).with_name('RESULTS.json')
    out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','cases':len(records),
                      'n4_ratio_lower_bound':records[-1]['gram_ratio_lower_bound_approx'],
                      'output':str(out)}))


if __name__=='__main__':
    run()
