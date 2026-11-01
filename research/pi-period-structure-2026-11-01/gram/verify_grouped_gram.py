"""Exact Gaussian-rational checks of the grouped Hardy Gram certificate.

Only the Python standard library is used. Rational imaginary nodes test the
kernel identities, not the transcendental value of pi or exp(period)=1.
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
    t=F(27,71)
    log_lower=2*(t+t**3/3+t**5/5)
    assert log_lower>F(4,5)
    phi_lower=F(23,20)-F(121,3834)
    assert phi_lower>F(11,10)
    checked = 0
    formula_cases = []
    for k in (1,2,3):
        for n in (1,2,3):
            radius=F(20)
            nodes=[QI(0,3*j) for j in range(k)]
            rows=[(j,s,0) for j in range(k) for s in range(n)]
            g=gram(rows,[0],nodes,[QI()]*k,radius)
            d=determinant(g)
            expected=confluent_formula(nodes,n,radius)
            assert d==expected and d.im==0 and d.re>0
            assert g==adjoint(g)
            for j in range(k):
                for h in range(k):
                    direct=direct_log_kernel(nodes[j],nodes[h].conjugate(),n-1,n-1,radius)
                    for s in range(n):
                        for ell in range(n):
                            assert direct[s,ell]==log_kernel(nodes[j],nodes[h].conjugate(),s,ell,radius)
                            checked += 1
            formula_cases.append({'K':k,'jets_per_center':n,'rows':k*n,'gram_det':str(d.re)})

    # Complete source-format packet: H=4,v0=2,w0=w1=1,theta=1/4.
    # Its rows are s=0,1,beta=0 at two centers. Columns 1,X,X^2,X^3
    # use no exponential Y factor, so no false periodicity at the rational
    # imaginary node 6i is assumed.
    radius=F(20);nodes=[QI(),QI(0,6)];rows=[(j,s,0) for j in range(2) for s in range(2)]
    types=list(range(4));lp=log_powers(1)
    records=[]
    for delta in (F(0),F(1,100),F(1,10),F(1,2)):
        shifts=[QI(),QI(0,delta)]
        g=gram(rows,types,nodes,shifts,radius)
        gd=determinant(g)
        assert gd.im==0 and gd.re>0
        assert all(determinant([row[:n] for row in g[:n]]).re>0 for n in range(1,5))
        a=[]
        for j,s,beta in rows:
            point=nodes[j]+shifts[j]
            a.append([sum(comb(alpha,d)*point**(alpha-d)*lp[d][s]
                          for d in range(s+1) if d<=alpha) for alpha in range(4)])
        ad=determinant(a)
        assert ad==(nodes[1]+shifts[1])**4
        # Coefficient matrix in the direct-sum Hardy basis: (a,d), a+d=alpha.
        labels=[(a,d) for a in types for d in range(4)]
        c=[[QI(comb(alpha,a)*radius**d if a+d==alpha else 0)
            for alpha in range(4)] for a,d in labels]
        cg=times(adjoint(c),c);cd=determinant(cg)
        z=[[comb(a,beta)*shifts[j]**(a-beta)
            *sum(comb(d,e)*nodes[j]**(d-e)*lp[e][s] for e in range(s+1) if e<=d)/radius**d
            for a,d in labels] for j,s,beta in rows]
        assert times(z,c)==a
        assert ad.norm2()<=gd.re*cd.re
        diagonal_product=F(1)
        for i in range(len(g)):
            assert g[i][i].im==0
            diagonal_product*=g[i][i].re
        assert gd.re<=diagonal_product
        if delta==F(1,100):
            assert gd.re*cd.re/ad.norm2()<F(9,4)
            assert diagonal_product/gd.re>10000
        records.append({'imaginary_shift':str(delta),'actual_det':str(ad.re),
                        'actual_det_squared':str(ad.norm2()),'gram_det':str(gd.re),
                        'coefficient_gram_det':str(cd.re),
                        'bound_squared_over_actual_squared':str(gd.re*cd.re/ad.norm2()),
                        'rowwise_bound_squared_over_grouped_bound_squared':str(diagonal_product/gd.re)})
    result={'status':'passed','confluent_gram_cases':formula_cases,
            'independent_kernel_coefficient_checks':checked,
            'nonzero_shift_certificates':records,
            'log_49_over_22_rational_lower':str(log_lower),
            'large_K_phi_14_rational_lower':str(phi_lower),
            'scope':'Exact kernel and grouped Cauchy-Binet checks at rational imaginary nodes; no claim exp(6i)=1, no pi approximation or interpolation-threshold conclusion.'}
    (Path(__file__).parent/'grouped_gram_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':run()
