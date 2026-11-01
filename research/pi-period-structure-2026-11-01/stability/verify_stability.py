"""Exact certificates for period-matrix perturbation; standard library only."""
from fractions import Fraction as F
from dataclasses import dataclass
from math import comb, factorial
from pathlib import Path
import json


@dataclass(frozen=True)
class QI:
    re: F = F(0)
    im: F = F(0)

    @staticmethod
    def of(x):
        return x if isinstance(x, QI) else QI(F(x))

    def __add__(self, x):
        x = self.of(x)
        return QI(self.re + x.re, self.im + x.im)

    __radd__ = __add__

    def __neg__(self):
        return QI(-self.re, -self.im)

    def __sub__(self, x):
        return self + -self.of(x)

    def __rsub__(self, x):
        return self.of(x) + -self

    def __mul__(self, x):
        x = self.of(x)
        return QI(self.re*x.re-self.im*x.im, self.re*x.im+self.im*x.re)

    __rmul__ = __mul__

    def __truediv__(self, x):
        x = self.of(x)
        d = x.re*x.re+x.im*x.im
        if not d:
            raise ZeroDivisionError
        return self*QI(x.re/d, -x.im/d)

    def __pow__(self, n):
        out = QI(F(1))
        for _ in range(n):
            out = out*self
        return out

    def __bool__(self):
        return bool(self.re or self.im)

    def majorant(self):
        return abs(self.re)+abs(self.im)


def matmul(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), QI())
             for j in range(len(b[0]))] for i in range(len(a))]


def identity(n):
    return [[QI(F(int(i == j))) for j in range(n)] for i in range(n)]


def inverse(a):
    n = len(a)
    b = [list(row)+list(ir) for row, ir in zip(a, identity(n))]
    for j in range(n):
        k = next(i for i in range(j, n) if b[i][j])
        b[j], b[k] = b[k], b[j]
        q = b[j][j]
        b[j] = [x/q for x in b[j]]
        for i in range(n):
            if i != j:
                q = b[i][j]
                b[i] = [x-q*y for x, y in zip(b[i], b[j])]
    return [row[n:] for row in b]


def determinant(a):
    n = len(a)
    b = [list(row) for row in a]
    result = QI(F(1))
    for j in range(n):
        candidates = [i for i in range(j, n) if b[i][j]]
        if not candidates:
            return QI()
        k = candidates[0]
        if k != j:
            b[j], b[k] = b[k], b[j]
            result = -result
        q = b[j][j]
        result = result*q
        for i in range(j+1, n):
            f = b[i][j]/q
            for k in range(j+1, n):
                b[i][k] = b[i][k]-f*b[j][k]
    return result


def norm_majorant(a):
    return max(sum(x.majorant() for x in row) for row in a)


def atan_interval(inv, n):
    s = sum((F((-1)**k, (2*k+1)*inv**(2*k+1)) for k in range(n)), F())
    next_term = F((-1)**n, (2*n+1)*inv**(2*n+1))
    return min(s, s+next_term), max(s, s+next_term)


def pi_interval(n=14):
    # Exact tangent identity and branch control for Machin's formula.
    tan2 = F(2,5)/(1-F(1,25))
    tan4 = 2*tan2/(1-tan2*tan2)
    assert (tan4-F(1,239))/(1+tan4*F(1,239)) == 1
    assert 4*(F(1,5)-F(1,375))-F(1,239)>0
    assert F(4,5)<1  # angle < 4/5 < 1 < pi/2; pi>2 by its arctangent integral
    a, b = atan_interval(5, n)
    c, d = atan_interval(239, n)
    return 16*a-4*d, 16*b-4*c


def exp_interval(x, n=22):
    s = sum((F(x)**k/F(factorial(k)) for k in range(n+1)), F())
    first = F(x)**(n+1)/factorial(n+1)
    return s, s+first/(1-F(x, n+2))


# A(r) for rows (s,j,beta)=(0,0,0),(0,1,0),(1,0,0),(1,1,0)
# and columns 1,X,X^3,Y X^2. Coefficients A_0,...,A_3.
COEFF = []
for k in range(4):
    a = [[QI() for _ in range(4)] for _ in range(4)]
    if k == 0:
        a[0][0] = a[1][0] = a[2][1] = a[3][1] = QI(F(1))
    if k == 1:
        a[1][1] = QI(F(1)); a[3][3] = QI(F(2))
    if k == 2:
        a[1][3] = a[3][3] = QI(F(1)); a[3][2] = QI(F(3))
    if k == 3:
        a[1][2] = QI(F(1))
    COEFF.append(a)


def evaluated(r):
    return [[sum((COEFF[k][i][j]*r**k for k in range(4)), QI())
             for j in range(4)] for i in range(4)]


def translated_coefficients(r):
    return [[[sum((COEFF[d][i][j]*comb(d,k)*r**(d-k)
                   for d in range(k,4)), QI())
              for j in range(4)] for i in range(4)] for k in range(4)]


def certificate(p, q, w):
    lo, hi = pi_interval()
    t = F(p, q)
    sigma = QI(F(0), 2*t)
    radius = 2*max(abs(lo-t), abs(hi-t))
    a = evaluated(sigma)
    b = inverse(a)
    assert matmul(b, a) == identity(4)
    shifted = translated_coefficients(sigma)
    assert shifted[0] == a
    for delta in [QI(F(1,7)), QI(F(0),F(-1,13)), QI(F(2,3),F(1,5))]:
        direct = evaluated(sigma+delta)
        rebuilt = [[sum((shifted[k][i][j]*delta**k for k in range(4)), QI())
                    for j in range(4)] for i in range(4)]
        assert direct == rebuilt
    norms = [norm_majorant(matmul(b, shifted[k])) for k in range(1,4)]
    eta = sum((norms[k-1]*radius**k for k in range(1,4)), F())
    assert 0 < eta < 1
    target = F(1,100) if q == 7 else F(1,100000)
    assert eta < target
    det = determinant(a)
    assert det == sigma**4*(1-sigma)
    lower = max(abs(det.re), abs(det.im))*(1-eta)**4
    assert lower > 1
    # Exact source packet: H=4w, v0=2w, w0=w1=w, theta=1/4.
    H, v0 = 4*w, 2*w
    packet = [(s,j,bb) for s in range(4) for j in range(2) for bb in range(3)
              if v0*s+4*w*bb < H]
    assert packet == [(0,0,0),(0,1,0),(1,0,0),(1,1,0)]
    assert all(w*h+w*alpha <= H for h,alpha in [(0,0),(0,1),(0,3),(1,2)])
    assert exp_interval(w-1)[1] < q < exp_interval(w)[0]
    assert max(abs(lo-t), abs(hi-t)) < F(1,q*q)
    return {
        'p':p, 'q':q, 'weight':w, 'rows':4, 'H':H,
        'eta_upper_simple':str(target), 'eta_exact':str(eta),
        'pi_disk_radius_exact':str(radius),
        'preconditioned_coefficient_norms':[str(x) for x in norms],
        'period_determinant_lower_bound_exceeds_one':True,
        'source_packet_and_q_weight_checked':True,
        'approximation_error_below_q_inverse_square':True,
        'scope':'Finite source-format certificate; no separated-weight threshold asserted.'
    }


def staircase_checks():
    results=[]
    # The exact staircase formula is a previously established determinant identity.
    # We check its counts and derive certified inverse-instability lower bounds.
    for n in [2,3,4,20,200,2000,20000]:
        H=4*n
        bs=[(H-s+3)//4 for s in range(H)]
        M=2*sum(bs); W=sum(b*b for b in bs)
        assert M == 4*n*(n+1)
        assert W == F(2,3)*n*(n+1)*(2*n+1)
        assert F(W,M) == F(2*n+1,6)
        denominator = 10 if n==200 else 100 if n==2000 else 1000 if n==20000 else None
        rec={'n':n,'H':H,'rows':M,'determinant_order':W,
             'order_per_row':str(F(W,M))}
        if denominator:
            rho=F(denominator+1,denominator)
            # Norm E >= rho^(W/M)-1 > 100, tested without logarithms or roots.
            assert rho**(2*n+1) > 101**6
            rec.update(relative_center_change=str(rho-1),
                       relative_matrix_perturbation_lower_bound=100)
        results.append(rec)
    return results


def actual_period_instability():
    # Fixed genuine approximation 355/113. The source staircase identity gives
    # det A(r)=r^W for arbitrary G, with H=20n and w1=5.
    lo, hi = pi_interval()
    t = F(355,113)
    assert lo < hi < t
    x = F(1,12000000)
    assert t/hi-1 > x
    n = 200000000
    k = 2*n+1
    # A finite binomial lower bound avoids forming a power with billions of digits.
    lower = sum((F(comb(k,j))*x**j for j in range(41)), F())
    assert lower > 101**6
    M=20*n*(n+1)
    W=F(10,3)*n*(n+1)*(2*n+1)
    assert W.denominator == 1 and F(W,M) == F(k,6)
    assert exp_interval(4)[1] < 113 < exp_interval(5)[0]
    return {'p':355,'q':113,'w0':'1/4','v0':'1','w1':5,
            'theta':'1/4','K':2,'n':n,'H':20*n,'rows':M,
            'determinant_order':int(W),'relative_center_change_lower':str(x),
            'binomial_terms_used':41,
            'relative_matrix_perturbation_strictly_exceeds':100,
            'method':'Exact staircase determinant formula plus Machin interval and finite positive binomial sum; enormous matrix is not materialized.'}


def collision_identity_counts():
    checks=[]
    # Exact finite counts; asymptotic formula is proved in the note.
    for m,w,H,v0,theta,K in [(1,[1],8,F(1),F(1,4),2),
                            (2,[1,2],12,F(1),F(1,2),2),
                            (3,[1,2,3],16,F(2),F(1,2),3)]:
        from itertools import product
        counts=[]
        for beta in product(*(range(int(H*theta/wi)+1) for wi in w)):
            rem=H-sum(wi*bi for wi,bi in zip(w,beta))/theta
            if rem>0:
                q=rem/v0
                per_center=-((-q.numerator)//q.denominator)
                counts.append(K*per_center)
        M=sum(counts); S=sum(x*x for x in counts)
        assert all(n%K==0 for n in counts)
        checks.append({'m':m,'H':H,'rows':M,'sum_n_squared':S,
                       'normalized':str(F(S,M*H)),
                       'limiting_normalization':str(F(2*K,1)/(v0*(m+2)))})
    return checks


def highdim_margins():
    records=[]
    for m, target in [(32,2),(40,11)]:
        theta=F(1,2); base=F(3,5); C=F(13,10)
        power=C**m; K=power.numerator//power.denominator
        w0=base**(-m); v0=2*K*theta**m*w0
        wmin=10**12
        radius_cost=F(14*K)/w0
        # ln(28K)<28K and every w_i>=wmin make the max norm radial.
        assert F(28*K,wmin)<radius_cost
        saving_lower=F(28,81)*(base/theta)**m/(m+2)
        arithmetic_upper=1-m*theta/(m+2)+F(25,9)/(theta*v0)+theta/wmin
        holomorphic_upper=radius_cost+F(25,36)/v0
        gap=saving_lower-arithmetic_upper-holomorphic_upper
        assert gap>target
        assert K*theta**m<1 and K*w0*theta**m/v0==F(1,2)
        records.append({'m':m,'K':K,'w_min_required':wmin,
                        'saving_rational_lower':str(saving_lower),
                        'arithmetic_rational_upper':str(arithmetic_upper),
                        'holomorphic_rational_upper':str(holomorphic_upper),
                        'gap_rational_lower':str(gap),'gap_strictly_exceeds':target,
                        'unknown_interpolation_thresholds_additionally_required':True})
    return records


def main():
    data={'status':'passed', 'period_certificates':[certificate(22,7,2),certificate(355,113,5)],
          'staircase_instability':staircase_checks(),
          'actual_period_instability':actual_period_instability(),
          'highdim_gap_certificates':highdim_margins(),
          'identity_type_counts':collision_identity_counts(),
          'scope':'Exact finite stability certificates and count/formula checks; no uniform inverse bound or pi bad-approximability result.'}
    (Path(__file__).parent/'stability_verification.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'status':'passed','period_certificates':2,
                      'staircase_counts':len(data['staircase_instability']),
                      'identity_count_cases':3,
                      'highdim_gap_certificates':2,
                      'eta_bounds':[x['eta_upper_simple'] for x in data['period_certificates']]}))


if __name__=='__main__': main()
