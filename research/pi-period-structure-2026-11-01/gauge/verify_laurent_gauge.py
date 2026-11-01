"""Exact Laurent gauge, arithmetic, and envelope certificates.

Python standard library only. General statements are proved in LAURENT_GAUGE.md.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, lcm
from pathlib import Path
import json


def ceilq(x):
    x = F(x)
    return -((-x.numerator) // x.denominator)


def ibinom(n, k):
    if n >= 0:
        return comb(n, k) if k <= n else 0
    return (-1) ** k * comb(k - n - 1, k)


def weight(a, w):
    return sum(x*y for x, y in zip(a, w))


def packet(H, w, theta, v0, K):
    rows = []
    for j in range(K):
        for beta in product(*(range(ceilq(H*theta/wi)) for wi in w)):
            for s in range(max(0, ceilq((H-weight(beta, w)/theta)/v0))):
                rows.append((j, s, beta))
    return sorted(rows)


def columns(H, w, w0=F(1)):
    return [(h, a) for a in product(*(range(int(H/wi)+1) for wi in w))
            for h in range(max(0, int((H-weight(a, w))/w0)+1))]


def mul(a, b, N):
    c = [F(0) for _ in range(N+1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b[:N+1-i]):
            c[i+j] += x*y
    return c


def matrix(rr, cc, centers, d):
    N = max(s for j, s, beta in rr)
    max_power = max(sum(a) for h, a in cc)
    powers = {}
    for j in set(j for j, s, beta in rr):
        for i, c in enumerate(centers):
            p = [F(j)*c] + [F((-1)**(k+1), k) for k in range(1, N+1)]
            powers[j, i, 0] = [F(1)] + [F(0)]*N
            for n in range(1, max_power+1):
                powers[j, i, n] = mul(powers[j, i, n-1], p, N)
    out = []
    for j, s, beta in rr:
        row = []
        for h, alpha in cc:
            if any(a < b for a, b in zip(alpha, beta)):
                row.append(F(0)); continue
            p = [F(ibinom(h-d, k)) for k in range(N+1)]
            factor = 1
            for i, (a, b) in enumerate(zip(alpha, beta)):
                factor *= comb(a, b)
                p = mul(p, powers[j, i, a-b], N)
            row.append(factor*p[s])
        out.append(row)
    return out


def pivots(A):
    a = [row[:] for row in A]
    r = 0; out = []
    for c in range(len(a[0])):
        k = next((i for i in range(r, len(a)) if a[i][c]), None)
        if k is None: continue
        a[r], a[k] = a[k], a[r]
        v = a[r][c]; a[r] = [x/v for x in a[r]]
        for i in range(r+1, len(a)):
            v = a[i][c]
            if v: a[i] = [x-v*y for x, y in zip(a[i], a[r])]
        out.append(c); r += 1
        if r == len(a): break
    return out


def det(A):
    a = [row[:] for row in A]; out = F(1)
    for c in range(len(a)):
        k = next((i for i in range(c, len(a)) if a[i][c]), None)
        if k is None: return F(0)
        if k != c: a[c], a[k] = a[k], a[c]; out = -out
        v = a[c][c]; out *= v
        for i in range(c+1, len(a)):
            x = a[i][c]/v
            if x:
                for j in range(c+1, len(a)):
                    a[i][j] -= x*a[c][j]
    return out


def gauge_case(H, w, theta, v0, K, centers):
    rr = packet(H, w, theta, v0, K); cc = columns(H, w)
    A = matrix(rr, cc, centers, 0)
    ix = pivots(A); assert len(ix) == len(rr)
    D = det([[r[i] for i in ix] for r in A]); assert D
    row_lookup = {r:i for i, r in enumerate(rr)}
    checked = 0; gauges = [1, 3, H+3]
    for d in gauges:
        B = matrix(rr, cc, centers, d)
        for i, (j, s, beta) in enumerate(rr):
            for c in range(len(cc)):
                value = sum(ibinom(-d, s-ell)*A[row_lookup[j, ell, beta]][c]
                            for ell in range(s+1))
                assert value == B[i][c]; checked += 1
        assert det([[r[i] for i in ix] for r in B]) == D
        # A fixed integral triangular mixture of the selected columns.
        mixed = [[r[ix[k]] + (r[ix[k-1]] if k else 0)
                  for k in range(len(ix))] for r in B]
        assert det(mixed) == D
    # A sparse source-row subpacket need not preserve its determinant.
    assert ibinom(-1, 1) == -1 and ibinom(0, 1) == 0
    return dict(H=H, weights=list(w), rows=len(rr), columns=len(cc),
                gauges=gauges, determinant=str(D), entries_checked=checked,
                same_minor_and_mixed_minor_unchanged=True)


def ga(a, b): return (a[0]+b[0], a[1]+b[1])
def gm(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def gmul(a, b, N):
    c = [(F(0), F(0)) for _ in range(N+1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b[:N+1-i]):
            c[i+j] = ga(c[i+j], gm(x, y))
    return c


def arithmetic_case():
    H=7; w=(2,3); theta=F(1,2); v0=F(1); K=2
    rr=packet(H,w,theta,v0,K); cc=columns(H,w)
    q=(7,13); p=(22,41)
    T=[ceilq(wi/(theta*v0)) for wi in w]
    L=[lcm(*range(1,t)) if t>1 else 1 for t in T]
    N=max(s for j,s,b in rr); checked=nonreal=negative=0
    for d in [1,3,10]:
        for j,s,beta in rr:
            for h,alpha in cc:
                if any(a<b for a,b in zip(alpha,beta)): continue
                poly=[(F(ibinom(h-d,k)),F(0)) for k in range(N+1)]
                factor=1
                for i,(a,b) in enumerate(zip(alpha,beta)):
                    g=[(F(0),F(0)) for _ in range(N+1)]
                    g[0]=(F(0),F(2*j*p[i],q[i]))
                    for k in range(1,min(T[i],N+1)):
                        g[k]=(F((-1)**(k+1),k),F(0))
                    for _ in range(a-b): poly=gmul(poly,g,N)
                    factor *= comb(a,b)*(q[i]*L[i])**(a-b)
                value=tuple(x*factor for x in poly[s])
                assert all(x.denominator==1 for x in value)
                checked += 1; nonreal += bool(value[1]); negative += (h<d)
    return dict(scaled_entries_checked=checked, nonreal_entries=nonreal,
                negative_exponent_entries=negative, q=list(q), L=L)


def envelope_cases():
    out=[]; point_checks=0
    for H,w0,R,b in product([F(3),F(8),F(24)], [F(1),F(2)],
                             [F(14),F(30)], [F(1),F(7),F(40)]):
        a=R/w0; d_star=H*max(a-b,F(0))/(2*R); d=d_star.numerator//d_star.denominator
        G=max(b,(a+b)/2)
        upper=max(R*d+H*b,R*H/w0-R*d)
        assert upper<=H*G+R and 0<=d<=H/w0
        real_upper=max(R*d_star+H*b,R*H/w0-R*d_star)
        assert real_upper==H*G
        # Convexity bound for all original integer h and an additional rational mesh.
        hs=set(F(h) for h in range(int(H/w0)+1))
        hs.update(H*F(k,12)/w0 for k in range(13))
        for h in hs:
            assert R*abs(h-d)+(H-w0*h)*b<=upper; point_checks+=1
        # Check the real minimizer against all integer gauges in the allowed interval.
        for dd in range(int(H/w0)+1):
            assert max(R*dd+H*b,R*H/w0-R*dd)>=H*G
        out.append(dict(H=str(H),w0=str(w0),R=str(R),b=str(b),
                        gauge=d,real_optimum=str(G),finite_exponent=str(upper)))
    return dict(cases=len(out),column_budget_points=point_checks,
                representative=out[-6:])


def radius_checks():
    checks=0
    for A,kappa,v0 in product([F(1,10),F(1,50)], [F(1),F(3,2)],
                              [F(2),F(5),F(10)]):
        d=1/(2*v0); u0=kappa/A
        assert kappa>A and kappa>2*d and u0>1
        # Before the kink derivative is positive; after it derivative is negative.
        for f in [F(1,4),F(1,2),F(3,4)]:
            u=1+(u0-1)*f
            assert 1<u<u0 and d/u>0; checks+=1
        for f in [F(1),F(3,2),F(2),F(10)]:
            u=u0*f
            assert d/u-A/2<0; checks+=1
    return checks


def main():
    recurrence_checks=0
    for n in range(-12,13):
        for k in range(1,16):
            assert k*ibinom(n,k)==(n-k+1)*ibinom(n,k-1)
            recurrence_checks+=1
    cases=[gauge_case(3,(1,),F(1,2),F(1),1,(F(3),)),
           gauge_case(4,(1,),F(1,2),F(1),2,(F(3),)),
           gauge_case(7,(2,3),F(1,2),F(1),2,(F(3),F(5))),
           gauge_case(5,(1,2),F(1,3),F(2),2,(F(3),F(5)))]
    result=dict(status='passed',gauge_cases=cases,
                integer_binomial_recurrences=recurrence_checks,
                arithmetic=arithmetic_case(),envelope=envelope_cases(),
                radius_derivative_checks=radius_checks(),
                sparse_source_counterexample={'ungauged':0,'gauge_d_1':-1},
                scope='Exact algebra/envelope checks; no independent source interpolation or pi bad-approximability claim.')
    (Path(__file__).with_name('laurent_gauge_verification.json')).write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
