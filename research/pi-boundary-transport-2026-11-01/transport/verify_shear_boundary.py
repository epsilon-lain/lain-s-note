"""Exact checks of the shear identity and the weighted boundary remainder.

Only Python standard-library rational arithmetic is used. Rational surrogate
nodes test polynomial identities, not periodicity or any approximation of pi.
"""
from fractions import Fraction as Q
from itertools import product
from math import comb, ceil, factorial
from pathlib import Path
import json


def multiindices(m, total):
    return [a for a in product(range(total + 1), repeat=m) if sum(a) <= total]


def ambient(m, n):
    return [(s, beta) for s in range(n + 1)
            for beta in multiindices(m, n - s)]


def mcomb(a, b):
    return product_number(comb(x, y) for x, y in zip(a, b))


def product_number(xs):
    ans = Q(1)
    for x in xs:
        ans *= x
    return ans


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def power(x, a):
    return product_number(xx ** aa for xx, aa in zip(x, a))


def tentry(source, target, lam):
    s, beta = source
    ell, alpha = target
    gamma = tuple(a - b for a, b in zip(alpha, beta))
    if min(gamma) < 0 or sum(gamma) != s - ell:
        return Q(0)
    return mcomb(alpha, beta) * power(lam, gamma)


def matrix(rows, cols, lam):
    return [[tentry(r, c, lam) for c in cols] for r in rows]


def matmul(a, b):
    if not b:
        return [[] for _ in a]
    n = len(b[0])
    out = [[Q(0) for _ in range(n)] for _ in a]
    for i, row in enumerate(a):
        for k, x in enumerate(row):
            if x:
                for j, y in enumerate(b[k]):
                    if y:
                        out[i][j] += x * y
    return out


def eye(n):
    return [[Q(i == j) for j in range(n)] for i in range(n)]


def determinant(a):
    a = [row[:] for row in a]
    out = Q(1)
    for k in range(len(a)):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return Q(0)
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            out *= -1
        p = a[k][k]
        out *= p
        for i in range(k + 1, len(a)):
            if not a[i][k]:
                continue
            x = a[i][k] / p
            for j in range(k + 1, len(a)):
                a[i][j] -= x * a[k][j]
    return out


def residual_formula(source, target, lam, v0, b, height):
    s, beta = source
    ell, alpha = target
    gamma = tuple(a - bb for a, bb in zip(alpha, beta))
    if min(gamma) < 0 or sum(gamma) != s - ell:
        return Q(0)
    c = tuple(bb - v0 for bb in b)
    gap = height - v0 * s - dot(b, beta)
    assert dot(c, gamma) >= gap
    alt = sum((-1) ** sum(kappa) * mcomb(gamma, kappa)
              for kappa in product(*(range(g + 1) for g in gamma))
              if dot(c, kappa) < gap)
    return mcomb(alpha, beta) * power(tuple(-x for x in lam), gamma) * alt


def polynomial_actual(row, center, lam, hdegree, alpha):
    """Direct expansion of z^h prod_i(z+u_i)^alpha_i at shifted nodes."""
    s, beta = row
    if any(bb > a for a, bb in zip(alpha, beta)):
        return Q(0)
    poly = [comb(hdegree, k) * center ** (hdegree-k)
            for k in range(hdegree + 1)]
    for a, bb, ll in zip(alpha, beta, lam):
        degree = a - bb
        base = center * (1 + ll)
        factor = [comb(a, bb) * comb(degree, k) * base ** (degree-k)
                  for k in range(degree + 1)]
        out = [Q(0) for _ in range(len(poly) + degree)]
        for i, x in enumerate(poly):
            for j, y in enumerate(factor):
                out[i+j] += x * y
        poly = out
    return poly[s] if s < len(poly) else Q(0)


def polynomial_sheared(row, center, lam, hdegree, alpha):
    s, beta = row
    if any(bb > a for a, bb in zip(alpha, beta)):
        return Q(0)
    degree = hdegree + sum(alpha) - sum(beta)
    if s > degree:
        return Q(0)
    return (mcomb(alpha, beta)
            * product_number((1+ll)**(a-bb) for ll,a,bb in zip(lam,alpha,beta))
            * comb(degree, s) * center**(degree-s))


def exponential_local(row, centers, h, alpha):
    """Direct local coefficient of e^(h z) prod_i(center_i+z+u_i)^alpha_i."""
    s,beta=row
    if any(bb>a for a,bb in zip(alpha,beta)):
        return Q(0)
    poly=[Q(1)]
    for center,a,bb in zip(centers,alpha,beta):
        degree=a-bb
        factor=[comb(a,bb)*comb(degree,k)*center**(degree-k)
                for k in range(degree+1)]
        out=[Q(0) for _ in range(len(poly)+degree)]
        for i,x in enumerate(poly):
            for j,y in enumerate(factor):
                out[i+j]+=x*y
        poly=out
    return sum(poly[k]*Q(h**(s-k),factorial(s-k))
               for k in range(min(s,len(poly)-1)+1))


def run():
    cases = [
        {'m':1, 'v0':1, 'b':(3,), 'height':7, 'lam':(Q(1,7),)},
        {'m':2, 'v0':1, 'b':(3,5), 'height':7, 'lam':(Q(1,7),Q(-1,11))},
        {'m':2, 'v0':1, 'b':(1,3), 'height':5, 'lam':(Q(1,7),Q(-1,11))},
        {'m':1, 'v0':2, 'b':(2,), 'height':9, 'lam':(Q(1,7),)},
    ]
    records = []
    substitution_checks = 0
    zero_cost_checks = 0
    one_dimensional_checks = 0
    boundary_checks = 0
    local_exponential_checks = 0
    for case in cases:
        m,v0,b,height,lam = (case[k] for k in ('m','v0','b','height','lam'))
        # Positive weights imply s+|beta| <= (H-1)/min(v0,b).
        n = (height-1)//min((v0,)+b)
        indices = ambient(m,n)
        inside = [r for r in indices if v0*r[0]+dot(b,r[1])<height]
        outside = [r for r in indices if r not in inside]
        c = tuple(bb-v0 for bb in b)
        neg = tuple(-x for x in lam)
        t = matrix(indices,indices,lam)
        tm = matrix(indices,indices,neg)
        assert matmul(t,tm)==eye(len(indices))
        assert matmul(t,t)==matrix(indices,indices,tuple(2*x for x in lam))
        u = matrix(inside,inside,neg)
        ui = matrix(inside,inside,lam)
        v = matrix(inside,outside,neg)
        assert matmul(ui,u)==eye(len(inside))
        assert determinant(u)==1
        assert all(x==0 for row in matrix(outside,inside,neg) for x in row)
        rem = matmul(ui,v)
        for i,source in enumerate(inside):
            gap = height-v0*source[0]-dot(b,source[1])
            for j,target in enumerate(outside):
                expected = residual_formula(source,target,lam,v0,b,height)
                assert rem[i][j]==expected
                boundary_checks += 1
                gamma = tuple(a-bb for a,bb in zip(target[1],source[1]))
                if min(gamma)>=0 and sum(gamma)==source[0]-target[0]:
                    if any(ci==0 and gi>0 for ci,gi in zip(c,gamma)):
                        assert expected==0
                        zero_cost_checks += 1
                    if m==1 and c[0]>0:
                        r = ceil(Q(gap,c[0]))-1
                        g = gamma[0]
                        assert g>r
                        closed = (comb(source[1][0]+g,source[1][0])
                                  *(-lam[0])**g*(-1)**r*comb(g-1,r))
                        assert closed==expected
                        one_dimensional_checks += 1
        # Direct polynomial substitution has a separate computation path.
        columns = [(h,alpha) for h in range(3) for alpha in multiindices(m,4)]
        for center in (Q(0),Q(3),Q(-2,3)):
            actual = [[polynomial_actual(r,center,lam,h,a) for h,a in columns]
                      for r in inside]
            sheared = [[polynomial_sheared(r,center,lam,h,a) for h,a in columns]
                       for r in indices]
            assert matmul(matrix(inside,indices,neg),sheared)==actual
            bi = [[polynomial_sheared(r,center,lam,h,a) for h,a in columns]
                  for r in inside]
            bo = [[polynomial_sheared(r,center,lam,h,a) for h,a in columns]
                  for r in outside]
            rhs = matmul(rem,bo) if outside else [[Q(0) for _ in columns] for _ in inside]
            rhs = [[x+y for x,y in zip(br,rr)] for br,rr in zip(bi,rhs)]
            assert matmul(ui,actual)==rhs, case
            substitution_checks += len(inside)*len(columns)
        # Formula (7a) directly uses local exponential coefficients. No
        # exponential periodicity is required at these surrogate centers.
        scale=tuple(1+x for x in lam)
        zeta=tuple(x/tt for x,tt in zip(lam,scale))
        for j in (0,1,2):
            omega=Q(3)
            actual=[[exponential_local(r,tuple(j*omega*tt for tt in scale),h,a)
                     *power(scale,r[1])/power(scale,a) for h,a in columns]
                    for r in inside]
            period=[[exponential_local(r,(j*omega,)*m,h,a) for h,a in columns]
                    for r in indices]
            assert actual==matmul(matrix(inside,indices,tuple(-x for x in zeta)),period)
            local_exponential_checks+=len(inside)*len(columns)
        records.append({'m':m,'v0':v0,'b':list(b),'H':height,
                        'ambient_rows':len(indices),'source_rows':len(inside),
                        'outside_rows':len(outside),'det_internal_transport':'1'})

    # Two-center square minors check the exact prefactor in (7b).
    determinant_records=[]
    source=[(0,(0,)),(1,(0,))]
    extra=[(0,(1,))]
    columns=[(0,(0,)),(0,(1,)),(0,(3,)),(1,(2,))]
    for omega in (Q(3),Q(6)):
        for ll in (Q(-1,11),Q(1,7),Q(1,1000)):
            tt=1+ll
            zz=ll/tt
            residual=[[residual_formula(r,t,(zz,),1,(3,),2) for t in extra]
                      for r in source]
            actual=[]
            corrected=[]
            for j in range(2):
                actual.extend([[exponential_local(r,(j*omega*tt,),h,a)
                                for h,a in columns] for r in source])
                inside=[[exponential_local(r,(j*omega,),h,a) for h,a in columns]
                        for r in source]
                outside=[[exponential_local(r,(j*omega,),h,a) for h,a in columns]
                         for r in extra]
                tail=matmul(residual,outside)
                corrected.extend([[x+y for x,y in zip(br,tr)]
                                  for br,tr in zip(inside,tail)])
            da=determinant(actual)
            db=determinant(corrected)
            exponent=sum(a[0] for h,a in columns)-2*sum(beta[0] for s,beta in source)
            assert da and da==tt**exponent*db
            determinant_records.append({'surrogate_omega':str(omega),'lambda':str(ll),
                                        'scaling_exponent':exponent,'det_actual':str(da),
                                        'det_corrected_period':str(db)})

    # Exact rational instance of the Chernoff estimate, e^a = 2.
    # This checks the finite positive series majorant, independently of signs.
    bound_checks=0
    c=(2,4)
    z=(Q(1,250),Q(1,500))
    scaled=tuple(zz*2**cc for zz,cc in zip(z,c))
    assert all(x<1 for x in scaled)
    for beta in multiindices(2,3):
        for gap in range(1,11):
            finite=sum(mcomb(tuple(bb+gg for bb,gg in zip(beta,gamma)),beta)
                       *power(z,gamma)
                       for gamma in multiindices(2,16) if dot(c,gamma)>=gap)
            bound=Q(1,2**gap)*product_number((1-zz)**(-bb-1)
                                            for zz,bb in zip(scaled,beta))
            assert finite<=bound
            bound_checks+=1
    result={'status':'passed','packets':records,
            'boundary_coefficient_checks':boundary_checks,
            'one_dimensional_closed_form_checks':one_dimensional_checks,
            'zero_cost_cancellation_checks':zero_cost_checks,
            'direct_polynomial_substitution_entry_checks':substitution_checks,
            'direct_local_exponential_entry_checks':local_exponential_checks,
            'raw_minor_prefactor_certificates':determinant_records,
            'positive_series_bound_checks':bound_checks,
            'scope':'Exact rational identities and analytic-majorant spot checks; no periodicity assertion at surrogate nodes, no high-dimensional determinant estimate or pi approximation conclusion.'}
    target=Path(__file__).with_name('shear_boundary_verification.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    run()
