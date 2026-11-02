"""Exact checks for weighted center coalescence and aggregate support.

Python standard library only. Finite examples do not prove the asymptotic
theorem or any statement about pi's rational approximation.
"""
from fractions import Fraction as F
from heapq import heappop, heappush
from itertools import product
from math import comb, factorial, prod
from pathlib import Path
import json


def ceilq(x):
    return -(-x.numerator // x.denominator)


def cheapest_prefix(weights, count):
    """C_w(n), 0<=n<=count, by a duplicate-free exact lattice heap."""
    zero = (0,) * len(weights)
    heap, seen, out = [(F(0), zero)], {zero}, [F(0)]
    for _ in range(count):
        cost, point = heappop(heap)
        out.append(out[-1] + cost)
        for i, w in enumerate(weights):
            nxt = point[:i] + (point[i] + 1,) + point[i+1:]
            if nxt not in seen:
                seen.add(nxt)
                heappush(heap, (cost+w, nxt))
    return out


def complete_rows(H, K, v0, weights, theta):
    spatial = [b for b in product(*(range(ceilq(theta*H/w)) for w in weights))
               if sum(w*x for w, x in zip(weights, b)) < theta*H]
    return [(j, s, b) for j in range(K)
            for s in range(ceilq(H/v0)) for b in spatial
            if v0*s + sum(w*x for w, x in zip(weights, b))/theta < H]


def lower_cost(rows, weights, K):
    times = sorted({s for j, s, b in rows})
    blocks = []
    for s in times:
        bs = [b for j, t, b in rows if j == 0 and t == s]
        blocks.append((s, len(bs), sum(sum(w*x for w, x in zip(weights, b)) for b in bs)))
    prefix = cheapest_prefix(weights, max((K*b for _, b, _ in blocks), default=0))
    terms = []
    for s, b, spatial_sum in blocks:
        assert prefix[b] == spatial_sum
        term = prefix[K*b] - K*spatial_sum
        assert term >= 0
        terms.append(term)
    return sum(terms), blocks


def scaled_entry(row, col, centers):
    """Ordinary coefficient entry multiplied by s!, at integer centers."""
    j, s, beta = row
    h, alpha = col
    if any(b > a for a, b in zip(alpha, beta)):
        return 0
    coeff = [1]
    factor = prod(comb(a, b) for a, b in zip(alpha, beta))
    for a, b, r in zip(alpha, beta, centers):
        part = [comb(a-b, k)*(j*r)**(a-b-k) for k in range(a-b+1)]
        nxt = [0]*(len(coeff)+len(part)-1)
        for p, x in enumerate(coeff):
            for q, y in enumerate(part):
                nxt[p+q] += x*y
        coeff = nxt
    return factor*sum(coeff[k]*(factorial(s)//factorial(s-k))*h**(s-k)
                      for k in range(min(s, len(coeff)-1)+1))


def determinant(matrix):
    a = [r[:] for r in matrix]
    previous, sign = 1, 1
    for k in range(len(a)-1):
        p = next((i for i in range(k, len(a)) if a[i][k]), None)
        if p is None:
            return 0
        if p != k:
            a[p], a[k] = a[k], a[p]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                num = a[i][j]*pivot-a[i][k]*a[k][j]
                assert num % previous == 0
                a[i][j] = num//previous
            a[i][k] = 0
        previous = pivot
    return sign*a[-1][-1]


def value(rows, cols, centers):
    return F(determinant([[scaled_entry(r, c, centers) for c in cols] for r in rows]),
             prod(factorial(s) for j, s, b in rows))


def trim(p):
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def interpolate(values):
    basis, out = [F(1)], [F(0)]
    diff = list(map(F, values))
    for k in range(len(values)):
        out += [F(0)]*max(0, len(basis)-len(out))
        for i, x in enumerate(basis):
            out[i] += diff[0]*x
        diff = [b-a for a, b in zip(diff, diff[1:])]
        nxt = [F(0)]*(len(basis)+1)
        for i, x in enumerate(basis):
            nxt[i] -= k*x/(k+1)
            nxt[i+1] += x/(k+1)
        basis = nxt
    return trim(out)


def tensor_polynomial(rows, cols):
    m = len(rows[0][2])
    D = tuple(sum(a[i] for h, a in cols)-sum(b[i] for j, s, b in rows) for i in range(m))
    assert all(d >= 0 for d in D)
    grids = [list(range(d+1)) for d in D]
    coeffs = {point: value(rows, cols, point) for point in product(*grids)}
    for i in range(m):
        other = [j for j in range(m) if j != i]
        nxt = {}
        for fixed in product(*(grids[j] for j in other)):
            vals = []
            for k in grids[i]:
                key = [0]*m
                key[i] = k
                for j, v in zip(other, fixed):
                    key[j] = v
                vals.append(coeffs[tuple(key)])
            p = interpolate(vals)
            p += [F(0)]*(len(grids[i])-len(p))
            for k in grids[i]:
                key = [0]*m
                key[i] = k
                for j, v in zip(other, fixed):
                    key[j] = v
                nxt[tuple(key)] = p[k]
        coeffs = nxt
    return D, {a: c for a, c in coeffs.items() if c}


def aggregate(poly, D):
    """Dictionary (zeta multi-index, omega exponent) -> rational coefficient."""
    result = {}
    for rho, c in poly.items():
        delta = tuple(d-r for d, r in zip(D, rho))
        assert all(x >= 0 for x in delta)
        for gamma in product(*(range(d+1) for d in delta)):
            key = (gamma, sum(rho))
            result[key] = result.get(key, F(0)) + c*(-1)**sum(gamma)*prod(
                comb(d, g) for d, g in zip(delta, gamma))
    return {k: c for k, c in result.items() if c}


def polynomial_checks():
    examples = [
        ('two_variable_saturation', F(3), 2, F(2), (F(1), F(2)), F(1,2),
         [(0,(0,0)), (1,(0,0)), (0,(1,0)), (0,(0,1)), (0,(2,0)), (0,(1,1))],
         {(2,1): F(-2)}),
        ('six_row_complete_cancellation', F(6), 2, F(2), (F(2),), F(1,4),
         [(0,(0,)), (1,(0,)), (2,(0,)), (0,(1,)), (1,(1,)), (2,(2,))],
         {(3,): F(3), (4,): F(1)})]
    records = []
    for name, H, K, v0, weights, theta, cols, expected in examples:
        rows = complete_rows(H, K, v0, weights, theta)
        assert len(rows) == len(cols) and len(set(cols)) == len(cols)
        assert all(h+sum(w*a for w,a in zip(weights, alpha)) <= H for h, alpha in cols)
        L, blocks = lower_cost(rows, weights, K)
        D, poly = tensor_polynomial(rows, cols)
        assert poly == expected
        for centers in ((-1,)*len(weights), tuple(3+2*i for i in range(len(weights)))):
            expected_value = sum(c*prod(r**a for r,a in zip(centers,rho)) for rho,c in poly.items())
            assert value(rows, cols, centers) == expected_value
        E = sum(w*d for w,d in zip(weights,D))
        assert min(sum(w*r for w,r in zip(weights,rho)) for rho in poly) == L
        agg = aggregate(poly,D)
        assert max(sum(w*g for w,g in zip(weights,gamma)) for gamma,power in agg) == E-L
        records.append({'name':name,'source_rows':[[j,s,list(b)] for j,s,b in rows],
                        'columns':[[h,list(a)] for h,a in cols], 'weights':list(map(str,weights)),
                        'D':D, 'L_H':str(L),'E_H':str(E),'aggregate_weight_cap':str(E-L),
                        'period_polynomial':{str(rho):str(c) for rho,c in poly.items()},
                        'aggregate_coefficients':[
                            {'zeta_power':g,'omega_power':p,'coefficient':str(c)} for (g,p),c in agg.items()]})
    return records


def ordering_checks():
    count = 0
    for weights in ((F(1),F(2)),(F(2),F(3)),(F(1),F(3),F(4)),
                    (F(1,2),F(2)),(F(2,3),F(5,2))):
        n = 40
        heap_prefix = cheapest_prefix(weights,n)
        cutoff = min(weights)*(n-1)
        costs = sorted(sum(w*x for w,x in zip(weights,a))
                       for a in product(*(range(int(cutoff/w)+1) for w in weights)))
        brute = [F(0)]
        for c in costs[:n]:
            brute.append(brute[-1]+c)
        assert heap_prefix == brute
        count += n+1
        for b in range(11):
            for K in range(1,5):
                assert heap_prefix[K*b] >= K*heap_prefix[b]
                count += 1
    one_dim = 0
    for w in (F(1),F(3,2),F(7)):
        pref = cheapest_prefix((w,),100)
        for b in range(21):
            for K in range(1,5):
                assert pref[K*b]-K*pref[b] == w*F(K*(K-1),2)*b*b
                one_dim += 1
    return {'exact_ordering_and_nonnegativity_checks':count,'one_dimensional_exact_identities':one_dim}


def root_interval(K,m):
    lo,hi=F(1),F(K)
    for _ in range(64):
        mid=(lo+hi)/2
        if mid**m<=K:lo=mid
        else:hi=mid
    assert lo**m <= K <= hi**m
    return lo,hi


def lattice_records():
    records=[]
    for weights,theta,K,v0,heights in (
            ((F(2),),F(1,2),3,F(2),(24,48,96)),
            ((F(1),F(2)),F(1,2),2,F(2),(24,48,96)),
            ((F(1),F(2),F(3)),F(2,3),3,F(2),(12,24,36))):
        m=len(weights)
        rootlo,roothi=root_interval(K,m)
        factor=F(m)*theta/(m+2)
        for H in heights:
            rows=complete_rows(F(H),K,v0,weights,theta)
            L,blocks=lower_cost(rows,weights,K)
            M=len(rows)
            assert M==K*sum(b for s,b,spatial in blocks)
            records.append({'m':m,'H':H,'K':K,'weights':list(map(str,weights)),
                            'theta':str(theta),'M_H':M,'L_H':str(L),
                            'exact_ratio':str(L/(M*H)),
                            'proved_limit_enclosure':[str((rootlo-1)*factor),str((roothi-1)*factor)]})
    return records


if __name__=='__main__':
    result={'status':'passed','polynomial_certificates':polynomial_checks(),
            'ordering_checks':ordering_checks(),'lattice_count_records':lattice_records(),
            'scope':'Exact finite weighted coalescence and reciprocal-aggregate support checks. '
                    'The infinite-height limit is proved in the accompanying note, not inferred '
                    'from these finite records. No determinant magnitude or pi assertion.'}
    Path(__file__).with_name('weighted_cutoff_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','polynomial_certificates':len(result['polynomial_certificates']),
                      'ordering_checks':result['ordering_checks'],
                      'lattice_count_records':len(result['lattice_count_records'])},indent=2))
