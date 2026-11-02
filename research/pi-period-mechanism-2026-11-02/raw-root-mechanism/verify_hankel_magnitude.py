"""Exact algebraic certificates for the two-frequency Hankel magnitude theorem.

Python standard library only. No floating-point root or asymptotic inference.
"""
from fractions import Fraction as F
from itertools import permutations
from math import factorial, comb, prod, lcm, gcd
from pathlib import Path
import json
import hashlib
from functools import lru_cache


def trim(p):
    p = list(map(F, p))
    while len(p) > 1 and not p[-1]:
        p.pop()
    return p


def add(a, b):
    out = [F(0)] * max(len(a), len(b))
    for i, x in enumerate(a): out[i] += x
    for i, x in enumerate(b): out[i] += x
    return trim(out)


def scale(a, c):
    return trim([x * c for x in a])


def mul(a, b, cap=None):
    size = len(a) + len(b) - 1
    if cap is not None: size = min(size, cap + 1)
    out = [F(0)] * size
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i+j < size: out[i+j] += x*y
    return trim(out)


def power(a, k):
    out = [F(1)]
    for _ in range(k): out = mul(out, a)
    return out


def divide(a, b):
    a, b = trim(a), trim(b)
    if len(a) < len(b):
        assert a == [0]
        return [F(0)]
    out = [F(0)] * (len(a)-len(b)+1)
    while a != [0] and len(a) >= len(b):
        shift = len(a)-len(b)
        c = a[-1]/b[-1]
        out[shift] += c
        for i, x in enumerate(b): a[i+shift] -= c*x
        a = trim(a)
    assert a == [0], a
    return trim(out)


def evaluate(p, x):
    value = F(0)
    for c in reversed(p): value = value*x+c
    return value


def shift(p, a):
    out = [F(0)]
    for i, c in enumerate(p): out = add(out, scale(power([a, 1], i), c))
    return out


def determinant(matrix):
    a = [list(map(F, row)) for row in matrix]
    out = F(1)
    for i in range(len(a)):
        j = next((j for j in range(i, len(a)) if a[j][i]), None)
        if j is None: return F(0)
        if j != i:
            a[j], a[i] = a[i], a[j]
            out = -out
        pivot = a[i][i]
        out *= pivot
        for j in range(i+1, len(a)):
            c = a[j][i]/pivot
            for k in range(i+1, len(a)): a[j][k] -= c*a[i][k]
    return out


def partitions(k, cap=None):
    if not k:
        yield ()
        return
    if cap is None: cap = k
    for j in range(min(k, cap), 0, -1):
        for tail in partitions(k-j, j): yield (j,)+tail


def cells(lam):
    for i, length in enumerate(lam):
        for j in range(length):
            hook = length-j+sum(other > j for other in lam[i+1:])
            yield i, j, hook


@lru_cache(None)
def tableau_paths(lam):
    if not lam: return 1
    total=0
    for i,length in enumerate(lam):
        if i+1<len(lam) and lam[i+1]==length: continue
        smaller=list(lam)
        smaller[i]-=1
        if smaller[-1]==0: smaller.pop()
        total+=tableau_paths(tuple(smaller))
    return total


def raw_moment(n, k):
    out = F(0)
    for lam in partitions(k):
        if len(lam) > n: continue
        term = F(factorial(k))
        for i, j, h in cells(lam):
            term *= F((n+j-i)**2, h*h*(2*n+j-i))
        out += term
    return out


def moment_rationals():
    x = [0, 1]
    x2 = power(x, 2)
    q1, q9, q25 = [-1, 0, 4], [-9, 0, 4], [-25, 0, 4]
    return {
        2: (x2, scale(q1, 4)),
        4: (scale(mul(x2, [-2, 0, 1]), 3), scale(mul(q9, q1), 16)),
        6: (scale(mul(x2, [-16, 0, 54, 0, -33, 0, 4]), 15),
            scale(mul(mul(q25, q9), power(q1, 2)), 64)),
    }


def symbolic_moment_checks():
    # This denominator contains every content denominator for partitions of <=6.
    common = [F(1)]
    for d in range(-5, 6): common = mul(common, power([d, 2], 2))
    raw_numerators = []
    path_checks=0
    for k in range(7):
        numerator = [F(0)]
        for lam in partitions(k):
            top, bottom, hook_product = [F(1)], [F(1)], 1
            for i, j, h in cells(lam):
                top = mul(top, power([j-i, 1], 2))
                bottom = mul(bottom, [j-i, 2])
                hook_product *= h
            assert tableau_paths(lam)==factorial(k)//hook_product
            assert factorial(k)%hook_product==0
            path_checks+=1
            contribution = mul(top, divide(common, bottom))
            numerator = add(numerator, scale(contribution, F(factorial(k), hook_product**2)))
        raw_numerators.append(numerator)
    for k, (top, bottom) in moment_rationals().items():
        centered = [F(0)]
        for j in range(k+1):
            centered = add(centered, scale(mul(power([0, -F(1,2)], k-j), raw_numerators[j]), comb(k,j)))
        assert mul(centered, bottom) == mul(common, top)
    # Degenerate short partitions are checked directly, not by canceling 0/0.
    small_checks = 0
    for n in range(1, 6):
        raw = [raw_moment(n, k) for k in range(7)]
        for k, (top, bottom) in moment_rationals().items():
            centered = sum(F(comb(k,j))*(-F(n,2))**(k-j)*raw[j] for j in range(k+1))
            assert centered == evaluate(top, n)/evaluate(bottom, n)
            small_checks += 1
    return {'universal_polynomial_moment_identities': 3, 'small_n_centered_moments': small_checks,
            'independent_cell_addition_path_checks':path_checks}


def independent_hilbert_series_checks():
    # Compute det H_n(t) through degree 6 directly by permutation expansion.
    # This route does not use Schur functions or hook lengths.
    checks = 0
    for n in range(1, 6):
        det_series = [F(0)]
        for perm in permutations(range(n)):
            inversions = sum(perm[i] > perm[j] for i in range(n) for j in range(i+1,n))
            term = [F((-1)**inversions)]
            for i, j in enumerate(perm):
                term = mul(term, [F(1, factorial(k)*(i+j+k+1)) for k in range(7)], cap=6)
            det_series = add(det_series, term)
        for k in range(7):
            value = det_series[k]*factorial(k)/det_series[0]
            assert value == raw_moment(n, k)
            checks += 1
    return checks


def lower_bound_certificate():
    common = mul(mul([-25,0,4], [-9,0,4]), power([-1,0,4],2))
    gap_numerator = scale(common, F(5,6))
    t = F(44,7)
    for k, (top, bottom) in moment_rationals().items():
        contribution = mul(top, divide(common, scale(bottom, F(1, {2:4,4:16,6:64}[k]))))
        coefficient = F((-1)**(k//2))*t**k/F(factorial(k)*{2:4,4:16,6:64}[k])
        gap_numerator = add(gap_numerator, scale(contribution, coefficient))
    numerator = scale(gap_numerator, 235298)
    expected = [44118375,0,-186040404,0,131969540,0,-35378248,0,2841056]
    assert numerator == trim(expected)
    shifted = shift(numerator, 3)
    assert shifted == list(map(F,[1908713103,11262098088,21947378508,21140470608,
                                    11464693580,3658868208,680567864,68185344,2841056]))
    assert all(c>0 for c in shifted)
    gap2 = evaluate(numerator,2)/(235298*evaluate(common,2))
    assert gap2 == F(41809379,1111783050)>0
    return {'shifted_positive_coefficients_ascending':[int(c) for c in shifted],
            'n2_gap_above_one_sixth':str(gap2),
            'universal_range':'All integers n>=2'}


def endpoint_checks():
    checks = 0
    examples = []
    for n in range(1, 9):
        A = prod(factorial(a) for a in range(n))
        B = prod(factorial(k) for k in range(2*n))
        H0 = determinant([[F(1,i+j+1) for j in range(n)] for i in range(n)])
        assert H0 == F(A**4,B)
        for z, E in ((F(1),F(1)),(F(2),F(3)),(F(-3,2),F(2,5))):
            mat = []
            for j in range(2):
                for s in range(n):
                    row=[]
                    for h in range(2):
                        for a in range(n):
                            value=sum(F(comb(a,k))* (j*z)**(a-k)*F(h**(s-k),factorial(s-k))
                                      for k in range(min(s,a)+1))
                            row.append(E**(j*h)*value)
                    mat.append(row)
            direct=determinant(mat)
            mu=[(E-1)/z]
            for k in range(1,2*n-1): mu.append(E/z-k*mu[-1]/z)
            hankel=z**(n*n)*determinant([[mu[i+j] for j in range(n)] for i in range(n)])/A**2
            assert direct==hankel
            checks+=1
            if E==1:
                # Integer formal Z makes the row-factorial-cleared value integral.
                assert (direct*A**2).denominator==1
                examples.append({'n':n,'P_at_one':str(direct),'row_factorial_multiplier':str(A**2)})
        assert sum(2*n-2*a+1 for a in range(n+1))==(n+1)**2
        assert all(h+2*a<=2*n for h in range(2) for a in range(n))
    return {'endpoint_hankel_identities':checks,'hilbert_product_identities':8,'examples':examples}


def atan_interval(inverse,terms):
    total=sum(F((-1)**k,(2*k+1)*inverse**(2*k+1)) for k in range(terms))
    next_term=F((-1)**terms,(2*terms+1)*inverse**(2*terms+1))
    return min(total,total+next_term),max(total,total+next_term)


def pi_check():
    a,b=atan_interval(5,20)
    c,d=atan_interval(239,6)
    lo,hi=16*a-4*d,16*b-4*c
    assert 3<lo<hi<F(22,7)
    return {'lower':str(lo),'upper':str(hi),'pi_between_three_and_twenty_two_sevenths':True}


def interval_mul(a,b):
    values=[x*y for x in a for y in b]
    return min(values),max(values)


def interval_evaluate(poly,x):
    value=(F(0),F(0))
    for c in reversed(poly):
        lo,hi=interval_mul(value,x)
        value=(lo+c,hi+c)
    return value


def least_denominator_diagnostic():
    # The companion's own verifier certifies these polynomial coefficients.
    # This supplement independently re-encloses their values at actual 2*pi*i.
    source=Path(__file__).resolve().parent.parent/'root-separation-next/raw_near_period_verification.json'
    data=json.loads(source.read_text())
    assert data['status']=='passed'
    a,b=atan_interval(5,80)
    c,d=atan_interval(239,30)
    x=(32*a-8*d,32*b-8*c)
    records=[]
    for case in data['cases']:
        n=case['n']
        values={int(k):F(v) for k,v in case['polynomial'].items()}
        denominator=lcm(*(v.denominator for v in values.values()))
        assert gcd(*(int(denominator*v) for v in values.values()))==1
        row_multiplier=prod(factorial(s) for s in range(n))**2
        assert row_multiplier%denominator==0
        q=[F(0)]*(max(values)-n+1)
        for k,v in values.items():
            assert k>=n and (k-n)%2==0
            q[k-n]=v*(-1)**((k-n)//2)
        lo,hi=interval_evaluate(q,x)
        lo,hi=interval_mul((lo,hi),(x[0]**n,x[1]**n))
        assert lo*hi>0
        absolute=(lo,hi) if lo>0 else (-hi,-lo)
        cleared=(denominator*absolute[0],denominator*absolute[1])
        assert cleared[0]>1
        # Integer outward enclosures keep this diagnostic compact.
        lower=cleared[0].numerator//cleared[0].denominator
        upper=-(-cleared[1].numerator//cleared[1].denominator)
        records.append({'n':n,'least_coefficient_denominator':str(denominator),
                        'row_multiplier_over_least_denominator':str(row_multiplier//denominator),
                        'primitive_integer_coefficients':True,
                        'cleared_actual_absolute_integer_enclosure':[str(lower),str(upper)],
                        'cleared_actual_absolute_exceeds_one':True})
    assert [r['n'] for r in records]==list(range(4,11))
    return {'source':'../root-separation-next/raw_near_period_verification.json',
            'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'status':'passed','cases':records,
            'scope':'Finite least-polynomial-denominator calculation for n=4..10 only; '
                    'not an asymptotic least-denominator theorem.'}


def main():
    out={'status':'passed',
         'endpoint_checks':endpoint_checks(),
         'moment_checks':symbolic_moment_checks(),
         'independent_hilbert_series_moments':independent_hilbert_series_checks(),
         'uniform_lower_certificate':lower_bound_certificate(),
         'machin_interval':pi_check(),
         'least_denominator_diagnostic':least_denominator_diagnostic(),
         'scope':'Two-frequency complete one-dimensional raw-column family only. '
                 'The actual-period magnitude has a uniform two-sided estimate. '
                 'The specified row-factorial clearing, not every reduced denominator, is analyzed. '
                 'No root-convergence or rational-approximation theorem is certified.'}
    Path(__file__).with_name('hankel_magnitude_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],
                      'endpoint_hankel_identities':out['endpoint_checks']['endpoint_hankel_identities'],
                      'universal_moment_identities':3,
                      'small_n_centered_moments':15,
                      'independent_hilbert_series_moments':out['independent_hilbert_series_moments'],
                      'least_denominator_cases':len(out['least_denominator_diagnostic']['cases']),
                      'uniform_cosine_lower_bound':'1/6 for every n>=2'},indent=2))


if __name__=='__main__': main()
