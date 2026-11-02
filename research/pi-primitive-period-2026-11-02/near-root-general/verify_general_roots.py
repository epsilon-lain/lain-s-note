"""Exact checks for general near-period roots; Python standard library only.

The all-n convergence theorem is proved in the note. Finite rational bounds
here require no floating-point roots or numerical value of pi.
"""
from fractions import Fraction as F
from math import factorial, comb, prod
from pathlib import Path
import json


def det(matrix):
    a = [list(map(F, row)) for row in matrix]
    if not a:
        return F(1)
    previous, sign = F(1), 1
    for k in range(len(a)-1):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return F(0)
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        diagonal = a[k][k]
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                a[i][j] = (a[i][j]*diagonal-a[i][k]*a[k][j])/previous
            a[i][k] = F(0)
        previous = diagonal
    return sign*a[-1][-1]


def trim(p):
    while len(p)>1 and not p[-1]:
        p.pop()
    return p


def evaluate(p, x):
    y = F(0)
    for c in reversed(p):
        y = y*x+c
    return y


def interpolate(values):
    out = [F(0)]*len(values)
    basis = [F(1)]
    differences = list(map(F, values))
    for k in range(len(values)):
        for j, c in enumerate(basis):
            out[j] += c*differences[0]
        differences = [b-a for a,b in zip(differences,differences[1:])]
        nxt = [F(0)]*(len(basis)+1)
        for j,c in enumerate(basis):
            nxt[j] -= c*k/(k+1)
            nxt[j+1] += c/(k+1)
        basis = nxt
    return trim(out)


def degree_bound(n):
    return n*n//2 if n%2==0 else (n*n-3)//2


def hankel_period(n, z):
    if z == 0:
        return F(0)
    z = F(z)
    mu = [F(0)]
    for k in range(1,2*n-1):
        mu.append(1/z-k*mu[-1]/z)
    A = prod(factorial(a) for a in range(n))
    return z**(n*n)*det([[mu[i+j] for j in range(n)] for i in range(n)])/A**2


def raw_period(n, z):
    rows = [(j,s) for j in range(2) for s in range(n)]
    cols = [(h,a) for h in range(2) for a in range(n)]
    # Direct ordinary Taylor coefficients, independent of the Hankel recurrence.
    return det([[sum(F(comb(a,k))*F(j*z)**(a-k)*F(h**(s-k),factorial(s-k))
                     for k in range(min(a,s)+1)) for h,a in cols] for j,s in rows])


def reflection_factor(n, z, sign):
    return det([[sum(F(comb(a,k))*F(z,2)**(a-k)*
                     (F(1,2)**(s-k)+sign*(-1)**a*F(-1,2)**(s-k))/factorial(s-k)
                     for k in range(min(a,s)+1)) for a in range(n)] for s in range(n)])


def polynomial_checks():
    records = []
    direct_checks = reflection_checks = 0
    for n in range(2,11):
        bound = degree_bound(n)
        p = interpolate([hankel_period(n,z) for z in range(bound+1)])
        assert len(p)-1 <= bound
        assert all(p[k]==0 for k in range(n))
        assert p[n] == (-1)**(n+1)*(n-1)
        assert all(c==0 for k,c in enumerate(p) if (k-n)%2)
        for z in (-1,1,2):
            value = evaluate(p,z)
            assert value == raw_period(n,z)
            direct_checks += 1
            assert value == reflection_factor(n,z,1)*reflection_factor(n,z,-1)
            reflection_checks += 1
        z=bound+2
        assert evaluate(p,z)==hankel_period(n,z)
        n0,n1=(n+1)//2,n//2
        assert n*(n-1)//2-comb(n0,2)-comb(n1,2)==n0*n1==n*n//4
        records.append({'n':n,'degree':len(p)-1,'proved_degree_bound':bound,
                        'exact_order_at_zero':n,'first_coefficient':str(p[n]),
                        'normalized_degree_bound':bound-n,
                        'polynomial':{str(k):str(c) for k,c in enumerate(p) if c}})
    return {'cases':records,'raw_determinant_checks':direct_checks,
            'reflection_factor_checks':reflection_checks}


def factorial_hankel_checks():
    count = 0
    for n in range(1,13):
        for alpha in (0,2):
            value=det([[factorial(i+j+alpha) for j in range(n)] for i in range(n)])
            expected=prod(factorial(j)*factorial(j+alpha) for j in range(n))
            assert value == expected
            count += 1
        if n>=2:
            A=prod(factorial(a) for a in range(n))
            modified=[[factorial(i+j) if i+j else 0 for j in range(n)] for i in range(n)]
            assert det(modified)==(1-n)*A*A
            count += 1
    return count


def finite_distance_bounds():
    records=[]
    scale=10**6
    for n in (8,12,20,32,48,64):
        A=prod(factorial(a) for a in range(n))
        B=prod(factorial(a) for a in range(2*n))
        D=degree_bound(n)-n
        assert D>0
        for k in (1,2,3):
            upper_t=F(44*k,7)
            v=upper_t**(n*n-n)*F(A*A,(n-1)*B)
            record={'n':n,'period_index':k,'normalized_degree_bound':D,
                    'rational_v_numerator_bits':v.numerator.bit_length(),
                    'rational_v_denominator_bits':v.denominator.bit_length(),
                    'v_less_than_one':v<1}
            if v<1:
                lo,hi=0,scale
                while hi-lo>1:
                    mid=(lo+hi)//2
                    if F(mid,scale)**D>=v:
                        hi=mid
                    else:
                        lo=mid
                q=F(hi,scale)
                assert q**D>=v and F(lo,scale)**D<v
                if q<1:
                    radius=upper_t*q/(1-q)
                    record.update({'certified_power_root_upper':str(q),
                                   'certified_root_distance_upper':str(radius),
                                   'rational_power_comparison_passed':True})
                else:
                    record['grid_too_coarse_for_radius']=True
            records.append(record)
    return records


def main():
    polynomials=polynomial_checks()
    result={'status':'passed','polynomial_checks':polynomials,
            'factorial_hankel_checks':factorial_hankel_checks(),
            'finite_root_distance_bounds':finite_distance_bounds(),
            'scope':'General convergence follows from the proof, factorial asymptotics and root-product lemma. '
                    'The verifier checks finite exact identities and rational sufficient distance bounds. '
                    'No assertion of imaginary/simple roots, optimal convergence rate, or rational approximation.'}
    Path(__file__).with_name('general_root_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','exact_period_polynomials':len(polynomials['cases']),
                      'raw_determinant_checks':polynomials['raw_determinant_checks'],
                      'reflection_factor_checks':polynomials['reflection_factor_checks'],
                      'factorial_hankel_checks':result['factorial_hankel_checks'],
                      'rational_root_distance_records':len(result['finite_root_distance_bounds']),
                      'certified_distance_bounds':sum('certified_root_distance_upper' in r
                                                       for r in result['finite_root_distance_bounds'])},indent=2))


if __name__=='__main__':
    main()
