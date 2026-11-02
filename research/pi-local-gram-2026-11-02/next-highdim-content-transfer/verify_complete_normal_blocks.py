"""Exact finite verification of complete normal-block products.

Only the Python standard library is used. All claims for arbitrary parameters
are proved in COMPLETE_NORMAL_BLOCK_PRODUCT.md; these checks audit identities.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, factorial, gcd, lcm, prod
from pathlib import Path
import json


def determinant(matrix):
    a = [list(map(F, row)) for row in matrix]
    if not a:
        return F(1)
    old, sign = F(1), 1
    for k in range(len(a)-1):
        pivot = next((i for i in range(k,len(a)) if a[i][k]), None)
        if pivot is None:
            return F(0)
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        diagonal = a[k][k]
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):
                a[i][j] = (a[i][j]*diagonal-a[i][k]*a[k][j])/old
            a[i][k] = F(0)
        old = diagonal
    return sign*a[-1][-1]


def A(n):
    return prod(factorial(j) for j in range(n))


def period_value(n, z):
    if not z:
        return F(0)
    z = F(z)
    mu = [F(0)]
    for k in range(1,2*n-1):
        mu.append(1/z-k*mu[-1]/z)
    return z**(n*n)*determinant([[mu[i+j] for j in range(n)]
                                for i in range(n)])/A(n)**2


def trim(p):
    while len(p)>1 and not p[-1]:
        p.pop()
    return p


def polynomial(n):
    bound = n*(n-1)
    differences = [period_value(n,z) for z in range(bound+1)]
    out = [F(0)]*(bound+1)
    basis = [F(1)]
    for k in range(bound+1):
        for j,c in enumerate(basis):
            out[j] += c*differences[0]
        differences = [b-a for a,b in zip(differences,differences[1:])]
        new = [F(0)]*(len(basis)+1)
        for j,c in enumerate(basis):
            new[j] -= c*F(k,k+1)
            new[j+1] += c/F(k+1)
        basis = new
    return trim(out)


def poly_product(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y:
                    out[i+j] += x*y
    return trim(out)


def evaluate(p,z):
    out=F(0)
    for c in reversed(p):
        out=out*z+c
    return out


def imaginary_evaluate_even(p,x):
    assert all(not c for k,c in enumerate(p) if k%2)
    return sum(c*(-1)**(k//2)*x**k for k,c in enumerate(p) if c)


def indices(N,d):
    return sorted((v for v in product(range(N),repeat=d) if sum(v)<N),
                  key=lambda v:(sum(v),v))


def packet(N,d):
    normals=indices(N,d)
    rows=[(j,s,v) for v in normals for j in (0,1)
          for s in range(2*(N-sum(v)))]
    columns=[(h,a,v) for v in normals for h in (0,1)
             for a in range(2*(N-sum(v)))]
    return rows,columns


def raw_entry(row,column,Z):
    j,s,nu=row
    h,a,kappa=column
    if any(k<v for k,v in zip(kappa,nu)):
        return F(0)
    factor=prod(comb(k,v) for k,v in zip(kappa,nu))
    degree=a+sum(kappa)-sum(nu)
    return factor*sum(F(comb(degree,k))*F(j*Z)**(degree-k)
                      *F(h**(s-k),factorial(s-k))
                      for k in range(min(s,degree)+1))


def binomial_safe(n,k):
    return comb(n,k) if n>=k>=0 else 0


def main():
    base_polys={2*k:polynomial(2*k) for k in range(1,5)}
    for n,p in base_polys.items():
        dn=A(n)**2//A(n//2)**4
        assert len(p)-1==n*n//2
        assert dn==lcm(*(c.denominator for c in p))
        assert dn*p[-1]==(-1)**(n//2)
        assert all(not c for k,c in enumerate(p) if k%2)
        assert all(evaluate(p,z)==period_value(n,z) for z in (F(-3,2),F(7,3)))
    records=[]
    direct_checks=value_denominator_checks=0
    for N,d in ((2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3),(4,1)):
        rows,columns=packet(N,d)
        H=4*N
        theta=F(1,2*N)
        assert F(2,N)/theta==4 and 2/theta==4*N
        assert 2*theta<1 and theta<1
        complete={(j,s,nu) for nu in product(range(N+1),repeat=d)
                  for last in (0,1) for j in (0,1) for s in range(2*N+1)
                  if 2*s+4*sum(nu)+4*N*last<H}
        assert set(rows)==complete
        assert len(rows)==len(columns)==4*comb(N+d,d+1)
        assert all(h+2*a+F(2,N)*sum(v)<=H for h,a,v in columns)
        # Closure of the monic integer triangular X_i -> U_i+T change.
        colset=set(columns)
        for h,a,kappa in columns:
            for lam in product(*(range(k+1) for k in kappa)):
                assert (h,a+sum(kappa)-sum(lam),lam) in colset
        combined=[F(1)]
        min_d=1
        powers={k:comb(N-k+d-1,d-1) for k in range(1,N+1)}
        for k,multiplicity in powers.items():
            assert multiplicity>0
            for _ in range(multiplicity):
                combined=poly_product(combined,base_polys[2*k])
            min_d *= (A(2*k)**2//A(k)**4)**multiplicity
        degree=4*binomial_safe(N+d,d+2)+2*comb(N+d,d+1)
        assert len(combined)-1==degree
        assert degree==2*sum(k*k*r for k,r in powers.items())
        assert min_d==lcm(*(c.denominator for c in combined))
        primitive=[int(c*min_d) for c in combined]
        assert all(c*min_d==v for c,v in zip(combined,primitive))
        assert abs(primitive[-1])==1 and gcd(*primitive)==1
        for Z in (F(1),F(-3,2),F(7,3)):
            raw=determinant([[raw_entry(row,col,Z) for col in columns] for row in rows])
            expected=prod(period_value(2*k,Z)**r for k,r in powers.items())
            assert raw==expected==evaluate(combined,Z) and raw
            direct_checks += 1
        for b in range(2,12):
            for a in range(6*b,7*b+1):
                if gcd(a,b)!=1:
                    continue
                value=imaginary_evaluate_even(primitive,F(a,b))
                assert value and value.denominator==b**degree
                value_denominator_checks += 1
        records.append({'N':N,'normal_dimension_d':d,'spatial_dimension_m':d+1,
                        'rows_and_columns':len(rows),'raw_block_multiplicities':powers,
                        'degree':degree,'exact_minimal_coefficient_denominator':str(min_d),
                        'primitive_unit_leading_coefficient':True,
                        'raw_budget_complete_source_and_triangular_closure':'passed'})
    count_checks=0
    for N in range(2,25):
        for d in range(1,18):
            powers=[comb(N-k+d-1,d-1) for k in range(1,N+1)]
            assert 4*sum(k*r for k,r in enumerate(powers,1))==4*comb(N+d,d+1)
            assert 2*sum(k*k*r for k,r in enumerate(powers,1))==(
                4*binomial_safe(N+d,d+2)+2*comb(N+d,d+1))
            count_checks += 1
    result={'status':'passed','base_even_polynomials':len(base_polys),
            'complete_multidimensional_packets':len(records),
            'independent_raw_determinant_checks':direct_checks,
            'exact_rational_value_denominator_checks':value_denominator_checks,
            'combinatorial_identity_checks':count_checks,'records':records,
            'scope':'All-parameter factorization and opposing limits are proved in the note. '
                    'Finite exact checks do not establish the asymptotics by extrapolation. '
                    'No critical-density transfer or pi approximation claim.'}
    Path(__file__).with_name('complete_normal_blocks_verification.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))


if __name__=='__main__':
    main()
