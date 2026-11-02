"""Independent modular audit of the complete raw normal-block product."""
from fractions import Fraction
from itertools import product
from math import comb, factorial, prod
from pathlib import Path
import json


def indices(d, bound):
    return sorted((x for x in product(range(bound+1), repeat=d) if sum(x)<=bound),
                  key=lambda x:(sum(x),x))


def determinant_mod(matrix, p):
    a=[list(row) for row in matrix]
    out=1
    for k in range(len(a)):
        pivot=next((j for j in range(k,len(a)) if a[j][k]%p), None)
        if pivot is None:
            return 0
        if pivot!=k:
            a[k],a[pivot]=a[pivot],a[k]
            out=-out
        diagonal=a[k][k]%p
        out=out*diagonal%p
        inv=pow(diagonal,-1,p)
        for j in range(k+1,len(a)):
            ratio=a[j][k]*inv%p
            for l in range(k+1,len(a)):
                a[j][l]=(a[j][l]-ratio*a[k][l])%p
    return out%p


def hankel_period(n,z,p):
    inv=pow(z,-1,p)
    mu=[0]
    for k in range(1,2*n-1):
        mu.append((1-k*mu[-1])*inv%p)
    A=prod(factorial(j) for j in range(n))%p
    H=[[mu[i+j] for j in range(n)] for i in range(n)]
    return pow(z,n*n,p)*determinant_mod(H,p)*pow(A*A,-1,p)%p


def raw_matrix(N,d,z,p):
    normal=indices(d,N-1)
    rows=[(nu,j,s) for nu in normal for j in (0,1) for s in range(2*(N-sum(nu)))]
    cols=[(kap,h,a) for kap in normal for h in (0,1) for a in range(2*(N-sum(kap)))]
    assert len(rows)==len(cols)==4*comb(N+d,d+1)
    full_source={(nu,j,s) for nu in indices(d,N) for j in (0,1) for s in range(2*N+1)
                 for beta_m in (0,1) if 2*s+4*sum(nu)+4*N*beta_m<4*N}
    assert set(rows)==full_source
    assert all(h+2*a+Fraction(2,N)*sum(kap)<4*N for kap,h,a in cols)
    result=[]
    for nu,j,s in rows:
        row=[]
        for kap,h,a in cols:
            if any(k<v for k,v in zip(kap,nu)):
                row.append(0)
                continue
            factor=prod(comb(k,v) for k,v in zip(kap,nu))
            power=a+sum(kap)-sum(nu)
            entry=sum(comb(power,l)*pow(j*z,power-l,p)*pow(h,s-l,p)*
                      pow(factorial(s-l),-1,p) for l in range(min(s,power)+1))
            row.append(factor*entry%p)
        result.append(row)
    return normal,result


def main():
    primes=(1000003,1000033)
    for p in primes:
        assert all(p%k for k in range(2,int(p**0.5)+1))
    records=[]
    for N,d in ((2,1),(2,2),(2,3),(3,1),(3,2),(3,3),(4,1)):
        normal=indices(d,N-1)
        degree=sum(2*(N-sum(nu))**2 for nu in normal)
        assert degree==4*comb(N+d,d+2)+2*comb(N+d,d+1)
        results=[]
        for p in primes:
            for z in (1,2,7):
                normal,matrix=raw_matrix(N,d,z,p)
                actual=determinant_mod(matrix,p)
                expected=prod(hankel_period(2*(N-sum(nu)),z,p) for nu in normal)%p
                assert actual==expected
                results.append({'prime':p,'center':z,'determinant':actual})
        assert any(r['determinant'] for r in results)
        records.append({'N':N,'d':d,'rows':4*comb(N+d,d+1),'degree':degree,'checks':results})
    out={'status':'passed','packets':len(records),'independent_raw_modular_checks':sum(len(r['checks']) for r in records),
         'records':records,'scope':'Modular identity checks supplement the general independent proof; no numerical asymptotic inference.'}
    Path(__file__).with_name('normal_block_review_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2))


if __name__=='__main__':
    main()
