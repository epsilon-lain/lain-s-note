"""Finite certificate for bounded integer mixing of scaled polynomial columns.

All decisions use a prime field. Nonzero residues certify characteristic-zero
nonvanishing; a formal-period evaluation at W=1 certifies a nonzero polynomial.
The general existence and height statements are proved in the accompanying note.
"""
from math import comb, isqrt
from random import Random
from pathlib import Path
import json


def determinant(matrix, p):
    a=[row[:] for row in matrix]
    result=1
    for c in range(len(a)):
        pivot=next((r for r in range(c,len(a)) if a[r][c]%p),None)
        if pivot is None:
            return 0
        if pivot!=c:
            a[pivot],a[c]=a[c],a[pivot]
            result=-result
        value=a[c][c]%p
        result=result*value%p
        inv=pow(value,-1,p)
        for r in range(c+1,len(a)):
            factor=a[r][c]*inv%p
            if factor:
                for j in range(c+1,len(a)):
                    a[r][j]=(a[r][j]-factor*a[c][j])%p
                a[r][c]=0
    return result%p


def multiply(A,C,p):
    return [[sum(x*y for x,y in zip(row,col))%p for col in zip(*C)]
            for row in A]


def polynomial_product(a,b,order,p):
    result=[0]*order
    for i,x in enumerate(a):
        for j,y in enumerate(b[:order-i]):
            result[i+j]=(result[i+j]+x*y)%p
    return result


def jet_matrix(rows,columns,center_step,q,p,order):
    logseries=[0]+[pow(k,-1,p)*(-1)**(k+1)%p for k in range(1,order)]
    powers={}
    for j in {r[0] for r in rows}:
        series=logseries[:]
        series[0]=j*center_step%p
        powers[j,0]=[1]+[0]*(order-1)
        for k in range(1,max(k for h,k in columns)+1):
            powers[j,k]=polynomial_product(powers[j,k-1],series,order,p)
    one_plus={h:[comb(h,s)%p if s<=h else 0 for s in range(order)]
              for h,k in columns}
    products={(j,h,k):polynomial_product(one_plus[h],powers[j,k],order,p)
              for j in {r[0] for r in rows} for h,k in columns}
    return [[0 if k<beta else comb(k,beta)*pow(q,k,p)*products[j,h,k-beta][s]%p
             for h,k in columns] for j,s,beta in rows]


def main():
    p=10009
    assert all(p%d for d in range(2,isqrt(p)+1)) and p%4==1
    nonresidue=next(a for a in range(2,p) if pow(a,(p-1)//2,p)==p-1)
    i=pow(nonresidue,(p-1)//4,p)
    assert i*i%p==p-1
    H,w1,theta_den,K,q=12,2,4,2,7
    columns=[(h,k) for k in range(H//w1+1) for h in range(H-w1*k+1)]
    rows=[(j,s,beta) for j in range(K) for beta in range(H//(w1*theta_den)+1)
          for s in range(H-w1*theta_den*beta)]
    N,M=len(columns),len(rows)
    # Original rational Gaussian center step 2*i*(22/7).
    A=jet_matrix(rows,columns,44*i*pow(7,-1,p)%p,q,p,H)
    # Formal exact-period matrix evaluated at W=1, after e^(hjW)=1 reduction.
    B=jet_matrix(rows,columns,1,q,p,H)
    random=Random(3112026)
    for attempt in range(1,21):
        C=[[random.randrange(3) for _ in range(M)] for _ in range(N)]
        da=determinant(multiply(A,C,p),p)
        db=determinant(multiply(B,C,p),p)
        if da and db:
            break
    assert da and db
    assert all(0<=v<=2 for row in C for v in row)
    result={'status':'passed','prime':p,'sqrt_minus_one_mod_prime':i,
            'parameters':{'H':H,'m':1,'K':K,'w0':1,'v0':1,'w1':w1,'theta':'1/4',
                          'rational_approximation':'22/7','log_truncation_T1':18},
            'rows':[list(r) for r in rows],'columns':[list(c) for c in columns],
            'matrix_shape':[M,N],'mixing_shape':[N,M],'mixing_entries':C,
            'original_mixed_determinant_mod_prime':da,
            'period_mixed_determinant_at_W1_mod_prime':db,
            'attempts':attempt,
            'scope':'Both determinants are nonzero in characteristic zero. The second certifies a nonzero Q[W] polynomial; its value at 2*pi*i is nonzero using known transcendence of pi. No approximation-exponent or badly-approximable conclusion.'}
    Path(__file__).with_name('mixed_columns_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','matrix_shape','mixing_shape','attempts','original_mixed_determinant_mod_prime','period_mixed_determinant_at_W1_mod_prime']}))


if __name__=='__main__':
    main()
