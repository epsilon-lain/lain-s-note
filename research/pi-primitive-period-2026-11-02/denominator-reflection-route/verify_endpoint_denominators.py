"""Exact certificates for endpoint content and least coefficient denominators.

Standard library only. The infinite-n assertions follow from the accompanying
proof, not the finite checks. All coefficient operations below are integral or
rational; no floating-point decision is used.
"""
from fractions import Fraction as F
from itertools import combinations
from math import factorial, comb, prod, gcd, lcm
from pathlib import Path
import json


def A(n):
    return prod(factorial(j) for j in range(n))


def trim(p):
    while len(p)>1 and p[-1]==0:
        p.pop()
    return p


def add(a,b,sign=1):
    out=a+[0]*max(0,len(b)-len(a))
    for i,x in enumerate(b):out[i]+=sign*x
    return trim(out)


def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y:out[i+j]+=x*y
    return trim(out)


def eval_poly(p,x):
    out=0
    for c in reversed(p):out=out*x+c
    return out


def det_poly(matrix):
    n=len(matrix)
    dp={0:[1]}
    for i in range(n):
        nxt={}
        for mask,value in dp.items():
            for j in range(n):
                if (mask>>j)&1:continue
                inversions=i-(mask&((1<<j)-1)).bit_count()
                key=mask|(1<<j)
                nxt[key]=add(nxt.get(key,[0]),mul(value,matrix[i][j]),(-1)**inversions)
        dp=nxt
    return dp[(1<<n)-1]


def det_integer(matrix):
    a=[r[:] for r in matrix]
    if not a:return 1
    previous,sign=1,1
    for k in range(len(a)-1):
        pivot=next((i for i in range(k,len(a)) if a[i][k]),None)
        if pivot is None:return 0
        if pivot!=k:a[pivot],a[k]=a[k],a[pivot];sign=-sign
        p=a[k][k]
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):
                value=a[i][j]*p-a[i][k]*a[k][j]
                assert value%previous==0
                a[i][j]=value//previous
            a[i][k]=0
        previous=p
    return sign*a[-1][-1]


def raw_scaled_det(n,Z):
    def entry(j,s,h,a):
        return sum(comb(a,k)*(factorial(s)//factorial(s-k))*h**(s-k)*(j*Z)**(a-k)
                   for k in range(min(a,s)+1))
    return det_integer([[entry(j,s,h,a) for h in range(2) for a in range(n)]
                        for j in range(2) for s in range(n)])


def moment_poly(k,endpoint=1):
    # t * integral exp(-v) (endpoint-t*v)^k dv.
    return [0]+[(-1)**j*comb(k,j)*factorial(j)*endpoint**(k-j) for j in range(k+1)]


def interpolate_poly_values(values):
    # Values are themselves polynomials; reconstruct the outer variable.
    count=len(values)
    width=max(map(len,values))
    output=[[F(0)]*width for _ in range(count)]
    for tdegree in range(width):
        dif=[F(v[tdegree] if tdegree<len(v) else 0) for v in values]
        basis=[F(1)]
        for k in range(count):
            for j,c in enumerate(basis):output[j][tdegree]+=c*dif[0]
            dif=[b-a for a,b in zip(dif,dif[1:])]
            nxt=[F(0)]*(len(basis)+1)
            for j,c in enumerate(basis):
                nxt[j]-=c*k/(k+1);nxt[j+1]+=c/(k+1)
            basis=nxt
    return [trim(p) for p in output]


def endpoint_split_checks():
    records=[]
    for n in range(2,7):
        values=[]
        for E in range(n+1):
            moments=[add([E*c for c in moment_poly(k,1)],moment_poly(k,0),-1)
                     for k in range(2*n-1)]
            values.append(det_poly([[moments[i+j] for j in range(n)] for i in range(n)]))
        split=interpolate_poly_values(values)
        for k,h in enumerate(split):
            ell=n-k;exponent=k*k+ell*ell
            factor=(-1)**ell*A(k)**2*A(ell)**2
            assert all(c==0 for c in h[:exponent])
            tail=[c/factor for c in h[exponent:]]
            assert tail and tail[0]==1
            assert all(c.denominator==1 for c in tail)
            linear=tail[1] if len(tail)>1 else 0
            assert linear==2*k*ell*(ell-k)
            records.append({'n':n,'k':k,'l':ell,'valuation':exponent,
                            'constant':str(tail[0]),'linear':str(linear),
                            'integer_cross_coefficients':True})
    return records


def partitions(size,cap=None):
    if not size:yield ();return
    if cap is None:cap=size
    for first in range(min(size,cap),0,-1):
        for tail in partitions(size-first,first):yield(first,)+tail


def schur_checks():
    count=0
    for k in range(1,6):
        for size in range(9):
            for lam in partitions(size):
                if len(lam)>k:continue
                full=lam+(0,)*(k-len(lam))
                aa=[full[i]+k-i-1 for i in range(k)]
                base=[k-i-1 for i in range(k)]
                numerator=det_integer([[factorial(a+j) for j in range(k)] for a in aa])
                denominator=det_integer([[factorial(a+j) for j in range(k)] for a in base])
                actual=F(numerator,denominator)
                dimension=F(1)
                for i,j in combinations(range(k),2):dimension*=F(aa[i]-aa[j],j-i)
                formula=dimension*prod(F(factorial(a),factorial(b)) for a,b in zip(aa,base))
                assert actual==formula and actual.denominator==1 and actual>0
                # Independent hook-content evaluation of the integral basis value.
                hook_dimension=F(1)
                for i,length in enumerate(lam):
                    for j in range(length):
                        hook=length-j+sum(x>j for x in lam[i+1:])
                        hook_dimension*=F(k+j-i,hook)
                assert hook_dimension==dimension and dimension.denominator==1
                count+=1
    return count


def denominator_checks():
    known={4:144,5:864,6:74649600,7:22394880000,
           8:758487711744000000,9:26759446470328320000000,
           10:489407583931698134384640000000000}
    records=[]
    for n in range(2,11):
        moments=[[0]]
        for k in range(1,2*n-1):moments.append(add([0,1],[0]+[-k*c for c in moments[-1]]))
        hankel=det_poly([[moments[i+j] for j in range(n)] for i in range(n)])
        integer_poly=list(reversed(hankel+[0]*(n*n+1-len(hankel))))
        integer_poly=trim(integer_poly)
        r,s=n//2,(n+1)//2
        rowden=A(n)**2;G=A(r)**2*A(s)**2
        assert rowden%G==0 and all(c%G==0 for c in integer_poly)
        common=rowden//G
        content=gcd(*integer_poly)
        least=rowden//gcd(rowden,content)
        polynomial=[F(c,rowden) for c in integer_poly]
        assert least==lcm(*(c.denominator for c in polynomial))
        assert polynomial[n]==(-1)**(n+1)*(n-1)
        if n%2==0:
            degree=2*r*r;leading=F((-1)**r,common)
            assert least==common
        else:
            degree=2*r*r+2*r-1
            c=4*r*(r+1);leading=F((-1)**(r+1)*c,common)
            lower=common//gcd(common,c)
            assert least%lower==0 and common%least==0
        assert len(polynomial)-1==degree and polynomial[-1]==leading
        for z in (-1,1,2):assert eval_poly(integer_poly,z)==raw_scaled_det(n,z)
        if n in known:assert least==known[n]
        records.append({'n':n,'degree':degree,'leading_coefficient':str(leading),
                        'valid_denominator':str(common),'least_denominator':str(least),
                        'raw_integer_content':str(content),'forced_content':str(G),
                        'raw_evaluations_checked':3})
    divisibilities=0
    for n in range(2,61):
        r,s=n//2,(n+1)//2;G=A(r)**2*A(s)**2
        assert A(n)**2%G==0
        for k in range(n+1):
            assert (A(k)**2*A(n-k)**2)%G==0
            divisibilities+=1
        if n%2 and n>=5:
            assert (A(n)**2//G)%(4*r*(r+1))==0
    return records,divisibilities


def main():
    split=endpoint_split_checks()
    schur=schur_checks()
    records,divisibilities=denominator_checks()
    out={'status':'passed','endpoint_split_checks':split,'laguerre_schur_integrality_checks':schur,
         'denominator_cases':records,'factorial_content_divisibility_checks':divisibilities,
         'scope':'Exact finite checks for the endpoint proof and coefficient-denominator bounds. '
                 'The all-n denominator asymptotic follows from the accompanying proof. '
                 'No rational approximation or pi bad-approximability conclusion.'}
    Path(__file__).with_name('endpoint_denominator_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','endpoint_split_checks':len(split),
                      'laguerre_schur_integrality_checks':schur,'denominator_cases':len(records),
                      'raw_determinant_evaluations':3*len(records),
                      'factorial_content_divisibility_checks':divisibilities},indent=2))


if __name__=='__main__':main()
