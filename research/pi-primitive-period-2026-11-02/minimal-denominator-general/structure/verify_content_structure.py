"""Independent integer-Hankel content checks; Python standard library only.

The odd content formula is checked finitely, not proved by this script.
The even formula follows from the separately symmetric expectation argument.
"""
from fractions import Fraction as F
from math import factorial, prod, gcd
from pathlib import Path
import json


def trim(a):
    while len(a)>1 and a[-1]==0:
        a.pop()
    return a


def sub(a,b):
    c=a+[0]*max(0,len(b)-len(a))
    for i,v in enumerate(b):c[i]-=v
    return trim(c)


def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y:c[i+j]+=x*y
    return trim(c)


def exactdiv(a,b):
    a=a[:]
    if len(a)<len(b):
        assert a==[0]
        return [0]
    q=[0]*(len(a)-len(b)+1)
    while len(a)>=len(b) and a!=[0]:
        shift=len(a)-len(b)
        assert a[-1]%b[-1]==0
        c=a[-1]//b[-1]
        q[shift]=c
        for i,v in enumerate(b):a[i+shift]-=c*v
        trim(a)
    assert a==[0]
    return trim(q)


def A(n):return prod(factorial(j) for j in range(n))


def content_checks(N=12):
    f=[[0]]
    for k in range(1,2*N-1):
        p=[-k*x for x in f[-1]]+[0]*(k+1-len(f[-1]))
        p[k]+=1
        f.append(trim(p))
    older=[[1]]*(2*N+1)
    previous=f
    records=[]
    for n in range(2,N+1):
        current=[]
        for offset in range(2*(N-n)+1):
            top=sub(mul(previous[offset],previous[offset+2]),
                    mul(previous[offset+1],previous[offset+1]))
            current.append(exactdiv(top,older[offset+2]))
        p=current[0]
        g=gcd(*p)
        r=n//2
        basic=A(r)**2*A(n-r)**2
        predicted=basic if n%2==0 else 4*r*(r+1)*basic
        assert g==abs(p[-1])==predicted
        assert all(c%basic==0 for c in p)
        if n%2==0:
            assert len(p)-1==n*n//2
            assert p[-1]==(-1)**r*basic
            assert A(n)**2%basic==0
        else:
            assert len(p)-1==(n*n-3)//2
            assert p[-1]==(-1)**(r+1)*predicted
        records.append({'n':n,'degree':len(p)-1,'content':str(g),
                        'leading_coefficient':str(p[-1]),
                        'balanced_divisor':str(basic),
                        'least_integer_coefficient_denominator':str(A(n)**2//gcd(A(n)**2,g)),
                        'odd_content_equality_scope':'Finite check only' if n%2 else 'Matches proved even formula'})
        older,previous=previous,current
    return records


def partitions(n,cap=None):
    if n==0:
        yield ()
        return
    if cap is None:cap=n
    for k in range(min(n,cap),0,-1):
        for tail in partitions(n-k,k):yield (k,)+tail


def det(a):
    a=[list(map(F,row)) for row in a]
    value=F(1)
    for i in range(len(a)):
        pivot=next((j for j in range(i,len(a)) if a[j][i]),None)
        if pivot is None:return F(0)
        if pivot!=i:
            a[i],a[pivot]=a[pivot],a[i]
            value=-value
        v=a[i][i]
        value*=v
        for j in range(i+1,len(a)):
            c=a[j][i]/v
            for k in range(i+1,len(a)):a[j][k]-=c*a[i][k]
    return value


def schur_checks():
    count=0
    for k in range(1,7):
        for degree in range(8):
            for lam in partitions(degree):
                if len(lam)>k:continue
                padded=lam+(0,)*(k-len(lam))
                exponents=[i+padded[k-1-i] for i in range(k)]
                vander=prod(exponents[j]-exponents[i] for i in range(k) for j in range(i+1,k))
                ratio=F(vander,A(k))
                assert ratio.denominator==1
                formula=ratio*prod(factorial(a)//factorial(i) for i,a in enumerate(exponents))
                direct=det([[factorial(a+j) for j in range(k)] for a in exponents])/A(k)**2
                assert direct==formula and direct.denominator==1
                count+=1
    return count


def factorial_divisibility_checks():
    count=0
    for n in range(2,41):
        r=n//2
        s=n-r
        basic=A(r)**2*A(s)**2
        assert A(n)**2%basic==0
        for k in range(n+1):
            assert (A(k)**2*A(n-k)**2)%basic==0
            count+=1
        if n%2 and r>=2:
            B=A(n)//(A(r)*A(s))
            explicit=factorial(r)*factorial(r+1)*prod(
                factorial(r+j)//factorial(j) for j in range(2,r+1))
            assert B==explicit
            assert B%(2*r*(r+1))==0
            assert (A(n)**2//basic)%(4*r*(r+1))==0
    return count


def main():
    records=content_checks()
    out={'status':'passed','content_records':records,
         'independent_factorial_determinant_schur_checks':schur_checks(),
         'balanced_factorial_divisibility_checks':factorial_divisibility_checks(),
         'scope':'Finite independent verification of the even theorem and balanced divisor. '
                 'Odd exact content agreement is experimental only.'}
    Path(__file__).with_name('content_structure_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','complete_polynomial_content_checks':len(records),
                      'schur_integrality_checks':out['independent_factorial_determinant_schur_checks'],
                      'balanced_factorial_divisibility_checks':out['balanced_factorial_divisibility_checks']},indent=2))


if __name__=='__main__':main()
