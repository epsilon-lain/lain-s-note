"""Exact coalescence valuation and jet-rank checks; Python standard library."""
from fractions import Fraction as Q
from math import comb, factorial
from random import Random
from itertools import product
from pathlib import Path
import json


def determinant(matrix):
    a=[[Q(x) for x in row] for row in matrix]
    value=Q(1)
    for k in range(len(a)):
        pivot=next((i for i in range(k,len(a)) if a[i][k]),None)
        if pivot is None:return Q(0)
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];value=-value
        p=a[k][k];value*=p
        for i in range(k+1,len(a)):
            c=a[i][k]/p
            for j in range(k+1,len(a)):a[i][j]-=c*a[k][j]
    return value


def rank(matrix):
    a=[[Q(x) for x in row] for row in matrix]
    row=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(row,len(a)) if a[i][col]),None)
        if pivot is None:continue
        a[row],a[pivot]=a[pivot],a[row]
        for i in range(row+1,len(a)):
            c=a[i][col]/a[row][col]
            for j in range(col,len(a[0])):a[i][j]-=c*a[row][j]
        row+=1
        if row==len(a):break
    return row


def multiply(a,b,N):
    out=[Q(0)]*(N+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<=N:out[i+j]+=x*y
    return out


def entry(s,beta,h,alpha,center):
    if alpha<beta:return Q(0)
    g=[Q(center)]+[Q((-1)**(i+1),i) for i in range(1,s+1)]
    value=[Q(1)]+[Q(0)]*s
    for _ in range(alpha-beta):value=multiply(value,g,s)
    return comb(alpha,beta)*sum(Q(comb(h,k))*value[s-k] for k in range(min(s,h)+1))


def rows(profile,K):
    return [(s,j,b) for s,B in enumerate(profile) for j in range(K) for b in range(B)]


def matrix(profile,K,cols,r):
    return [[entry(s,b,h,a,j*r) for h,a in cols] for s,j,b in rows(profile,K)]


def collapsed(profile,K,cols):
    return [[entry(s,a,h,alpha,0) for h,alpha in cols]
            for s,B in enumerate(profile) for a in range(K*B)]


def polynomial_from_values(values):
    """Newton forward differences at 0,1,...,d, converted exactly."""
    differences=list(map(Q,values));out=[Q(0)]*len(values);basis=[Q(1)]
    for k in range(len(values)):
        for i,b in enumerate(basis):out[i]+=differences[0]*b/factorial(k)
        differences=[differences[i+1]-differences[i] for i in range(len(differences)-1)]
        nxt=[Q(0)]*(len(basis)+1)
        for i,b in enumerate(basis):nxt[i]-=k*b;nxt[i+1]+=b
        basis=nxt
    while len(out)>1 and not out[-1]:out.pop()
    return out


def evaluate(poly,x):
    return sum(a*x**i for i,a in enumerate(poly))


def check_case(profile,K,cols):
    rr=rows(profile,K);assert len(rr)==len(cols)
    degree=sum(a for h,a in cols)-sum(b for s,j,b in rr)
    assert degree>=0
    pol=polynomial_from_values([determinant(matrix(profile,K,cols,r)) for r in range(degree+1)])
    W=K*(K-1)//2*sum(B*B for B in profile)
    assert all(not c for c in pol[:W])
    constant=1
    for B in profile:
        for i in range(K):
            for j in range(i+1,K):constant*= (j-i)**(B*B)
    lead=pol[W] if W<len(pol) else Q(0)
    assert lead==constant*determinant(collapsed(profile,K,cols))
    for x in [Q(-3,2),Q(2,3),Q(degree+2)]:
        assert evaluate(pol,x)==determinant(matrix(profile,K,cols,x))
    val=next((i for i,c in enumerate(pol) if c),None)
    return {'profile':profile,'K':K,'rows':len(rr),'universal_lower_order':W,
            'actual_order':val,'polynomial_coefficients':[str(c) for c in pol]},pol


def toeplitz_rank_checks():
    # Actual full row packet: H=4,v0=2,w1=w0=1,theta=1/4,K=2.
    # Rows beta=0,s=0,1; columns 1,X,X^3,Y X^2.
    coeff=[
        [[1,0,0,0],[1,0,0,0],[0,1,0,0],[0,1,0,0]],
        [[0,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,2]],
        [[0,0,0,0],[0,0,0,1],[0,0,0,0],[0,0,3,1]],
        [[0,0,0,0],[0,0,1,0],[0,0,0,0],[0,0,0,0]],
    ]
    for r in [Q(0),Q(1),Q(-2),Q(3,2)]:
        expected=matrix([1,1],2,[(0,0),(0,1),(0,3),(1,2)],r)
        assert expected==[[sum(coeff[k][i][j]*r**k for k in range(4)) for j in range(4)] for i in range(4)]
    ranks=[]
    for t in range(5):
        T=[[Q(0)]*(4*(t+1)) for _ in range(4*(t+1))]
        for block_i in range(t+1):
            for block_j in range(block_i+1):
                k=block_i-block_j
                if k>=len(coeff):continue
                for i in range(4):
                    for j in range(4):T[4*block_i+i][4*block_j+j]=coeff[k][i][j]
        ranks.append(rank(T))
    assert ranks==[2,5,8,12,16]
    deficiencies=[4*(t+1)-r for t,r in enumerate(ranks)]
    assert deficiencies==[2,3,4,4,4]
    return {'ranks':ranks,'deficiencies':deficiencies,'local_invariant_orders':[0,0,1,3]}


def genuine_cancellation():
    # A full source-format row packet, without claiming its unknown weight thresholds.
    H=6;v0=2;w0=w1=1;theta=Q(1,6);K=2
    rr=[(s,j,b) for s in range(4) for j in range(K) for b in range(2)
        if v0*s+w1*b/theta<H]
    assert rr==rows([1,1,1],2)
    cols=[(0,1),(1,0),(2,1),(2,2),(2,4),(3,2)]
    assert all(w0*h+w1*a<=H for h,a in cols)
    case,p=check_case([1,1,1],K,cols)
    assert p==[Q(0)]*7+[Q(5)]
    terms=[];grouped={}
    for aa in product(range(5),repeat=3):
        A=[[entry(s,a,h,alpha,0) for h,alpha in cols]
           for s in range(3) for a in [0,aa[s]]]
        value=determinant(A)
        if not value:continue
        degree=sum(aa);grouped[degree]=grouped.get(degree,Q(0))+value
        terms.append({'types':aa,'degree':degree,'coefficient':str(value)})
    assert len(terms)==11
    assert grouped=={6:Q(0),7:Q(5),8:Q(0)}
    assert [t['coefficient'] for t in terms if t['degree']==6]==['-12','16','-4']
    return {'case':case,'columns':cols,'nonzero_terms':terms,
            'grouped_coefficients':{str(k):str(v) for k,v in grouped.items()},
            'smallest_individually_nonzero_degree':6,'summed_order':7}


def main():
    cases=[]
    explicit,pol=check_case([1,1],2,[(0,0),(0,1),(0,3),(1,2)])
    assert pol==[Q(0),Q(0),Q(0),Q(0),Q(1),Q(-1)]
    cases.append(explicit)
    for profile,K in [([2,1],2),([2,1],3),([1,1,1],2)]:
        cols=[(s,a) for s,B in enumerate(profile) for a in range(K*B)]
        case,p=check_case(profile,K,cols);assert case['actual_order']==case['universal_lower_order']
        assert sum(bool(c) for c in p)==1;cases.append(case)
    rng=Random(20261101);pool=[(h,a) for h in range(3) for a in range(5)]
    nonzero=zero=0
    for _ in range(32):
        cols=rng.sample(pool,4)
        case,p=check_case([1,1],2,cols)
        if case['actual_order'] is None:zero+=1
        else:nonzero+=1
    out={'status':'passed','exact_cases':cases,'additional_minor_checks':32,
         'additional_nonzero_minors':nonzero,'additional_zero_minors':zero,
         'rank_filtration':toeplitz_rank_checks(),
         'genuine_cancellation':genuine_cancellation(),
         'scope':'Universal one-dimensional coalescence divisibility; not small-error vanishing at omega, and no pi approximation conclusion.'}
    (Path(__file__).parent/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
