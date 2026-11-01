"""Exact finite checks for lower-grid Pascal base packing.
The theorem proof is in BUCKET_BASE_PACKING.md. No pi claim is tested.
"""
from fractions import Fraction as F
from itertools import product, combinations
from math import comb, prod, ceil, floor
from pathlib import Path
import json, random


def below(a,b): return all(x<=y for x,y in zip(a,b))
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def binomial(a,b):
    return prod(comb(x,y) for x,y in zip(a,b)) if below(b,a) else 0


def rank_q(columns):
    pivots={}
    for column in columns:
        v=list(map(F,column))
        for i,p in sorted(pivots.items()):
            if v[i]:
                q=v[i];v=[x-q*y for x,y in zip(v,p)]
        lead=next((i for i,x in enumerate(v) if x),None)
        if lead is not None:
            q=v[lead];pivots[lead]=[x/q for x in v]
    return len(pivots)


def determinant(a):
    a=[row[:] for row in a]
    n=len(a);sign=1;previous=1
    if not n:return 1
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None:return 0
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];sign=-sign
        p=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value=a[i][j]*p-a[i][k]*a[k][j]
                assert value%previous==0
                a[i][j]=value//previous
            a[i][k]=0
        previous=p
    return sign*a[-1][-1]


def determinant_q(a):
    a=[list(map(F,row)) for row in a];out=F(1)
    for k in range(len(a)):
        pivot=next((i for i in range(k,len(a)) if a[i][k]),None)
        if pivot is None:return F(0)
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];out=-out
        p=a[k][k];out*=p
        for i in range(k+1,len(a)):
            q=a[i][k]/p
            for j in range(k+1,len(a)):a[i][j]-=q*a[k][j]
    return out


def greedy_mod(points,source,prime=1000003):
    pivots={};chosen=[]
    for a in points:
        v=[binomial(a,b)%prime for b in source]
        for i,p in sorted(pivots.items()):
            if v[i]:
                q=v[i];v=[(x-q*y)%prime for x,y in zip(v,p)]
        lead=next((i for i,x in enumerate(v) if x),None)
        if lead is not None:
            inverse=pow(v[lead],-1,prime)
            pivots[lead]=[(x*inverse)%prime for x in v]
            chosen.append(a)
            if len(chosen)==len(source):break
    return chosen


def pack(points,source,k):
    for seed in range(100):
        remaining=points[:];random.Random(seed).shuffle(remaining)
        groups=[]
        for _ in range(k):
            base=greedy_mod(remaining,source)
            if len(base)!=len(source):break
            groups.append(base);used=set(base)
            remaining=[a for a in remaining if a not in used]
        if len(groups)==k:return groups
    raise AssertionError('Search did not find a certificate; theorem not inferred from search.')


def matching(source,target):
    owner={}
    def visit(i,seen):
        for j,a in enumerate(target):
            if j in seen or not below(source[i],a):continue
            seen.add(j)
            if j not in owner or visit(owner[j],seen):
                owner[j]=i;return True
        return False
    for i in range(len(source)):assert visit(i,set())
    return [(source[i],target[j]) for j,i in sorted(owner.items())]


def small_checks():
    source=[(0,0),(1,0),(0,1)];t=F(3,2);points=set()
    for b in source:
        ranges=[range(floor(t*x),ceil(t*(x+1))) for x in b]
        points.update(product(*ranges))
    points=sorted(points)
    shadow_checks=0
    for mask in range(1<<len(source)):
        s=[b for i,b in enumerate(source) if mask>>i&1]
        shadow=[a for a in points if any(below(b,a) for b in s)]
        assert len(shadow)>=t**2*len(s)
        shadow_checks+=1
    rank_checks=0
    for mask in range(1<<len(points)):
        e=[a for i,a in enumerate(points) if mask>>i&1]
        rank=rank_q([[binomial(a,b) for b in source] for a in e])
        assert len(points)-len(e)>=2*(len(source)-rank)
        rank_checks+=1
    groups=pack(points,source,2)
    ds=[determinant([[binomial(a,b) for a in base] for b in source]) for base in groups]
    assert all(ds)
    return {'source':source,'target':points,'shadow_checks':shadow_checks,
            'rank_inequality_checks':rank_checks,'bases':groups,'determinants':ds}


def bucket_case(c,L,N,t,k):
    m=len(c);v0=F(1);theta=F(1,2);H=v0*N+L
    cap=t*(L+sum(c));weights=tuple(theta*(x+v0) for x in c)
    assert t**m>=k and N>=cap/min(c)
    source=[b for b in product(*(range(ceil(L/x)) for x in c))
            if dot(c,b)<L and sum(b)<=N]
    points=[a for a in product(*(range(ceil(cap/x)) for x in c))
            if dot(c,a)<cap and sum(a)<=N]
    groups=pack(points,source,k)
    assert len(set(a for base in groups for a in base))==k*len(source)
    ds=[];delta=[0]*m;permutations=[]
    for base in groups:
        d=determinant([[binomial(a,b) for a in base] for b in source])
        assert d
        ds.append(d)
        pairs=matching(source,base)
        permutations.append(pairs)
        for b,a in pairs:
            assert sum(a)<=N and N-sum(a)<=N-sum(b)
            assert dot(weights,a)<=H
            for i in range(m):delta[i]+=a[i]-b[i]
    # Independent coefficient determinant identity for rational shifts.
    zeta=tuple(F(1,7+4*i) for i in range(m))
    factored=prod(ds)*prod((-z)**d for z,d in zip(zeta,delta))
    assert factored
    # One smaller exact rational transport determinant per group.
    for base,d in zip(groups,ds):
        exponents=[sum(a[i] for a in base)-sum(b[i] for b in source) for i in range(m)]
        assert all(x>=0 for x in exponents)
        if len(source)<=8:
            matrix=[[binomial(a,b)*prod((-zeta[i])**(a[i]-b[i]) for i in range(m))
                     if below(b,a) else F(0) for a in base] for b in source]
            assert determinant_q(matrix)==d*prod((-z)**e for z,e in zip(zeta,exponents))
    return {'costs':list(c),'L':L,'N':N,'H':str(H),'dilation':str(t),'K':k,
            'source_count':len(source),'target_count':len(points),'bases':groups,
            'pascal_determinants':[str(x) for x in ds],
            'aggregate_increment':delta,'aggregate_weighted_increment':str(dot(weights,delta)),
            'all_target_types_within_column_budget':True,
            'nonzero_exact_coefficient_numerator':str(factored.numerator),
            'nonzero_exact_coefficient_denominator':str(factored.denominator)}


if __name__=='__main__':
    result={'status':'passed','small_exhaustive':small_checks(),
            'bucket_certificates':[bucket_case((2,5),20,21,F(3,2),2),
                                   bucket_case((3,7),13,16,F(2),3)],
            'scope':'Lower-grid theorem and local conserved-order bucket certificates only. Global low-order gluing and nonzero complete period minors are not established.'}
    Path(__file__).with_name('bucket_packing_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'shadow_checks':result['small_exhaustive']['shadow_checks'],
                      'rank_checks':result['small_exhaustive']['rank_inequality_checks'],
                      'buckets':[{'rows':x['source_count'],'targets':x['target_count'],'K':x['K']} for x in result['bucket_certificates']]}))
