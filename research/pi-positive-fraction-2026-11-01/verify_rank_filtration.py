"""Exact finite checks for boundary rank filtration and determinant jets."""
from fractions import Fraction as F
from itertools import product
from math import comb, ceil, floor
from pathlib import Path
import json


def rank(a):
    if not a or not a[0]:return 0
    a=[[F(x) for x in row] for row in a];r=0
    for c in range(len(a[0])):
        p=next((i for i in range(r,len(a)) if a[i][c]),None)
        if p is None:continue
        a[r],a[p]=a[p],a[r];v=a[r][c]
        a[r]=[x/v for x in a[r]]
        for i in range(r+1,len(a)):
            q=a[i][c]
            if q:a[i]=[x-q*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r


def inverse(a):
    n=len(a);a=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(a)]
    for c in range(n):
        p=next(i for i in range(c,n) if a[i][c]);a[c],a[p]=a[p],a[c]
        v=a[c][c];a[c]=[x/v for x in a[c]]
        for i in range(n):
            if i!=c:
                q=a[i][c]
                if q:a[i]=[x-q*y for x,y in zip(a[i],a[c])]
    return [row[n:] for row in a]


def mul(p,q,D):
    return [sum(p[i]*q[k-i] for i in range(k+1)) for k in range(D+1)]


def det_series(a,D):
    n=len(a);dp={0:[F(1)]+[F(0)]*D}
    for row in range(n):
        nxt={}
        for mask,v in dp.items():
            for c in range(n):
                if mask>>c&1:continue
                target=mask|(1<<c)
                value=mul(v,a[row][c],D)
                sign=(-1)**((mask>>(c+1)).bit_count())
                if target not in nxt:nxt[target]=[F(0)]*(D+1)
                nxt[target]=[x+sign*y for x,y in zip(nxt[target],value)]
        dp=nxt
    return dp[(1<<n)-1]


def multi(m,N):
    return [x for x in product(range(N+1),repeat=m) if sum(x)<=N]


def residual(source,target,xi,v0,vi,H):
    s,b=source;t,a=target;g=tuple(x-y for x,y in zip(a,b))
    if min(g)<0 or sum(g)!=s-t:return 0,0
    c=tuple(v-v0 for v in vi);gap=H-v0*s-sum(v*x for v,x in zip(vi,b))
    if sum(x*y for x,y in zip(c,g))<gap:return 0,0
    alt=sum((-1)**sum(k)*__import__('math').prod(comb(x,y) for x,y in zip(g,k))
            for k in product(*(range(x+1) for x in g))
            if sum(x*y for x,y in zip(c,k))<gap)
    scalar=F(alt)
    for x,y,z,n in zip(a,b,xi,g):scalar*=comb(x,y)*(-z)**n
    return sum(g),scalar


def test_jets():
    records=[];support_checks=0
    for m,vi,H in ((1,(3,),4),(2,(2,3),4)):
        v0=1;xi=tuple(F((-1)**i,i+3) for i in range(m))
        ambient=[(s,b) for s in range(H) for b in multi(m,H-1-s)]
        inside=[r for r in ambient if r[0]+sum(v*x for v,x in zip(vi,r[1]))<H]
        outside=[r for r in ambient if r not in inside]
        n=len(inside)
        B=[[F(int(i==j)*(i+1)+1) for j in range(n)] for i in range(n)]
        Bi=inverse(B)
        O=[[F((i+2)**(j%4+1)+(j+1)**(i%3+1)) for j in range(n)] for i in range(len(outside))]
        for D in (1,2,3):
            E=[[[F(0)]*(D+1) for _ in range(n)] for _ in range(n)]
            active=[i for i,r in enumerate(inside) if H-r[0]-sum(v*x for v,x in zip(vi,r[1]))<=max(v-1 for v in vi)*D]
            for i,r in enumerate(inside):
                for q,t in enumerate(outside):
                    degree,value=residual(r,t,xi,v0,vi,H)
                    if value and degree<=D:
                        assert i in active
                        for j in range(n):E[i][j][degree]+=value*O[q][j]
                    support_checks+=1
            reduced=[[[sum(E[i][k][d]*Bi[k][j] for k in range(n)) for d in range(D+1)]
                      for j in active] for i in active]
            db=det_series([[[x]+[F(0)]*D for x in row] for row in B],D)[0]
            for eta in range(D+1):
                full=[[[B[i][j] if d==0 else eta*E[i][j][d] for d in range(D+1)]
                       for j in range(n)] for i in range(n)]
                small=[[[F(i==j) if d==0 else eta*reduced[i][j][d] for d in range(D+1)]
                        for j in range(len(active))] for i in range(len(active))]
                assert det_series(full,D)==[db*x for x in det_series(small,D)]
            records.append({'m':m,'H':H,'D':D,'source_rows':n,'strip_rows':len(active),'eta_values':D+1})
    return records,support_checks


def test_rank_blocks():
    count=0;cap_changes=0
    for v0,v1,theta in ((1,3,F(1,2)),(2,5,F(1,2)),(1,10,F(1,3)),(3,4,F(1,2))):
        c=v1-v0;w=theta*v1
        for H in range(1,25):
            A=floor(F(H)/w)
            for N in range(ceil(F(H,v0))):
                b=min(N+1,ceil(F(H-v0*N,c)))
                targets=list(range(b,min(N,A)+1))
                matrix=[[F(comb(a,beta)*(-1)**(b-beta-1)*comb(a-beta-1,b-beta-1))*F(-2,7)**(a-beta)
                         for a in targets] for beta in range(b)]
                expected=min(b,len(targets))
                assert rank(matrix)==expected
                uncapped=min(b,N+1-b)
                if H>=F(v1*c,v0):assert expected==uncapped
                cap_changes+=expected!=uncapped
                count+=1
    return count,cap_changes


def exact_integral_ratio(v0,v1):
    v0=F(v0);v1=F(v1);c=v1-v0
    a=1/v1;b=2/(v1+v0);d=1/v0
    first=(v1*(b*b-a*a)/2-(b-a))/c
    second=((d-b)-v0*(d*d-b*b)/2)/c
    return (first+second)*2*v0*v1


def test_counts():
    records=[]
    for v0,v1 in ((1,3),(2,5),(1,10),(3,4),(1,101)):
        limit=F(v1-v0,v1+v0)
        assert exact_integral_ratio(v0,v1)==limit
        for H in (100,400,1600):
            w=F(v1,2);cap=floor(F(H)/w);total=0;rank_total=0
            for N in range(ceil(F(H,v0))):
                b=min(N+1,ceil(F(H-v0*N,v1-v0)))
                total+=b;rank_total+=min(b,max(0,min(N,cap)-b+1))
            records.append({'v0':v0,'v1':v1,'H':H,'source_rows_per_center':total,'boundary_rank_per_center':rank_total,'ratio':str(F(rank_total,total)),'proved_limit':str(limit)})
    return records


if __name__=='__main__':
    jets,support=test_jets();blocks,cap=test_rank_blocks()
    result={'status':'passed','determinant_jet_reductions':jets,'boundary_support_checks':support,'exact_cauchy_rank_blocks':blocks,'small_height_cap_changes':cap,'rank_counts':test_counts(),'scope':'Exact finite rank and formal truncated determinant identities. Infinite ratios follow from the written piecewise-integral proof, not the displayed finite counts. Rank of coefficient transport is not automatically rank on compatible period columns; no norm or approximation theorem is certified.'}
    Path(__file__).with_name('rank_filtration_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:(len(v) if isinstance(v,list) else v) for k,v in result.items()},indent=2))
