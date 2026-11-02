"""Exact center-polynomial coalescence checks; standard library only."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, factorial, prod, ceil
from pathlib import Path
import json


@lru_cache(None)
def stirling(n,h):
    if n==h==0:return 1
    if n==0 or h==0 or h>n:return 0
    return h*stirling(n-1,h)+stirling(n-1,h-1)


def rows(H,K,v0,v):
    ans=[]
    for beta in product(*(range(ceil(H/x)) for x in v)):
        rem=H-sum(x*b for x,b in zip(v,beta))
        if rem<=0:continue
        for s in range(ceil(rem/v0)):
            for j in range(K):ans.append((j,s,beta))
    return sorted(ans,key=lambda r:(r[1],r[0],r[2]))


def columns(H,w0,w):
    ans=[]
    for alpha in product(*(range(int(H/x)+1) for x in w)):
        rem=H-sum(x*a for x,a in zip(w,alpha))
        if rem<0:continue
        ans.extend((h,alpha) for h in range(int(rem/w0)+1))
    return sorted(ans,key=lambda c:(sum(c[1]),c[1],c[0]))


def entry(row,col,Z,xi,shifted):
    j,s,beta=row;h,alpha=col
    if any(a<b for a,b in zip(alpha,beta)):return 0
    out=0
    for d in product(*(range(a-b+1) for a,b in zip(alpha,beta))):
        n=s-sum(d)
        if n<0:continue
        time=factorial(h)*stirling(n,h) if shifted else h**n
        value=(factorial(s)//factorial(n))*time
        for a,b,dd,x in zip(alpha,beta,d,xi):
            value*=comb(a,b)*comb(a-b,dd)*(j*Z*x)**(a-b-dd)
        out+=value
    return out


def det(a):
    a=[r[:] for r in a];n=len(a);sgn=1;prev=1
    for k in range(n-1):
        p=next((i for i in range(k,n) if a[i][k]),None)
        if p is None:return 0
        if p!=k:a[k],a[p]=a[p],a[k];sgn=-sgn
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                num=a[i][j]*pivot-a[i][k]*a[k][j]
                assert num%prev==0
                a[i][j]=num//prev
            a[i][k]=0
        prev=pivot
    return sgn*a[-1][-1]


def pivots(a):
    a=[list(map(F,r)) for r in a];out=[];k=0
    for c in range(len(a[0])):
        p=next((i for i in range(k,len(a)) if a[i][c]),None)
        if p is None:continue
        a[k],a[p]=a[p],a[k]
        for i in range(k+1,len(a)):
            r=a[i][c]/a[k][c]
            if r:
                for j in range(c,len(a[0])):a[i][j]-=r*a[k][j]
        out.append(c);k+=1
        if k==len(a):break
    return out


def interpolate(values):
    out=[F(0)]*len(values);basis=[F(1)];diff=list(map(F,values))
    for k in range(len(values)):
        for i,x in enumerate(basis):out[i]+=x*diff[0]
        diff=[b-a for a,b in zip(diff,diff[1:])]
        nxt=[F(0)]*(len(basis)+1)
        for i,x in enumerate(basis):nxt[i]-=x*k/(k+1);nxt[i+1]+=x/(k+1)
        basis=nxt
    while len(out)>1 and not out[-1]:out.pop()
    return out


def evaluate(poly,x):
    val=F(0)
    for c in reversed(poly):val=val*x+c
    return val


def block_counts(H,v0,v):
    ans=[]
    for nu in product(*(range(ceil(H/x)) for x in v[:-1])):
        rem=H-sum(x*a for x,a in zip(v[:-1],nu))
        if rem<=0:continue
        for s in range(ceil(rem/v0)):
            n=ceil((rem-v0*s)/v[-1])
            if n>0:ans.append(n)
    return ans


def check(H,K,v0,w0,w,theta,xi,shifted,staircase=False):
    v=[x/theta for x in w];rr=rows(H,K,v0,v)
    if staircase:
        assert len(w)==1
        cc=[(s,(a,)) for s in range(ceil(H/v0))
            for a in range(K*ceil((H-v0*s)/v[0]))]
    else:
        pool=columns(H,w0,w)
        pp=pivots([[entry(r,c,1,xi,shifted) for c in pool] for r in rr])
        assert len(pp)==len(rr)
        cc=[pool[i] for i in pp]
    assert all(w0*h+sum(x*a for x,a in zip(w,alpha))<=H for h,alpha in cc)
    degree=sum(sum(a) for h,a in cc)-sum(sum(b) for j,s,b in rr)
    denom=prod(factorial(s) for j,s,b in rr)
    def exact(Z):return F(det([[entry(r,c,Z,xi,shifted) for c in cc] for r in rr]),denom)
    poly=interpolate([exact(i) for i in range(degree+1)])
    assert any(poly)
    assert evaluate(poly,-1)==exact(-1)
    assert evaluate(poly,degree+2)==exact(degree+2)
    ns=block_counts(H,v0,v)
    forced=comb(K,2)*sum(n*n for n in ns)
    assert len(rr)==K*sum(ns)
    val=next(i for i,x in enumerate(poly) if x)
    assert val>=forced
    if staircase:assert val==forced and len(poly)-1==forced
    E=len(rr)*H-sum(sum(x*b for x,b in zip(w,beta)) for j,s,beta in rr)
    assert min(w)*forced<=E
    return dict(m=len(w),K=K,H=H,rows=len(rr),basis='shifted' if shifted else 'raw',
                staircase=staircase,direction=xi,forced_order=forced,actual_order=val,
                degree=len(poly)-1,polynomial={str(i):str(x) for i,x in enumerate(poly) if x},
                selected_columns=cc)


def main():
    cases=[]
    configs=[(3,2,F(1),F(1,4),[F(1,2)],F(1,4),[1]),
             (3,2,F(1),F(1,4),[F(1,2),F(1)],F(1,4),[1,2]),
             (5,2,F(2),F(1,2),[F(1,2),F(1)],F(1,4),[2,3])]
    for cfg in configs:
        for shifted in (False,True):cases.append(check(*cfg,shifted))
    cases.append(check(3,2,F(1),F(1,4),[F(1,2)],F(1,4),[1],True,True))
    cases.append(check(4,3,F(1),F(1,4),[F(1,2)],F(1,8),[1],True,True))
    counts=[]
    for v in ([F(2)],[F(2),F(4)],[F(2),F(4),F(4)]):
        for H in (24,48,96):
            K=2;ns=block_counts(H,F(1),v);M=K*sum(ns);V=comb(K,2)*sum(n*n for n in ns)
            counts.append(dict(m=len(v),H=H,rows=M,forced_order=V,ratio=str(F(V,M*H)),
                               proved_limit=str(F(K-1,(len(v)+2)*v[-1]))))
    out=dict(status='passed',arithmetic='exact rational polynomial interpolation',
             complete_packets=cases,fixed_weight_count_checks=counts,
             scope='Center coalescence order and necessary degree costs. No lower bound at the actual period or pi approximation claim.')
    Path(__file__).with_name('coalescence_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','complete_packets':len(cases),
                      'sharp_staircases':2,'count_checks':len(counts),
                      'orders':[(x['forced_order'],x['actual_order']) for x in cases]},indent=2))


if __name__=='__main__':main()
