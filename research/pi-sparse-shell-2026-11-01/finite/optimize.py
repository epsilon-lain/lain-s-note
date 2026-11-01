"""Exact finite convex-cost flow and anchor Schur projection experiments."""
from fractions import Fraction as F
from collections import Counter, deque
from itertools import product
from math import comb
from pathlib import Path
import json, random, sys, time
import matrix_core as mc

P=1_000_000_009
mc.PRIME=P
DEN=10**12
CNUM=173_286_795_140

def ceilf(x):
    x=F(x); return -((-x.numerator)//x.denominator)

def rows(H,weights,K,theta,v0):
    return [(j,s,b) for j in range(K)
            for b in product(*(range(int(H*theta/w)+1) for w in weights))
            for s in range(max(0,ceilf((H-sum(w*x for w,x in zip(weights,b))/theta)/v0)))]

def weight(a,weights): return sum(x*w for x,w in zip(a,weights))

def cbounds(N=36):
    lo=sum((F(1,2*(2*k+1)*3**(2*k+1)) for k in range(N)),F())
    # ln(2)/4 = (1/2) atanh(1/3); positive tail bounded geometrically.
    hi=lo+F(9,16*(2*N+1)*3**(2*N+1))
    return lo,hi

def positive_pair(pair,allow_zero=True):
    a,b=pair;lo,hi=cbounds()
    if a==b==0:return allow_zero
    return a+b*(lo if b>=0 else hi)>0

class Flow:
    def __init__(self):self.names=[];self.adj=[];self.edges=[]
    def node(self,name):
        self.names.append(name);self.adj.append([]);return len(self.names)-1
    def edge(self,u,v,cap,a=0,b=0):
        i=len(self.edges)
        self.edges.extend([[u,v,cap,a,b],[v,u,0,-a,-b]])
        self.adj[u].append(i);self.adj[v].append(i+1)
        return i
    def dist(self,start=None):
        n=len(self.names); big=10**100
        ds=[0]*n if start is None else [big]*n
        pairs=[(0,0)]*n;prev=[None]*n
        queue=deque(range(n)) if start is None else deque([start])
        inside=[True]*n if start is None else [False]*n
        if start is not None:ds[start]=0;inside[start]=True
        while queue:
            u=queue.popleft();inside[u]=False
            for i in self.adj[u]:
                _,v,cap,a,b=self.edges[i]
                if not cap:continue
                nd=ds[u]+a*DEN+b*CNUM
                if nd<ds[v]:
                    ds[v]=nd;prev[v]=i;pairs[v]=(pairs[u][0]+a,pairs[u][1]+b)
                    if not inside[v]:queue.append(v);inside[v]=True
        return prev,pairs
    def send(self,s,t,count):
        for _ in range(count):
            prev,_=self.dist(s)
            assert prev[t] is not None
            v=t
            while v!=s:
                i=prev[v];self.edges[i][2]-=1;self.edges[i^1][2]+=1;v=self.edges[i][0]
    def dual(self):
        _,p=self.dist()
        for u,v,cap,a,b in self.edges:
            if cap:assert positive_pair((a+p[u][0]-p[v][0],b+p[u][1]-p[v][1]))
        return p

def make_flow(rr,H,weights,seed=0):
    ff=Flow();source=ff.node(['source']);sink=ff.node(['sink'])
    types=[a for a in product(*(range(H//w+1) for w in weights)) if weight(a,weights)<=H]
    fixed=Counter(b for j,s,b in rr if j==0)
    free=[i for i,(j,s,b) in enumerate(rr) if j!=0]
    slots={};tnodes={a:ff.node(['type',list(a)]) for a in types}
    order=list(free);random.Random(seed).shuffle(order)
    choices={}
    for i in order:
        j,s,beta=rr[i];r=ff.node(['row',i]);ff.edge(source,r,1)
        options=[a for a in types if all(x>=b for x,b in zip(a,beta))]
        random.Random(seed*103+i).shuffle(options)
        for a in options:
            key=(j,s,a)
            if key not in slots:
                slots[key]=ff.node(['slot',j,s,list(a)])
                ff.edge(slots[key],tnodes[a],1)
            e=ff.edge(r,slots[key],1,2*(weight(a,weights)-weight(beta,weights)))
            choices[i,a]=e
    for a in types:
        for k in range(1,len(free)+1):ff.edge(tnodes[a],sink,1,0,2*fixed[a]+2*k-1)
    ff.send(source,sink,len(free));dual=ff.dual()
    assignment=[b for j,s,b in rr]
    for (i,a),e in choices.items():
        if ff.edges[e][2]==0:assignment[i]=a
    assert len({(rr[i][0],rr[i][1],a) for i,a in enumerate(assignment)})==len(rr)
    return ff,dual,assignment

def eliminate(matrix):
    aa=[r[:] for r in matrix];piv=[];det=1;r=0
    for c in range(len(aa[0])):
        z=next((i for i in range(r,len(aa)) if aa[i][c]),None)
        if z is None:continue
        if z!=r:aa[z],aa[r]=aa[r],aa[z];det=-det
        v=aa[r][c];det=det*v%P;iv=pow(v,-1,P)
        aa[r]=[x*iv%P for x in aa[r]]
        for k in range(r+1,len(aa)):
            f=aa[k][c]
            if f:aa[k]=[(x-f*y)%P for x,y in zip(aa[k],aa[r])]
        piv.append(c);r+=1
        if r==len(aa):break
    return r,piv,det%P,aa[:r]

def project(matrix,anchor_echelon,anchor_pivots):
    freecols=[i for i in range(len(matrix[0])) if i not in set(anchor_pivots)]
    out=[]
    for line in matrix:
        v=line[:]
        for row,piv in zip(anchor_echelon,anchor_pivots):
            f=v[piv]
            if f:v=[(x-f*y)%P for x,y in zip(v,row)]
        assert not any(v[p] for p in anchor_pivots)
        out.append([v[c] for c in freecols])
    return out,freecols

def matrix_certificate(cc,rr,assignment):
    full=mc.vectors(cc,rr,assignment,W=2)
    anchor=[full[i] for i,r in enumerate(rr) if r[0]==0]
    other=[full[i] for i,r in enumerate(rr) if r[0]!=0]
    ar,ap,ad,ae=eliminate(anchor);assert ar==len(anchor)
    quotient,freecols=project(other,ae,ap)
    qr,qp,qd,_=eliminate(quotient)
    return {'anchor_rank':ar,'anchor_pivots':ap,'anchor_pivot_det':ad,
            'quotient_rank':qr,'quotient_minor_det':qd,
            'minor_indices':ap+[freecols[i] for i in qp],
            'det_in_this_column_order':ad*qd%P}

def main():
    fixtures=[
        {'H':8,'weights':[1,5],'theta':'1/4','rows':48,'columns':55,'cost_parts':['25/192','33/16']},
        {'H':16,'weights':[1,5],'theta':'1/3','rows':208,'columns':262,'cost_parts':['3/52','3831/1664']},
        {'H':12,'weights':[2,5],'theta':'1/4','rows':64,'columns':73,'cost_parts':['5/32','761/384']},
        {'H':8,'weights':[1,2,5],'theta':'1/3','rows':68,'columns':108,'cost_parts':['13/136','551/272']},
    ]
    for x in fixtures:x.update({'K':2,'v0':'1/2','w0':'1'})
    results=[]
    for old in fixtures:
        H=old['H'];ww=old['weights'];rr=rows(H,ww,old['K'],F(old['theta']),F(old['v0']))
        cc=mc.columns(len(ww),H,ww,F(old['w0']))
        for seed in range(20):
            t=time.time();ff,dual,aa=make_flow(rr,H,ww,seed)
            cost=mc.cost_parts(rr,aa,H,ww);cert=matrix_certificate(cc,rr,aa)
            print(H,ww,'seed',seed,'cost',list(map(str,cost)),'rank',cert['anchor_rank']+cert['quotient_rank'],'/',len(rr),'seconds',round(time.time()-t,2),flush=True)
            if cert['anchor_rank']+cert['quotient_rank']==len(rr):break
        result={k:v for k,v in old.items() if k in ('H','weights','K','theta','v0','w0','rows','columns')}
        result.update({'seed':seed,'assignment':aa,'cost_parts':list(map(str,cost)),
                       'old_cost_parts':old['cost_parts'],'certificate':cert,
                       'dual_potentials':[[name,list(pair)] for name,pair in zip(ff.names,dual)]})
        results.append(result)
        Path(__file__).with_name('optimal_fixtures.json').write_text(json.dumps(results,separators=(',',':'))+'\n')

if __name__=='__main__':main()
