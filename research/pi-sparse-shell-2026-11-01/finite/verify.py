"""Independent certificate checks; no optimization or modular-rank-zero pruning."""
from fractions import Fraction as F
from collections import Counter
from itertools import product
from math import comb, factorial, isqrt
from functools import lru_cache
from pathlib import Path
import json
import matrix_core as mc

P=1_000_000_009
I=430_477_711
mc.PRIME=P

def weight(a,ww):return sum(x*w for x,w in zip(a,ww))

def rows(H,ww,K,theta,v0):
    rr=[]
    for j in range(K):
        for b in product(*(range(int(H*theta/w)+1) for w in ww)):
            bound=(H-weight(b,ww)/theta)/v0
            bound=F(bound)
            for s in range(max(0,-((-bound.numerator)//bound.denominator))):rr.append((j,s,b))
    return rr

def log2_quarter_interval(N=36):
    lo=sum((F(1,2*(2*k+1)*3**(2*k+1)) for k in range(N)),F())
    hi=lo+F(9,16*(2*N+1)*3**(2*N+1))
    return lo,hi

def nonnegative(a,b):
    if a==b==0:return True
    lo,hi=log2_quarter_interval()
    return a+b*(lo if b>=0 else hi)>0

def flow_certificate(x,rr,aa):
    """Construct every residual arc from witness; check exact affine dual."""
    H=x['H'];ww=x['weights'];free=[i for i,r in enumerate(rr) if r[0]!=0]
    types=[a for a in product(*(range(H//w+1) for w in ww)) if weight(a,ww)<=H]
    fixed=Counter(b for j,s,b in rr if j==0)
    used=Counter(aa[i] for i in free)
    def key(name):return json.dumps(name,separators=(',',':'))
    potentials={key(name):tuple(p) for name,p in x['dual_potentials']}
    assert len(potentials)==len(x['dual_potentials'])
    expected_nodes={key(['source']),key(['sink'])}
    expected_nodes.update(key(['type',list(a)]) for a in types)
    slots=set();checks=0;flow_balance=Counter();total_a=total_b=0
    def edge(uname,vname,flow,a,b):
        nonlocal checks,total_a,total_b
        assert flow in (0,1)
        u,v=key(uname),key(vname)
        pu,pv=potentials[u],potentials[v]
        if flow<1:
            assert nonnegative(a+pu[0]-pv[0],b+pu[1]-pv[1]);checks+=1
        if flow>0:
            assert nonnegative(-a+pv[0]-pu[0],-b+pv[1]-pu[1]);checks+=1
            flow_balance[u]-=flow;flow_balance[v]+=flow
            total_a+=a;total_b+=b
    for i in free:
        j,s,beta=rr[i];r=['row',i];expected_nodes.add(key(r))
        edge(['source'],r,1,0,0)
        for a in types:
            if not all(u>=v for u,v in zip(a,beta)):continue
            slot=['slot',j,s,list(a)];slots.add((j,s,a));expected_nodes.add(key(slot))
            edge(r,slot,int(aa[i]==a),2*(weight(a,ww)-weight(beta,ww)),0)
    chosen={(rr[i][0],rr[i][1],aa[i]) for i in free}
    assert len(chosen)==len(free)
    for j,s,a in slots:edge(['slot',j,s,list(a)],['type',list(a)],int((j,s,a) in chosen),0,0)
    for a in types:
        for k in range(1,len(free)+1):
            edge(['type',list(a)],['sink'],int(k<=used[a]),0,2*fixed[a]+2*k-1)
    assert set(potentials)==expected_nodes
    assert flow_balance.pop(key(['source']))==-len(free)
    assert flow_balance.pop(key(['sink']))==len(free)
    assert all(v==0 for v in flow_balance.values())
    counts=Counter(aa)
    total_b+=sum(n*n for n in fixed.values())
    assert total_b==sum(n*n for n in counts.values())
    cp=[F(total_a,len(rr)*H),F(total_b,len(rr)*H)]
    assert list(map(str,cp))==x['cost_parts']
    return {'dual_nodes':len(potentials),'checked_residual_arcs':checks,
            'unnormalized_cost':[total_a,total_b],'counts':[[list(a),n] for a,n in sorted(counts.items())]},cp

def aggregate_certificate(x,rr,aa):
    """A still weaker relaxation: ignore individual row feasibility entirely."""
    H=x['H'];ww=x['weights']
    fixed=Counter(b for j,s,b in rr if j==0)
    free=Counter(a for (j,s,b),a in zip(rr,aa) if j!=0)
    types=[a for a in product(*(range(H//w+1) for w in ww)) if weight(a,ww)<=H]
    last=[(2*weight(a,ww),2*fixed[a]+2*free[a]-1,a) for a in types if free[a]]
    nxt=[(2*weight(a,ww),2*fixed[a]+2*free[a]+1,a) for a in types]
    # Every used marginal is strictly smaller than every unused marginal.
    # Checking only last/next suffices because each type's increments increase.
    for a,b,ty in last:
        for aa,bb,ty2 in nxt:
            assert (aa-a,bb-b)!=(0,0) and nonnegative(aa-a,bb-b)
            # An elementary two-term lower bound already suffices for all cases.
            assert all(aa-a+(bb-b)*endpoint>0 for endpoint in (F(14,81),F(25,144)))
    # Pick extremal pairs numerically only for readable reporting, verify exactly.
    hi=max(last,key=lambda t:t[0]+0.17328679514*t[1])
    lo=min(nxt,key=lambda t:t[0]+0.17328679514*t[1])
    assert all(nonnegative(hi[0]-a,hi[1]-b) for a,b,ty in last)
    assert all(nonnegative(a-lo[0],b-lo[1]) for a,b,ty in nxt)
    return {'candidate_type_count':len(types),'last_selected_marginal':list(hi[:2]),
            'first_unselected_marginal':list(lo[:2]),'strict_gap':list((lo[0]-hi[0],lo[1]-hi[1])),
            'unique_optimal_aggregate_counts':True,'all_pairs_checked':len(last)*len(nxt)}

def determinant(matrix):
    a=[r[:] for r in matrix];d=1;n=len(a)
    assert all(len(r)==n for r in a)
    for c in range(n):
        k=next((i for i in range(c,n) if a[i][c]),None)
        if k is None:return 0
        if k!=c:a[k],a[c]=a[c],a[k];d=-d
        v=a[c][c];d=d*v%P
        for i in range(c+1,n):
            f=a[i][c]*pow(v,-1,P)%P
            if f:
                for j in range(c+1,n):a[i][j]=(a[i][j]-f*a[c][j])%P
    return d%P

def schur_certificate(full,anchor_count,cert):
    piv=cert['anchor_pivots'];ix=cert['minor_indices']
    assert len(piv)==anchor_count and len(set(ix))==len(full)==len(ix)
    assert ix[:anchor_count]==piv
    # Exact modular row elimination restricted to the stored anchor pivots.
    anchor=[r[:] for r in full[:anchor_count]];anchor_det=1
    for k,p in enumerate(piv):
        r=next(i for i in range(k,anchor_count) if anchor[i][p])
        if r!=k:anchor[r],anchor[k]=anchor[k],anchor[r];anchor_det=-anchor_det
        v=anchor[k][p];anchor_det=anchor_det*v%P
        anchor[k]=[z*pow(v,-1,P)%P for z in anchor[k]]
        for i in range(k+1,anchor_count):
            f=anchor[i][p]
            anchor[i]=[(x-f*y)%P for x,y in zip(anchor[i],anchor[k])]
    assert anchor_det==cert['anchor_pivot_det']!=0
    quotient=[]
    for r in full[anchor_count:]:
        r=r[:]
        for p,ar in zip(piv,anchor):
            f=r[p];r=[(x-f*y)%P for x,y in zip(r,ar)]
        assert not any(r[p] for p in piv)
        quotient.append([r[i] for i in ix[anchor_count:]])
    qd=determinant(quotient)
    assert qd==cert['quotient_minor_det']!=0
    direct=determinant([[r[i] for i in ix] for r in full])
    assert direct==cert['det_in_this_column_order']==anchor_det*qd%P!=0
    return direct

def original_matrix(cc,rr,centers):
    N=max(r[1] for r in rr);degree=max(sum(a) for h,a in cc)
    log=[0]+[(-1)**(k+1)*pow(k,-1,P)%P for k in range(1,N+1)]
    powers={}
    for j in {r[0] for r in rr}:
        for i,R in enumerate(centers):
            f=log[:];f[0]=j*R%P;powers[j,i]=mc.powers(f,degree,N)
    @lru_cache(None)
    def poly(j,d):
        f=[1]+[0]*N
        for i,k in enumerate(d):f=mc.mul(f,powers[j,i][k],N)
        return f
    matrix=[]
    for j,s,b in rr:
        line=[]
        for h,alpha in cc:
            if any(x>y for x,y in zip(b,alpha)):line.append(0);continue
            factor=1
            for a,bb in zip(alpha,b):factor=factor*comb(a,bb)%P
            f=poly(j,tuple(a-bb for a,bb in zip(alpha,b)))
            line.append(factor*sum(comb(h,k)*f[s-k] for k in range(min(h,s)+1))%P)
        matrix.append(line)
    return matrix

def atan_bounds(d,N):
    s=sum((F((-1)**k,(2*k+1)*d**(2*k+1)) for k in range(N)),F())
    t=F((-1)**N,(2*N+1)*d**(2*N+1))
    return min(s,s+t),max(s,s+t)

def exp_bounds(x,N=60):
    x=F(x);lo=sum((x**k/factorial(k) for k in range(N+1)),F())
    return lo,lo+x**(N+1)/factorial(N+1)/(1-x/(N+2))

def main():
    assert all(P%d for d in range(2,isqrt(P)+1))
    assert I*I%P==P-1
    lo,hi=log2_quarter_interval();assert F(14,81)<lo<hi<F(25,144)
    al,au=atan_bounds(5,30);bl,bu=atan_bounds(239,10)
    pil,piu=16*al-4*bu,16*au-4*bl
    approximants={1:(6,2),2:(22,7),5:(355,113)}
    for w,(p,q) in approximants.items():
        assert exp_bounds(w-1)[1]<q<exp_bounds(w)[0]
        assert max(abs(pil-F(p,q)),abs(piu-F(p,q)))<F(1,q*q)
    results=[]
    fixtures=json.loads(Path(__file__).with_name('optimal_fixtures.json').read_text())
    for x in fixtures:
        H=x['H'];ww=x['weights'];K=x['K'];theta=F(x['theta']);v0=F(x['v0']);w0=F(x['w0'])
        assert v0<w0 and K*theta**len(ww)<1 and K*w0/v0*theta**len(ww)<1
        rr=rows(H,ww,K,theta,v0);cc=mc.columns(len(ww),H,ww,w0)
        aa=list(map(tuple,x['assignment']));M=len(rr);anchor_count=sum(r[0]==0 for r in rr)
        assert len(rr)==x['rows'] and len(cc)==x['columns'] and len(aa)==M
        assert len(set(rr))==len(rr) and len(set(cc))==len(cc)
        assert all(v0*s+weight(beta,ww)/theta<H for j,s,beta in rr)
        assert all(w0*h+weight(alpha,ww)<=H for h,alpha in cc)
        assert all(all(a>=b for a,b in zip(ty,r[2])) and weight(ty,ww)<=H for r,ty in zip(rr,aa))
        assert all(ty==r[2] for r,ty in zip(rr,aa) if r[0]==0)
        fc,cost=flow_certificate(x,rr,aa)
        ag=aggregate_certificate(x,rr,aa)
        test=mc.vectors(cc,rr,aa,W=2)
        cert=x['certificate'];assert cert['anchor_rank']==anchor_count and cert['quotient_rank']==M-anchor_count
        td=schur_certificate(test,anchor_count,cert)
        centers=[2*I*p*pow(q,-1,P)%P for p,q in (approximants[w] for w in ww)]
        original=original_matrix(cc,rr,centers)
        assert original[:anchor_count]==test[:anchor_count]
        ix=cert['minor_indices'];od=determinant([[r[i] for i in ix] for r in original]);assert od
        # F0=17 gives Ti=34wi, so every required original Gi jet equals log.
        assert max(r[1] for r in rr)<34 and all(34*w>max(r[1] for r in rr) for w in ww)
        old=list(map(F,x['old_cost_parts']))
        assert nonnegative(old[0]-cost[0],old[1]-cost[1]) and cost!=old
        bbar=F(sum(weight(r[2],ww) for r in rr),M*H)
        assert nonnegative(1-bbar-cost[0],-cost[1])
        result={k:x[k] for k in ('H','weights','K','theta','v0','w0','rows','columns','cost_parts','old_cost_parts')}
        result.update({'flow':fc,'aggregate_optimality':ag,'anchor_det':cert['anchor_pivot_det'],
                       'quotient_det':cert['quotient_minor_det'],'translated_det':td,'original_same_minor_det':od,
                       'bbar':str(bbar),'global_optimality':True,'strict_improvement':True})
        print(result,flush=True);results.append(result)
    # Separate short explicit witness in SIMPLE_OPTIMUM.md; all column indices 0..47.
    H=8;ww=[1,5];rr=rows(H,ww,2,F(1,4),F(1,2));cc=mc.columns(2,H,ww,F(1))
    aa=[]
    for j,s,b in rr:
        if j==0:a=b
        elif b==(1,0):a=(2,0)
        elif s<=4:a=(0,0)
        elif s<=11:a=(1,0)
        elif s==12:a=(2,0)
        else:a=(3,0)
        aa.append(a)
    assert all(all(a>=b for a,b in zip(ty,r[2])) and weight(ty,ww)<=H for r,ty in zip(rr,aa))
    assert len(set((r[0],r[1],ty) for r,ty in zip(rr,aa)))==48
    simple_test=determinant([r[:48] for r in mc.vectors(cc,rr,aa,W=2)])
    centers=[2*I*p*pow(q,-1,P)%P for p,q in (approximants[w] for w in ww)]
    simple_orig=determinant([r[:48] for r in original_matrix(cc,rr,centers)])
    assert simple_test==146986659 and simple_orig==918208015
    assert list(map(str,mc.cost_parts(rr,aa,H,ww)))==['13/96','63/32']
    Path(__file__).with_name('verification.json').write_text(json.dumps({'status':'passed','prime':P,'sqrt_minus_one':I,'fixtures':results,
        'explicit_H8_witness':{'test_det':simple_test,'original_det':simple_orig,'columns':list(range(48))},
        'scope':'Exact finite optimum among all a>=beta, wa<=H with ell=s, d=k=0, anchor fixed. No symbolic-zero claim inferred from modular rank. Unknown separated-weight thresholds and asymptotic claims are not certified.'},indent=2)+'\n')

if __name__=='__main__':main()
