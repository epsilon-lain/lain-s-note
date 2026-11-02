"""Bounded period-aware selection in the exact rational/log kernel.

The five feature weights and all rational sample locations are fixed before
quality evaluation. Machin intervals are used only after candidate generation.
"""
from free_rational_basis import *
from time import monotonic

def eval_poly(coeffs,x):
    z=F(0)
    for a in reversed(coeffs):z=z*x+a
    return z

def saturated(m,K,d):
    c=4*m+1;mu=[]
    for j in range(d+1):
        aa=bb=F(0)
        for k in range(K+1):
            a,b=rational_moment(2*m+j+k,2*m+k,c+k)
            f=comb(c+k-1,k);aa+=f*a;bb+=f*b
        mu.append((aa,bb))
    rows=[]
    for col in range(2):
        v=[x[col] for x in mu];L=lcm(*(x.denominator for x in v))
        row=[int(x*L) for x in v];g=reduce(gcd,map(abs,row))
        rows.append([x//g for x in row])
    mat=DomainMatrix([[ZZ(x) for x in r] for r in rows],(2,d+1),ZZ)
    diag,left,right=smith_normal_decomp(mat)
    assert diag==left*mat*right
    rank=sum(diag.to_Matrix()[i,i]!=0 for i in range(2))
    kernel=right.to_Matrix()[:,rank:].T
    assert sp.Matrix(rows)*kernel.T==sp.zeros(2,d+1-rank)
    # Coefficient preconditioning changes the basis, not its saturation.
    kernel=kernel.lll()
    return mu,[[int(x) for x in r] for r in kernel.tolist()]

def run():
    start=monotonic();lo,hi=rational_pi();summaries=[];candidates=[]
    for m in range(1,4):
        c,K,d=4*m+1,2*m,3*m+1
        mu,basis=saturated(m,K,d);dimension=len(basis)
        count=max(dimension+2,c+2,9)
        points=[F(j,count+1) for j in range(1,count+1)]
        base=[0]*(2*m)+[(-1)**j*comb(2*m,j) for j in range(2*m+1)]
        pdata=[]
        for h in basis:
            n=mul(base,h);beta,alpha=pi_period(n,c)
            samples=[eval_poly(n,x)/(1+x*x)**c for x in points]
            pdata.append((samples,beta,alpha))
        weights=[0,8*m,12*m,16*m,20*m]
        outputs=[];zero=0
        for power in weights:
            lam=F(1,2**power)
            features=[s+[lam*al,lam*be] for s,be,al in pdata]
            den=lcm(*(x.denominator for row in features for x in row))
            embedded=[[int(x*den) for x in row]+h for row,h in zip(features,basis)]
            reduced=sp.Matrix(embedded).lll()
            for row in reduced.tolist():
                h=[int(x) for x in row[-(d+1):]]
                assert any(h)
                # Recompute every embedding coordinate from its guarded H.
                n=mul(base,h);beta,alpha=pi_period(n,c)
                fs=[eval_poly(n,x)/(1+x*x)**c for x in points]+[lam*alpha,lam*beta]
                assert [F(x) for x in row[:len(fs)]]==[x*den for x in fs]
                assert all(sum(F(x)*mu[j][k] for j,x in enumerate(h))==0 for k in range(2))
                if not alpha:zero+=1;continue
                content=reduce(gcd,map(abs,h));h=[x//content for x in h]
                beta/=content;alpha/=content
                pq=-4*beta/alpha;q=pq.denominator
                up=max(abs(lo-pq),abs(hi-pq))*q*q
                lower=max(F(0),lo-pq,pq-hi)*q*q
                assert up-lower<F(1,10**30)
                record=dict(m=m,K=K,d=d,c=c,weight_power=power,
                            H=[str(x) for x in h],alpha=str(alpha),beta=str(beta),
                            p=str(pq.numerator),q=str(q),Q_interval=[str(lower),str(up)],
                            displayed_Q=float(up))
                outputs.append(record);candidates.append(record)
        best=min(outputs,key=lambda x:F(x['Q_interval'][1])) if outputs else None
        summaries.append(dict(m=m,K=K,d=d,kernel_dimension=dimension,
                              rational_sample_count=count,five_weight_powers=weights,
                              nonzero_pi_rows=len(outputs),zero_pi_rows=zero,
                              best={k:best[k] for k in ['p','q','displayed_Q','weight_power']} if best else None))
    out=dict(status='passed exact bounded controls',elapsed_seconds=monotonic()-start,
             generation='Saturated two-row moment kernel; rational sample features and dyadic alpha/beta costs; guarded exact integer LLL; no pi or target in selection.',
             embedding='Clear a common denominator D of each feature block, append unscaled original H coefficients as guards; LLL rows reconstruct all features exactly.',
             summaries=summaries,candidates=candidates,
             scope='Fixed m=1,2,3 and five predeclared weights only; no asymptotic conclusion inferred.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='candidates'},indent=2))

if __name__=='__main__':run()
