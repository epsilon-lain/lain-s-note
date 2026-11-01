"""Exact finite-jet shear certificates, Python standard library only.

General statements are proved in TAIL_CANCELLATION.md; finite checks below
cover nonzero minors, formal tail parameters, threshold equality/failure,
and Gaussian-integral arithmetic normalization. No floating-point ranks.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, lcm
from pathlib import Path
import json


def ceilq(x):
    x=F(x); return -((-x.numerator)//x.denominator)

def wt(a,w): return sum(x*y for x,y in zip(a,w))

def rows(H,w,theta,v0,K):
    out=[]
    for j in range(K):
        for beta in product(*(range(ceilq(H*theta/wi)) for wi in w)):
            for s in range(max(0,ceilq((H-wt(beta,w)/theta)/v0))):
                out.append((j,s,beta))
    return sorted(out,key=lambda x:(x[0],x[1],x[2]))

def columns(H,w,w0):
    return [(h,a) for a in product(*(range(int(H/wi)+1) for wi in w))
            for h in range(max(0,int((H-wt(a,w))/w0)+1))]

def add(A,B):
    out=dict(A)
    for k,x in B.items():
        out[k]=out.get(k,F())+x
        if not out[k]: del out[k]
    return out

def mul(A,B,N):
    out={}
    for (t,z),a in A.items():
        for (u,y),b in B.items():
            if t+u<=N: out[t+u,z+y]=out.get((t+u,z+y),F())+a*b
    return {k:x for k,x in out.items() if x}

def powpoly(A,e,N):
    out={(0,0):F(1)}
    for _ in range(e):out=mul(out,A,N)
    return out

def matrix(cc,rr,centers,tails):
    N=max(s for j,s,b in rr)
    maxpow=max(sum(a) for h,a in cc)
    cache={}
    for j in set(r[0] for r in rr):
        for i,c in enumerate(centers):
            p={(0,0):j*c}
            p.update({(k,0):F((-1)**(k+1),k) for k in range(1,N+1)})
            p.update({(k,1):x for k,x in tails[i].items() if k<=N and x})
            for d in range(maxpow+1):cache[j,i,d]=powpoly(p,d,N)
    out=[]
    for j,s,beta in rr:
        row=[]
        for h,alpha in cc:
            if any(a<b for a,b in zip(alpha,beta)):
                row.append({});continue
            p={(k,0):F(comb(h,k)) for k in range(min(h,N)+1)}
            factor=1
            for i,(a,b) in enumerate(zip(alpha,beta)):
                factor*=comb(a,b);p=mul(p,cache[j,i,a-b],N)
            row.append({z:factor*x for (t,z),x in p.items() if t==s and x})
        out.append(row)
    return out

def shear(rr,tails):
    N=max(s for j,s,b in rr)
    out=[]
    for j,s,beta in rr:
        line=[]
        for jj,ell,a in rr:
            if jj!=j or ell>s or any(x<y for x,y in zip(a,beta)):
                line.append({});continue
            p={(0,0):F(1)};factor=1
            for ai,bi,tau in zip(a,beta,tails):
                factor*=comb(ai,bi)
                p=mul(p,powpoly({(k,1):x for k,x in tau.items()},ai-bi,N),N)
            line.append({z:factor*x for (t,z),x in p.items() if t==s-ell and x})
        out.append(line)
    return out

def evaluate(A,z):return [[sum(c*z**k for k,c in p.items()) for p in row] for row in A]

def rankcols(A):
    a=[r[:] for r in A];r=0;piv=[]
    for c in range(len(a[0])):
        k=next((i for i in range(r,len(a)) if a[i][c]),None)
        if k is None:continue
        a[r],a[k]=a[k],a[r]
        v=a[r][c];a[r]=[x/v for x in a[r]]
        for i in range(r+1,len(a)):
            v=a[i][c]
            if v:a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        piv.append(c);r+=1
        if r==len(a):break
    return r,piv

def det(A):
    a=[r[:] for r in A];d=F(1)
    for c in range(len(a)):
        k=next((i for i in range(c,len(a)) if a[i][c]),None)
        if k is None:return F()
        if k!=c:a[c],a[k]=a[k],a[c];d=-d
        v=a[c][c];d*=v
        for i in range(c+1,len(a)):
            x=a[i][c]/v
            if x:
                for j in range(c+1,len(a)):a[i][j]-=x*a[c][j]
    return d

def period_translation(cc,rr,centers,tails,W=F(2)):
    """Expand constants and formal z-scaled tails around a common period symbol value.

    W=2 is an exact rational evaluation of the algebraic period-reduced formula;
    no claim exp(2)=1 is used. The underlying identity is polynomial in W.
    """
    N=max(s for j,s,b in rr); out=[]
    for j,s,beta in rr:
        line=[]
        L={(0,0):j*W}
        L.update({(k,0):F((-1)**(k+1),k) for k in range(1,N+1)})
        for h,alpha in cc:
            value={}
            if all(x>=y for x,y in zip(alpha,beta)):
                for a in product(*(range(b,al+1) for b,al in zip(beta,alpha))):
                    factor=1
                    P={(k,0):F(comb(h,k)) for k in range(min(h,N)+1)}
                    P=mul(P,powpoly(L,sum(alpha)-sum(a),N),N)
                    for i,(ai,bi,al) in enumerate(zip(a,beta,alpha)):
                        factor*=comb(ai,bi)*comb(al,ai)
                        shift={(0,0):j*(centers[i]-W)}
                        shift.update({(k,1):v for k,v in tails[i].items()})
                        P=mul(P,powpoly(shift,ai-bi,N),N)
                    for (t,z),x in P.items():
                        if t==s:value[z]=value.get(z,F())+factor*x
            line.append({z:x for z,x in value.items() if x})
        out.append(line)
    return out


def shear_case(H,w,theta,v0,K,centers):
    w0=F(1);rr=rows(H,w,theta,v0,K);cc=columns(H,w,w0)
    T=[ceilq(wi/(theta*v0)) for wi in w]
    # Nontrivial generic rational tails beginning exactly at the threshold.
    tails=[{ti:F(i+2,3),ti+1:F(-i-1,5)} for i,ti in enumerate(T)]
    A=matrix(cc,rr,centers,tails);U=shear(rr,tails);A0=evaluate(A,0)
    assert period_translation(cc,rr,centers,tails)==A
    M=len(rr);nontrivial=0;coefficient_checks=0
    for i in range(M):
        assert U[i][i]=={0:F(1)}
        assert all(not U[i][k] for k in range(i+1,M))
        nontrivial+=sum(bool(U[i][k]) for k in range(i))
        for col in range(len(cc)):
            v={}
            for k in range(M):
                for z,x in U[i][k].items():v[z]=v.get(z,F())+x*A0[k][col]
            v={k:x for k,x in v.items() if x}
            assert v==A[i][col];coefficient_checks+=len(set(v)|set(A[i][col]))
    rank,ix=rankcols(A0);assert rank==M
    D=det([[r[i] for i in ix] for r in A0]);assert D
    degree_bound=sum(max((max(p,default=0) for p in row),default=0) for row in A)
    # A determinant has degree at most this sum. Degree+1 exact values certify
    # the entire determinant polynomial, not just a few numerical samples.
    for z in range(degree_bound+1):
        Az=evaluate(A,z)
        assert det([[r[i] for i in ix] for r in Az])==D
    # Actual Gi cutoff has tau_i=-sum_{k>=Ti} (-1)^(k+1)t^k/k.
    N=max(s for j,s,b in rr)
    cutoff=[{k:-F((-1)**(k+1),k) for k in range(ti,N+1)} for ti in T]
    G=matrix(cc,rr,centers,cutoff)
    assert det([[r[i] for i in ix] for r in evaluate(G,1)])==D
    return {'H':H,'weights':list(w),'theta':str(theta),'v0':str(v0),'K':K,
            'rows':M,'columns':len(cc),'thresholds':T,'equality_threshold':any(v0*t==wi/theta for wi,t in zip(w,T)),
            'off_diagonal_shear_entries':nontrivial,'formal_coefficient_checks':coefficient_checks,
            'minor_columns':ix,'determinant':str(D),'cutoff_determinant_equal':True,
            'full_constant_and_tail_translation_checked':True,
            'determinant_polynomial_degree_bound':degree_bound,
            'exact_interpolation_values_checked':degree_bound+1}

# Gaussian-rational pairs, exact.
def ga(a,b):return (a[0]+b[0],a[1]+b[1])
def gm(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def gscale(a,c):return (a[0]*c,a[1]*c)
def gmul(A,B,N):
    out=[(F(),F()) for _ in range(N+1)]
    for i,a in enumerate(A):
        for j,b in enumerate(B[:N+1-i]):out[i+j]=ga(out[i+j],gm(a,b))
    return out

def gaussian_integrality_case():
    H=7;w=(2,3);theta=F(1,2);v0=F(1);K=2
    rr=rows(H,w,theta,v0,K);cc=columns(H,w,F(1));N=max(s for j,s,b in rr)
    q=(7,13);p=(22,41);T=[ceilq(wi/(theta*v0)) for wi in w]
    L=[lcm(*range(1,t)) if t>1 else 1 for t in T]
    scale=[qi*li for qi,li in zip(q,L)];checked=nonreal=0
    for j,s,beta in rr:
        for h,alpha in cc:
            if any(a<b for a,b in zip(alpha,beta)):continue
            poly=[(F(comb(h,k)),F()) if k<=h else (F(),F()) for k in range(N+1)]
            factor=F(1)
            for i,(a,b) in enumerate(zip(alpha,beta)):
                f=[(F(),F()) for _ in range(N+1)];f[0]=(F(),F(2*j*p[i],q[i]))
                for k in range(1,min(T[i],N+1)):f[k]=(F((-1)**(k+1),k),F())
                for _ in range(a-b):poly=gmul(poly,f,N)
                factor*=comb(a,b)*scale[i]**(a-b)
            value=gscale(poly[s],factor)
            assert all(x.denominator==1 for x in value)
            checked+=1;nonreal+=bool(value[1])
    return {'rows':len(rr),'columns':len(cc),'q':list(q),'p':list(p),'L':L,'thresholds':T,
            'nonzero_support_entries_checked':checked,'entries_with_nonzero_imaginary_part':nonreal,
            'all_scaled_entries_gaussian_integral':True}


def run():
    # These exact examples are about the shear identity, not unknown interpolation thresholds.
    specs=[(3,(1,),F(1,2),F(1),1,(F(3),)),
           (4,(1,),F(1,2),F(1),2,(F(3),)),
           (7,(2,3),F(1,2),F(1),2,(F(3),F(5))),
           (5,(1,2),F(1,3),F(2),2,(F(3),F(5)))]
    cases=[shear_case(*s) for s in specs]
    # Minimal failure when the start order is below the weight threshold.
    rr=rows(2,(1,),F(1,2),F(1),1);cc=[(0,(0,)),(0,(1,))]
    A=matrix(cc,rr,(F(3),),[{1:F(1)}])
    assert [det(evaluate(A,z)) for z in (0,1,2)]==[1,2,3]
    # Completeness matters even at the valid threshold T=2.
    sparse=[(0,2,(0,))]
    B=matrix([(0,(1,))],sparse,(F(3),),[{2:F(1,2)}])
    assert evaluate(B,0)==[[F(-1,2)]] and evaluate(B,1)==[[0]]
    result={'status':'passed','cases':cases,'below_threshold_determinant':'1+z',
            'incomplete_packet_log_entry':'-1/2','incomplete_packet_truncated_entry':'0',
            'arithmetic':gaussian_integrality_case(),
            'scope':'Exact finite checks of algebraic identities; not a proof of pi bad approximability or of manuscript Theorem 2.1.'}
    (Path(__file__).parent/'tail_cancellation_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':run()
