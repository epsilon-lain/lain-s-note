from fractions import Fraction as F
from math import isqrt, gcd
import json
from pathlib import Path


def atan_interval(x,n):
    terms=[(-1)**k*x**(2*k+1)/F(2*k+1) for k in range(n)]
    total=sum(terms,F(0)); nxt=(-1)**n*x**(2*n+1)/F(2*n+1)
    return min(total,total+nxt),max(total,total+nxt)


def pi_interval(n=25):
    a,b=atan_interval(F(1,5),n); c,d=atan_interval(F(1,239),n)
    return 16*a-4*d,16*b-4*c


def ceil(x): return -((-x.numerator)//x.denominator)
def floor(x): return x.numerator//x.denominator

def energy_bounds(p,q,s,a,b):
    u,v=p-q*b,p-q*a
    lo=F(0) if u<=0<=v else min(u*u,v*v)
    hi=max(u*u,v*v)
    return F(q*q)/s+s*lo,F(q*q)/s+s*hi


def spectral_candidates(a,b,s,L):
    # Exact complete enumeration whenever L/s<1; height bound comes from q^2/s<=L.
    assert L/s<1
    qmax=isqrt((s*L).numerator//(s*L).denominator)
    ans=[]; unknown=[]; count=0
    for q in range(1,qmax+1):
        for p in range(floor(a*q)-1,ceil(b*q)+2):
            count+=1
            lo,hi=energy_bounds(p,q,s,a,b)
            if hi<=L: ans.append((p,q,lo,hi))
            elif lo<=L: unknown.append((p,q,lo,hi))
    return qmax,count,ans,unknown


a,b=pi_interval(); s=F(4000000); L=F(7,1000)
qmax,count,ans,unknown=spectral_candidates(a,b,s,L)
assert [(p,q) for p,q,_,_ in ans]==[(355,113)]
assert not unknown
p,q,lo,hi=ans[0]
assert gcd(p,q)==1 and hi/2<F(35,10000)
# Below energy 1, independent integer directions cannot coexist. Check all enumerated
# integer pairs at several rational slopes/scales directly (no pi is needed).
tests=0
for alpha in [F(1,3),F(7,5),F(22,7),F(355,113),F(141421,100000)]:
    for scale in [F(1),F(3),F(11),F(100),F(10000)]:
        low=[]
        for q1 in range(1,min(200,isqrt(scale.numerator//scale.denominator))+1):
            center=floor(alpha*q1)
            for p1 in range(center-1,center+3):
                ev=F(q1*q1)/scale+scale*(p1-alpha*q1)**2
                if ev<1: low.append((p1,q1,ev))
        for p1,q1,e1 in low:
            for p2,q2,e2 in low:
                assert (p1*q2-p2*q1)==0
                tests+=1
# Exact determinant/variance identity for a finite weighted Fourier state.
pairs=[(3,1),(22,7),(355,113)]; weights=[F(1,10),F(2,10),F(7,10)]
A=sum(w*q*q for w,(p,q) in zip(weights,pairs))
B=sum(w*p*p for w,(p,q) in zip(weights,pairs))
C=sum(w*p*q for w,(p,q) in zip(weights,pairs))
Det=sum(weights[i]*weights[j]*(pairs[i][0]*pairs[j][1]-pairs[j][0]*pairs[i][1])**2
        for i in range(len(pairs)) for j in range(i+1,len(pairs)))
assert A*B-C*C==Det
alpha=F(314159,100000)
D=sum(w*(p-alpha*q)**2 for w,(p,q) in zip(weights,pairs))
assert A*D==A*A*(alpha-C/A)**2+Det
out={"status":"PASS","pi_interval_width":str(b-a),"scale":str(s),"threshold":str(L),
     "complete_height_bound":qmax,"integer_pairs_enumerated":count,
     "certified_directions":[{"p":p,"q":q,"lambda_lower":str(lo),"lambda_upper":str(hi),
                              "quality_upper_from_AMGM":str(hi/2)}],
     "unresolved":len(unknown),"rank_one_pair_checks":tests,"variance_identity":"PASS",
     "warning":"355/113 is a known benchmark; no new approximation or asymptotic theorem."}
Path(__file__).with_name('checks.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='certified_directions'},indent=2))
print('lambda approximately',float(lo),float(hi),'quality upper',float(hi/2))
