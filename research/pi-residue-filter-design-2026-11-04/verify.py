"""Independent standard-library replay; does not import the generators."""
from fractions import Fraction as F
from functools import lru_cache
from math import gcd,lcm,isqrt,comb,factorial
from itertools import permutations
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent

def ps(N):
    return [p for p in range(2,N+1) if all(p%d for d in range(2,isqrt(p)+1))]

@lru_cache(None)
def cf(s,t,n):
    if n==0:return F()
    c=F(s,t);d=F(2*n-1)
    for j in range(n-1,0,-1):d=2*j-1+j*j*c*c/d
    return c/d

@lru_cache(None)
def u_mod(s,t,n,p):
    # Independent coefficient formula for (1-2tz-s^2z^2)^(-1/2).
    D=s*s+t*t;inv4=pow(4,-1,p)
    return sum(comb(n,2*j)*comb(2*j,j)*pow(D,j,p)*pow(inv4,j,p)*pow(t,n-2*j,p)
               for j in range(n//2+1))%p

@lru_cache(None)
def W_exact(s,t,n):
    return sum(comb(n,2*j)*comb(2*j,j)*(s*s+t*t)**j*(2*t)**(n-2*j)
               for j in range(n//2+1))

def mod(r,p):
    assert r.denominator%p
    return r.numerator*pow(r.denominator,-1,p)%p

def det(M):
    size=len(M);v=F()
    for perm in permutations(range(size)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(size) for j in range(i+1,size))
        x=F(sign)
        for i,j in enumerate(perm):x*=M[i][j]
        v+=x
    return v

def cofactor_filter(params,H):
    r=H+2
    M=[[F(params[j][2]) for j in range(r)]]+[[cf(*params[j][:2],a) for j in range(r)] for a in range(1,H+1)]
    v=[(-1)**j*det([row[:j]+row[j+1:] for row in M]) for j in range(r)]
    assert any(v)
    D=lcm(*(x.denominator for x in v));C=[int(x*D) for x in v]
    g=gcd(*C);C=[c//g for c in C]
    K=sum(c*params[j][3] for j,c in enumerate(C))
    if K<0:C=[-c for c in C];K=-K
    assert K>0
    return C,K

def series_interval(k,N):
    total=F();power=F(1,k)
    for j in range(N):
        total+=power/(2*j+1)
        power/=-k*k
    other=total+power/(2*N+1)
    return min(total,other),max(total,other)

def quality(r,lo,hi):
    assert not lo<=r<=hi
    return min(abs(r-lo),abs(r-hi))*r.denominator**2,max(abs(r-lo),abs(r-hi))*r.denominator**2

def main():
    data=json.loads((ROOT/'results.json').read_text())
    data['rows']=sum((json.loads((ROOT/p).read_text()) for p in data['synchronous_files']),[])
    data['asynchronous_rows']=json.loads((ROOT/data['asynchronous_file']).read_text())
    params=data['parameters'];N=data['max_n'];primes=ps(2*N+3)
    a,b=F(2,5),F(3,7)
    assert (a+b)/(1-a*b)==1 and 0<a<1 and 0<b<1
    assert 2*a/(1-a*a)==F(20,21)
    assert 2*b/(1-b*b)==F(21,20)
    assert (b-a)/(1+a*b)==F(1,41)
    filters={}
    for row in data['filters']:
        C,K=cofactor_filter(params,row['H'])
        assert (C,K)==(row['coefficients'],row['K'])
        filters[row['H']]=(C,K)
    # Check the general Frobenius residue identity before using any filter.
    residue_cases=0
    for s,t,_,_ in params:
        for n in range(1,N+1):
            for p in primes:
                aa,bb=divmod(n,p)
                if not(p>2 and aa>=1 and p*p>n and s*t%p):continue
                if not u_mod(s,t,aa,p) or not u_mod(s,t,bb,p):continue
                assert mod(p*cf(s,t,n),p)==((-1)**((p-1)//2)*mod(cf(s,t,aa),p))%p
                residue_cases+=1
    A,B=series_interval(5,250),series_interval(239,70)
    lo,hi=16*A[0]-4*B[1],16*A[1]-4*B[0]
    assert [str(lo),str(hi)]==data['pi_interval'] and hi-lo<F(1,10**330)
    ref=quality(F(355,113),lo,hi)
    cancelled=0;removed=0
    for row in data['rows']:
        H,n=row['H'],row['n'];C,K=filters[H]
        values=[cf(s,t,n) for s,t,_,_ in params[:H+2]]
        D=lcm(*(v.denominator for v in values))
        P=4*sum(c*v.numerator*(D//v.denominator) for c,v in zip(C,values));Q=K*D
        r=F(P,Q);gg=gcd(abs(P),Q)
        assert str(P)==row['cost']['raw_numerator'] and str(Q)==row['cost']['raw_denominator']
        assert str(D)==row['cost']['common_denominator'] and str(gg)==row['cost']['actual_gcd']
        assert r==F(int(row['p']),int(row['q'])) and r.denominator==int(row['q'])
        # All new factors above n must come from K and the normalized
        # recurrence coefficients; they cannot arise from lcm(1,...,n).
        budget=K*lcm(*range(1,n+1))*lcm(*(W_exact(s,t,n) for s,t,_,_ in params[:H+2]))
        assert budget%r.denominator==0
        baseline=4*(cf(2,5,n)+cf(3,7,n))
        assert baseline==F(row['baseline'])
        qualified=[];eligible=[];exceptions=[];rm=[]
        for p in primes:
            if p==2:continue
            aa,bb=divmod(n,p)
            if not(1<=aa<=H and p*p>n):continue
            eligible.append(p)
            if K%p==0 or any(s*t%p==0 for s,t,_,_ in params[:H+2]):
                exceptions.append({'p':p,'reason':'parameter_or_period_denominator'});continue
            if any(not u_mod(s,t,aa,p) or not u_mod(s,t,bb,p) for s,t,_,_ in params[:H+2]):
                exceptions.append({'p':p,'reason':'normalized_coefficient_zero'});continue
            qualified.append(p);assert r.denominator%p
            if baseline.denominator%p==0:rm.append(p)
        assert qualified==row['certified_cancelled_primes'] and eligible==row['eligible_primes']
        assert exceptions==row['exceptions'] and rm==row['removed_relative_to_baseline']
        cancelled+=len(qualified);removed+=len(rm)
        big=r.denominator
        for p in primes:
            if p>n:break
            while big%p==0:big//=p
        assert str(big)==row['q_part_with_all_prime_factors_above_n']
        assert big.bit_length()==row['q_part_above_n_bits'] and r.denominator.bit_length()==row['q_bits']
        qlo,qhi=quality(r,lo,hi);blo,bhi=quality(baseline,lo,hi)
        assert (qhi<blo)==row['strictly_improves_baseline']
        assert (qlo>bhi)==row['strictly_worsens_baseline']
        assert row['strictly_worsens_baseline']
        assert (qlo>ref[1])==row['strictly_worse_than_355_over_113']
    cancelled_async=0
    for row in data['asynchronous_rows']:
        values=[cf(2,5,row['n']),cf(3,7,row['m']),cf(1,1,row['ell'])]
        assert row['ell']==(19*row['n']+9)//10
        D=lcm(*(v.denominator for v in values))
        P=4*(35*values[0]+35*values[1]-29*values[2])*D
        assert P.denominator==1
        r=F(P,6*D)
        assert r==F(int(row['p']),int(row['q']))
        blo,bhi=quality(4*(values[0]+values[1]),lo,hi)
        qlo,qhi=quality(r,lo,hi)
        assert (qhi<blo)==row['strictly_improves_baseline']
        assert (qlo>bhi)==row['strictly_worsens_baseline']
        assert row['strictly_worsens_baseline']
        assert qlo>ref[1]
        ns=[row['n'],row['m'],row['ell']];qualified=[]
        for p in primes:
            if p==2 or not(max(ns)<2*p and p<=min(ns) and p*p>max(ns)):continue
            if 6%p==0 or any(s*t%p==0 for s,t,_,_ in params[:3]):continue
            if any(not u_mod(s,t,1,p) or not u_mod(s,t,nn-p,p)
                   for (s,t,_,_),nn in zip(params[:3],ns)):continue
            assert r.denominator%p
            qualified.append(p)
        assert qualified==row['certified_cancelled_primes']
        cancelled_async+=len(qualified)
    control=json.loads((ROOT/'common_denominator_results.json').read_text())
    for row in control['rows']:
        n,N2=row['n'],row['N'];q=list(map(F,row['Q_coefficients']))
        Ps=[list(map(F,P)) for P in row['P_coefficients']]
        assert N2==3*n and len(q)==N2+1 and q[0]==1
        for c,P in zip([F(2,5),F(3,7),F(1)],Ps):
            fs=[F((-1)**j,2*j+1)*c**(2*j+1) for j in range(N2+n)]
            for j in range(N2+n):
                value=sum((q[i]*fs[j-i] for i in range(min(j,N2)+1)),F())
                assert value==(P[j] if j<N2 else 0)
        r=F(2,3)*sum((w*sum(P,F()) for w,P in zip([35,35,-29],Ps)),F())/sum(q,F())
        assert r==F(int(row['p']),int(row['q']))
        assert quality(r,lo,hi)[0]>1
    orbit=json.loads((ROOT/'orbit_control_results.json').read_text())
    v=[list(map(F,row)) for row in orbit['polynomial_weights']]
    angles=orbit['angles'];fs=[]
    for row in angles:
        A,B=row['A'],row['B'];D=A*A+B*B
        assert D==5**row['k'] and gcd(A,B)==1
        coefficients=[F(A*B,D)*F(4**j*factorial(j)**2,factorial(2*j+1))*F(B*B,D)**j for j in range(4)]
        assert coefficients==list(map(F,row['coefficients']))
        fs.append(coefficients)
    for n in range(1,4):assert sum((v[j][0]*fs[j][n]+v[j][1]*fs[j][n-1] for j in range(3)),F())==0
    assert sum((sum(v[j],F())*angles[j]['k'] for j in range(3)),F())==0
    assert sum((sum(v[j],F())*angles[j]['w'] for j in range(3)),F())==-1
    assert sum((sum(row,F()) for row in v),F())==0
    r=sum((v[j][0]*fs[j][0] for j in range(3)),F())
    assert r==F(orbit['candidate'])
    product=1
    for p,e in orbit['denominator_factors'].items():
        p=int(p)
        assert all(p%d for d in range(2,isqrt(p)+1))
        product*=p**e
    assert product==r.denominator and quality(r,lo,hi)[0]>10**26
    report={'status':'PASS','independent_residue_identity_cases':residue_cases,
         'independent_integer_filters':len(filters),'synchronous_candidates':len(data['rows']),
         'asynchronous_candidates':len(data['asynchronous_rows']),
         'qualified_prime_cancellations':cancelled,'removed_from_old_baseline':removed,
         'asynchronous_qualified_prime_cancellations':cancelled_async,
         'exact_remaining_denominator_budget_checks':len(data['rows']),
         'shared_denominator_controls':len(control['rows']),
         'gaussian_orbit_controls':1,
         'all_candidates_fail_to_improve_existing_quality':True,
         'scope':'Independent implementation, not an external review or a formal proof of the new theorem.'}
    (ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
