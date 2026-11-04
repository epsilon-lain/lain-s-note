"""Predesign prime-band cancellation using low-order arctan approximants.

Candidates and integer coefficients are generated before a pi interval is read.
SymPy is used only for small exact nullspaces.  The verifier uses the standard
library, a continued fraction, cofactors, and a separate series certificate.
"""
from fractions import Fraction as F
from math import gcd, lcm, isqrt
from pathlib import Path
from decimal import Decimal, localcontext
import json
import sympy as sp

ROOT=Path(__file__).resolve().parent
# atan(s/t) = alpha*theta + beta4*pi/4, theta=atan(2/5).
PARAMETERS=[(2,5,1,0),(3,7,-1,1),(1,1,0,1),
            (20,21,2,0),(21,20,-2,2),(1,41,-2,1)]
MAX_N=80

def primes(limit):
    return [p for p in range(2,limit+1) if all(p%d for d in range(2,isqrt(p)+1))]

def sequences(s,t,L):
    U,V=[1,-t],[0,-s]
    W=[1,2*t] # W_n=2^n u_n; u_n=(-1)^n U_n/n!.
    for n in range(1,L):
        U.append(-(2*n+1)*t*U[-1]+n*n*s*s*U[-2])
        V.append(-(2*n+1)*t*V[-1]+n*n*s*s*V[-2])
        a=2*(2*n+1)*t*W[-1]+4*n*s*s*W[-2]
        assert a%(n+1)==0
        W.append(a//(n+1))
    return [F(v,u) for v,u in zip(V,U)],W

def make_filter(H,tables):
    count=H+2
    M=sp.Matrix([[PARAMETERS[j][2] for j in range(count)]]+
       [[sp.Rational(tables[j][0][a]) for j in range(count)] for a in range(1,H+1)])
    null=M.nullspace()
    assert len(null)==1
    scale=lcm(*(x.q for x in null[0]))
    C=[int(x*scale) for x in null[0]]
    content=gcd(*C);C=[c//content for c in C]
    K=sum(C[j]*PARAMETERS[j][3] for j in range(count))
    assert K!=0
    if K<0:C=[-c for c in C];K=-K
    assert gcd(*C)==1
    assert sum(C[j]*PARAMETERS[j][2] for j in range(count))==0
    for a in range(1,H+1):assert sum((C[j]*tables[j][0][a] for j in range(count)),F())==0
    return C,K

def atan_bounds(k,N):
    a=sum((F((-1)**j,(2*j+1)*k**(2*j+1)) for j in range(N)),F())
    b=a+F((-1)**N,(2*N+1)*k**(2*N+1))
    return min(a,b),max(a,b)

def pi_bounds():
    a,b=atan_bounds(5,250),atan_bounds(239,70)
    # Exact tangent check plus 0<4 atan(1/5)-atan(1/239)<pi/2.
    v=F(1,5)
    for _ in range(2):v=2*v/(1-v*v)
    assert (v-F(1,239))/(1+v*F(1,239))==1
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]

def score(r,bounds):
    lo,hi=bounds;q=r.denominator
    low=F(0) if lo<=r<=hi else min(abs(lo-r),abs(hi-r))
    high=max(abs(lo-r),abs(hi-r))
    return q*q*low,q*q*high

def display(r):
    with localcontext() as c:
        c.prec=16
        return str(Decimal(r.numerator)/Decimal(r.denominator))

def reduce_cost(values,C,K):
    D=lcm(*(v.denominator for v in values))
    P=4*sum(c*v.numerator*(D//v.denominator) for c,v in zip(C,values))
    Q=K*D
    gg=gcd(abs(P),Q)
    return F(P,Q),{'common_denominator':str(D),'raw_numerator':str(P),
                   'raw_denominator':str(Q),'actual_gcd':str(gg),
                   'raw_denominator_bits':Q.bit_length()}

def main():
    tables=[sequences(s,t,2*MAX_N+3) for s,t,_,_ in PARAMETERS]
    ps=primes(2*MAX_N+3)
    filters=[];rows=[]
    for H in range(1,5):
        C,K=make_filter(H,tables)
        filters.append({'H':H,'angle_count':H+2,'coefficients':C,'K':K,
          'coefficient_max_bits':max(abs(c).bit_length() for c in C),'K_bits':K.bit_length()})
        for n in range(H+2,MAX_N+1):
            values=[tables[j][0][n] for j in range(H+2)]
            value,cost=reduce_cost(values,C,K)
            baseline=4*(tables[0][0][n]+tables[1][0][n])
            eligible=[];qualified=[];exceptions=[];removed=[]
            for p in ps:
                if p==2:continue
                a,b=divmod(n,p)
                if not(1<=a<=H and p*p>n):continue
                eligible.append(p)
                if K%p==0 or any((s*t)%p==0 for s,t,_,_ in PARAMETERS[:H+2]):
                    exceptions.append({'p':p,'reason':'parameter_or_period_denominator'});continue
                if any(tables[j][1][a]%p==0 or tables[j][1][b]%p==0 for j in range(H+2)):
                    exceptions.append({'p':p,'reason':'normalized_coefficient_zero'});continue
                assert value.denominator%p!=0
                qualified.append(p)
                if baseline.denominator%p==0:removed.append(p)
            large=value.denominator
            for p in ps:
                if p>n:break
                while large%p==0:large//=p
            rows.append({'H':H,'n':n,'p':str(value.numerator),'q':str(value.denominator),
                         'q_bits':value.denominator.bit_length(),
                         'source_q_bits':[v.denominator.bit_length() for v in values],
                         'baseline':str(baseline),'cost':cost,'eligible_primes':eligible,
                         'certified_cancelled_primes':qualified,'exceptions':exceptions,
                         'removed_relative_to_baseline':removed,
                         'q_part_with_all_prime_factors_above_n':str(large),
                         'q_part_above_n_bits':large.bit_length()})
    # Distinct-degree H=1 control. The fixed rational m-grid is only a finite
    # proposal rule, not a proof of asymptotic analytic alignment.
    async_rows=[]
    C,K=filters[0]['coefficients'],filters[0]['K']
    # This is a fixed rational proposal grid, not an assertion of asymptotic matching.
    for n in range(2,MAX_N+1):
        m0=25*n//24
        for m in range(max(1,m0-2),m0+4):
            if (n-m)%2==0:continue
            ell=(19*n+9)//10
            vals=[tables[0][0][n],tables[1][0][m],tables[2][0][ell]]
            value,cost=reduce_cost(vals,C,K)
            cancelled_async=[]
            for p in ps:
                if p==2 or not(max(n,m,ell)<2*p and p<=min(n,m,ell) and p*p>max(n,m,ell)):continue
                if K%p==0 or any((s*t)%p==0 for s,t,_,_ in PARAMETERS[:3]):continue
                orders=[n,m,ell]
                if any(tables[j][1][1]%p==0 or tables[j][1][orders[j]-p]%p==0 for j in range(3)):continue
                assert value.denominator%p!=0
                cancelled_async.append(p)
            async_rows.append({'n':n,'m':m,'ell':ell,'p':str(value.numerator),
                               'q':str(value.denominator),'cost':cost,
                               'certified_cancelled_primes':cancelled_async,
                               'baseline':str(4*(vals[0]+vals[1]))})
    bounds=pi_bounds() # Only now is pi used for independent post-generation scoring.
    ref=score(F(355,113),bounds)
    for row in rows+async_rows:
        r=F(int(row['p']),int(row['q']));lo,hi=score(r,bounds)
        blo,bhi=score(F(row['baseline']),bounds)
        assert lo>0 and blo>0
        row['quality_display']=display((lo+hi)/2)
        row['baseline_quality_display']=display((blo+bhi)/2)
        row['quality_ratio_display']=display(hi/blo)
        row['strictly_improves_baseline']=hi<blo
        row['strictly_worsens_baseline']=lo>bhi
        assert row['strictly_improves_baseline'] or row['strictly_worsens_baseline']
        row['strictly_worse_than_355_over_113']=lo>ref[1]
    summaries=[]
    for f in filters:
        local=[r for r in rows if r['H']==f['H']]
        best=min(local,key=lambda row:score(F(int(row['p']),int(row['q'])),bounds)[1])
        sample=next(row for row in local if row['n']==MAX_N)
        summaries.append({'H':f['H'],'rows':len(local),'best_n':best['n'],
          'best_fraction':str(F(int(best['p']),int(best['q']))),
          'best_quality_display':best['quality_display'],
          'baseline_improvements':sum(r['strictly_improves_baseline'] for r in local),
          'verified_prime_cancellations':sum(len(r['certified_cancelled_primes']) for r in local),
          'removed_from_baseline':sum(len(r['removed_relative_to_baseline']) for r in local),
          'sample_n80_q_bits':sample['q_bits'],'sample_n80_large_part_bits':sample['q_part_above_n_bits'],
          'sample_n80_cancelled_primes':sample['certified_cancelled_primes']})
    out={'schema':1,'source_head':'98a2b9fa0601da361cc76744c7c601b6d4b2af69',
         'parameters':[list(v) for v in PARAMETERS],'max_n':MAX_N,
         'candidate_generation_reads_pi':False,'filters':filters,'rows':rows,
         'asynchronous_rows':async_rows,'pi_interval':[str(x) for x in bounds],
         'summary':summaries,'asynchronous_summary':{'rows':len(async_rows),
             'improvements':sum(r['strictly_improves_baseline'] for r in async_rows),
             'qualified_prime_cancellations':sum(len(r['certified_cancelled_primes']) for r in async_rows)},
         'all_new_candidates_worse_than_reference':all(r['strictly_worse_than_355_over_113'] for r in rows+async_rows),
         'scope':'New coefficient selection and conditional prime-band cancellation; no successful critical pi construction.'}
    (ROOT/'data').mkdir(exist_ok=True)
    files=[]
    for H in range(1,5):
        path='data/synchronous-h'+str(H)+'.json'
        (ROOT/path).write_text(json.dumps([r for r in rows if r['H']==H],ensure_ascii=False,indent=2)+'\n')
        files.append(path)
    async_path='data/asynchronous.json'
    (ROOT/async_path).write_text(json.dumps(async_rows,ensure_ascii=False,indent=2)+'\n')
    del out['rows'];del out['asynchronous_rows']
    out['synchronous_files']=files;out['asynchronous_file']=async_path
    (ROOT/'results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'filters':filters,'summary':summaries,
                     'asynchronous_summary':out['asynchronous_summary'],
                     'all_new_candidates_worse_than_reference':out['all_new_candidates_worse_than_reference']},indent=2))

if __name__=='__main__':main()
