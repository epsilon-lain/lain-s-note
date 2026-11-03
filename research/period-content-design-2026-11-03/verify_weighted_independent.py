"""Independent stdlib-only replay of the matched-weight finite certificates.

No generator, SymPy, Gram projection, Gauss reduction or original moment code
is imported.  Frozen source bases are treated as untrusted data and checked.
"""
from fractions import Fraction as F
from math import comb, factorial, gcd, isqrt, lcm
from pathlib import Path
import argparse
import hashlib
import json
import time


def rank(rows):
    a=[list(map(F,r)) for r in rows]
    piv=0
    for j in range(len(a[0])):
        k=next((k for k in range(piv,len(a)) if a[k][j]),None)
        if k is None:continue
        a[piv],a[k]=a[k],a[piv]
        lead=a[piv][j]
        a[piv]=[v/lead for v in a[piv]]
        for k in range(piv+1,len(a)):
            lead=a[k][j]
            if lead:a[k]=[v-lead*u for v,u in zip(a[k],a[piv])]
        piv+=1
        if piv==len(a):break
    return piv


def poly_from_primitive(vector,c):
    # vector=(r_1,...,r_(2c-1),alpha), with R(0)=0.
    n=[F(0)]*(2*c+1)
    for k,r in enumerate(vector[:-1],1):
        n[k-1]+=k*r
        n[k+1]+=(k-2*c+2)*r
    for j in range(c):n[2*j]+=comb(c-1,j)*vector[-1]
    return n


def poly_from_H(h,m,c):
    n=[F(0)]*(2*c+1)
    for k,v in enumerate(h):
        for j in range(2*m+1):
            n[2*m+k+j]+=(-1)**j*comb(2*m,j)*v
    return n


def check_identity(vector,h,m,c):
    n=poly_from_primitive(vector,c)
    assert n==poly_from_H(h,m,c)
    assert not any(n[:2*m])
    for k in range(2*m):
        assert sum(comb(j,k)*n[j] for j in range(k,len(n)))==0
    return n


def independent_moments(m,d):
    """Expand (1-x*x)^k, unlike the generator's positive (1+x)^k sum."""
    c=4*m+1;s=2*c;K=4*c
    rho=F(s+K+1,2*(K+2))
    assert 0<rho<1
    co=[F(comb(s+k-1,k),2**(s+k)) for k in range(K+1)]
    co.append(F(comb(s+K,K+1),2**(s+K+1))/(1-rho))
    def beta_integer(a,b):
        # 1 / ((a+b-1) choose (a-1)) / (a+b-1).
        return F(1,(a+b-1)*comb(a+b-2,a-1))
    out=[]
    for ell in range(2*d+1):
        v=F(0)
        for k,a in enumerate(co):
            moment=sum(((-1)**j*comb(k,j)*beta_integer(4*m+ell+2*j+1,4*m+1)
                        for j in range(k+1)),F(0))
            assert moment>0
            v+=a*moment
        assert v>0
        out.append(v)
    return out


def dot(a,b,mu):
    # First convolve, then integrate; no original inner-product helper is used.
    product=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):product[i+j]+=x*y
    return sum((v*u for v,u in zip(product,mu)),F(0))


def atan_interval(q,n):
    assert n%2==0
    term=F(1,q);total=F(0)
    for j in range(n):
        total+=term/(2*j+1)
        term/=-q*q
    return total,total+term/(2*n+1)


def pi_interval(n):
    a,b=atan_interval(5,n);u,v=atan_interval(239,n)
    return 16*a-4*v,16*b-4*u


def quality_interval(p,q,lo,hi):
    r=F(p,q)
    assert lo<hi
    if r<lo:return (lo-r)*q*q,(hi-r)*q*q
    assert r>hi
    return (r-hi)*q*q,(r-lo)*q*q


def sqrt_interval(x,scale=10**15):
    k=isqrt(x.numerator*scale*scale//x.denominator)
    return str(F(k,scale)),str(F(k+1,scale))


def run(folder,outpath):
    start=time.monotonic()
    path=folder/'weighted_results.json'
    data=json.loads(path.read_text())
    assert data['source_commit']=='4d84aa6c383e5eed6b0dc447946d689bd3bbe010'
    assert data['generation_reads_pi'] is False
    # Exact rational tangent verification of Machin's angle identity.
    tan=F(1,5)
    tan=2*tan/(1-tan*tan)
    tan=2*tan/(1-tan*tan)
    assert (tan-F(1,239))/(1+tan*F(1,239))==1
    plo,phi=pi_interval(180)
    tighterlo,tighterhi=pi_interval(182)
    assert plo<tighterlo<tighterhi<phi
    sources={}
    source_hashes={}
    for forced,name in [(False,'hilbert_quotient_lift.json'),(True,'hilbert_forced_period_lift.json')]:
        p=Path(__file__).resolve().parents[1]/'pi-period-quotient-design-2026-11-02'/'period-main-kernel'/name
        source_hashes[name]=hashlib.sha256(p.read_bytes()).hexdigest()
        sources[forced]={b['m']:b for b in json.loads(p.read_text())['per_m']}
    blocks=[];count=0
    for block in data['blocks']:
        m=block['m'];c=4*m+1;degree=block['degree'];forced=block['forced']
        assert block['majorant_truncation']==4*c
        src=sources[forced][m]
        assert src['c']==c
        mu=independent_moments(m,degree)
        assert list(map(F,block['moment_vector']))==mu
        zH=[list(map(F,v)) for v in src['zero_period_H']]
        zV=[list(map(F,v)) for v in src['zero_period_lift_coordinates']]
        C=[list(map(F,v)) for v in src['period_lift_H']]
        CV=[list(map(F,v)) for v in src['period_lift_coordinates']]
        H=[list(map(F,v)) for v in block['optimized_H']]
        V=[list(map(F,v)) for v in block['optimized_primitive']]
        P=src['period_basis']
        assert len(zH)==len(zV)==4*m
        assert rank(zV)==4*m
        # Zero-period constraints have precisely the advertised nullity.
        constraint_columns=[]
        for j in range(2*c):
            e=[0]*(2*c);e[j]=1;n=poly_from_primitive(e,c)
            constraint_columns.append(n[:2*m]+[
                sum(comb(k,h)*n[k] for k in range(h,len(n))) for h in range(2*m)
            ]+[F(j==2*c-1),F(j<2*c-1)])
        constraints=list(map(list,zip(*constraint_columns)))
        assert rank(constraints)==4*m+2
        for h,v in zip(zH,zV):
            check_identity(v,h,m,c)
            assert v[-1]==0 and sum(v[:-1])==0
        for j in range(2):
            check_identity(CV[j],C[j],m,c)
            check_identity(V[j],H[j],m,c)
            assert V[j][-1]==CV[j][-1]==P[0][j]
            assert sum(V[j][:-1])==sum(CV[j][:-1])==P[1][j]
            assert all(dot(H[j],z,mu)==0 for z in zH)
            # Together with the complete-rank constraints above, equal periods
            # prove the optimized representative lies in the same affine fiber.
        S=[[dot(a,b,mu) for b in H] for a in H]
        assert S==[list(map(F,row)) for row in block['schur_metric']]
        assert S[0][0]>0 and S[0][0]*S[1][1]-S[0][1]**2>0
        seen=set()
        for row in block['candidates']:
            assert row['m']==m and row['c']==c and row['forced']==forced
            p,q=row['p'],row['q'];alpha=row['alpha'];z=row['z']
            R=row['R_coefficients'];ih=row['H_coefficients'];iv=R+[alpha]
            assert all(isinstance(a,int) for a in iv+ih)
            check_identity(iv,ih,m,c)
            assert z==sum(R) and gcd(*iv)==1
            rat=F(-4*z,2**(c-1)*alpha)
            assert rat==F(p,q) and p>0 and q>0 and gcd(p,q)==1
            assert row['raw_fraction_gcd']==gcd(abs(4*z),abs(2**(c-1)*alpha))
            coord=row['quotient_coordinates'];L=row['clearing_denominator'];content=row['removed_coordinate_content']
            rv=[sum(coord[j]*V[j][k] for j in range(2)) for k in range(2*c)]
            rh=[sum(coord[j]*H[j][k] for j in range(2)) for k in range(degree+1)]
            assert L==lcm(*(a.denominator for a in rv))
            before=[int(a*L) for a in rv]
            assert all(a*L==b for a,b in zip(rv,before))
            assert content==gcd(*before) and content>0
            assert iv==[a//content for a in before]
            assert ih==[int(a*L/content) for a in rh]
            assert all(a*L/content==b for a,b in zip(rh,ih))
            assert rv[-1].denominator==1 and sum(rv[:-1]).denominator==1
            force=src['forced_period_factor']
            assert rv[-1].numerator%force==0 and sum(rv[:-1]).numerator%force==0
            norm=dot(rh,rh,mu)
            assert norm==F(row['rational_norm_before_clear'])
            certificate=16*F(q**4,alpha**2)*dot(ih,ih,mu)
            assert certificate==F(row['quality_squared_upper'])
            assert certificate==16*F(q**4,rv[-1]**2)*norm
            qlo,qhi=quality_interval(p,q,plo,phi)
            assert qlo==F(row['quality_lower']) and qhi==F(row['quality_upper'])
            tlo,thi=quality_interval(p,q,tighterlo,tighterhi)
            assert 0<qlo<tlo<thi<qhi and qhi*qhi<certificate
            assert (p,q) not in seen
            seen.add((p,q));count+=1
        candidates=block['candidates']
        assert all(F(a['quality_squared_upper'])<=F(b['quality_squared_upper']) for a,b in zip(candidates,candidates[1:]))
        chosen=block['certificate_selected'];best=block['postcheck_best']
        assert chosen==candidates[0]
        assert best==min(candidates,key=lambda r:F(r['quality_upper']))
        lo,hi=sqrt_interval(F(chosen['quality_squared_upper']))
        blocks.append({'m':m,'forced':forced,'verified_candidates':len(candidates),
            'certificate_selected':f"{chosen['p']}/{chosen['q']}",
            'selected_quality_lower':chosen['quality_lower'],'selected_quality_upper':chosen['quality_upper'],
            'certificate_bound_sqrt_lower':lo,'certificate_bound_sqrt_upper':hi,
            'postcheck_best':f"{best['p']}/{best['q']}",
            'postcheck_best_quality_upper':best['quality_upper'],
            'zero_period_rank':len(zH),'moments_recomputed':len(mu)})
        print('verified',m,forced,len(candidates),flush=True)
    assert count==119
    report={'status':'passed','total_candidates':count,'blocks':blocks,
        'checks':['primitive numerator polynomial identity','endpoint zeros to order2m',
            'exact source zero-period ranks and affine fibers','weighted zero-period orthogonality',
            'independent alternating-beta moment reconstruction','Schur metric positivity',
            'coefficient denominator clearing and integer content','actual fully reduced fraction and gcd',
            'strict rational Machin quality intervals','strict weighted quality-square bound',
            'predeclared certificate selection and postcheck ranking'],
        'independence':'Python standard library only; no generator, projection solver, original moment helper, or SymPy imports.',
        'data_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'frozen_source_sha256':source_hashes,
        'elapsed_seconds':time.monotonic()-start,
        'scope':'Finite replay only. Orthogonality proves rational minimum within each frozen period fiber, not global small-quality construction or an infinite theorem.'}
    outpath.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','candidates':count,'elapsed_seconds':report['elapsed_seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--folder',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('independent_weighted_verification.json'))
    args=parser.parse_args();run(args.folder,args.output)
