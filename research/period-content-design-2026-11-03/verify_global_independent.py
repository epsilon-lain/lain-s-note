"""Independent replay of global minima of the rational weight certificate.

Imports only the independently written stdlib verifier, never the generator.
The binary form is integrated directly from target-period polynomial lifts.
"""
from fractions import Fraction as F
from math import gcd,lcm
from pathlib import Path
import hashlib
import json
import time
from verify_weighted_independent import (check_identity,independent_moments,dot,
    pi_interval,quality_interval,sqrt_interval)


def cf_candidates(theta):
    # Euclidean integer remainders with explicit unimodular matrix products.
    n,d=theta.numerator,theta.denominator
    mat=[[1,0],[0,1]]
    out=[];digits=[]
    while d:
        a,r=divmod(n,d);digits.append(a)
        mat=[[mat[0][0]*a+mat[0][1],mat[0][0]],
             [mat[1][0]*a+mat[1][1],mat[1][0]]]
        p,q=mat[0][0],mat[1][0]
        assert abs(mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0])==1
        if p>0 and q>0:out.append((p,q))
        n,d=d,r
    assert len(digits)==1 or digits[-1]>=2
    assert F(*out[-1])==theta
    return out


def run():
    start=time.monotonic()
    folder=Path(__file__).resolve().parent
    source_path=folder/'weighted_results.json'
    prior=json.loads(Path(__file__).with_name('independent_weighted_verification.json').read_text())
    assert prior['status']=='passed'
    assert hashlib.sha256(source_path.read_bytes()).hexdigest()==prior['data_sha256']
    data=json.loads(source_path.read_text())
    blocks={(b['m'],b['forced']):b for b in data['blocks']}
    selected_path=folder/'global_selected.json'
    selected=json.loads(selected_path.read_text())
    pi_lo,pi_hi=pi_interval(182)
    summaries=[];forms={};moment_cache={}
    for row in selected:
        m,forced=row['m'],row['forced'];c=4*m+1
        block=blocks[m,forced]
        name='hilbert_forced_period_lift.json' if forced else 'hilbert_quotient_lift.json'
        pth=Path(__file__).resolve().parents[1]/'pi-period-quotient-design-2026-11-02'/'period-main-kernel'/name
        assert hashlib.sha256(pth.read_bytes()).hexdigest()==prior['frozen_source_sha256'][name]
        original=next(b for b in json.loads(pth.read_text())['per_m'] if b['m']==m)
        P=original['period_basis'];det=P[0][0]*P[1][1]-P[0][1]*P[1][0]
        H=[list(map(F,r)) for r in block['optimized_H']]
        V=[list(map(F,r)) for r in block['optimized_primitive']]
        def target_lift(alpha,z):
            # Solve P v=(alpha,z) by Cramer's rule, without using stored S or T.
            v=[F(alpha*P[1][1]-P[0][1]*z,det),
               F(P[0][0]*z-alpha*P[1][0],det)]
            h=[v[0]*H[0][k]+v[1]*H[1][k] for k in range(len(H[0]))]
            vec=[v[0]*V[0][k]+v[1]*V[1][k] for k in range(len(V[0]))]
            assert vec[-1]==alpha and sum(vec[:-1])==z
            check_identity(vec,h,m,c)
            return h,vec
        hp,vp=target_lift(0,-2**(4*m))
        hq,vq=target_lift(4,0)
        if m not in moment_cache:moment_cache[m]=independent_moments(m,block['degree'])
        mu=moment_cache[m]
        binary=[[dot(hp,hp,mu),dot(hp,hq,mu)],
                [dot(hq,hp,mu),dot(hq,hq,mu)]]
        assert binary==[list(map(F,r)) for r in row['binary_form']]
        a,b,d=binary[0][0],binary[0][1],binary[1][1]
        theta=-b/a;epsilon=d-b*b/a
        assert a==F(row['a']) and theta==F(row['theta']) and epsilon==F(row['epsilon'])
        assert a>0 and epsilon>0 and theta>F(1,2)
        if m in forms:assert forms[m]==binary
        forms[m]=binary
        candidates=cf_candidates(theta)
        assert len(candidates)==row['all_convergents_checked']
        def objective(p,q):
            result=q*q*(a*p*p+2*b*p*q+d*q*q)
            assert result==a*q*q*(p-theta*q)**2+epsilon*q**4
            return result
        choices=[(objective(p,q),p,q) for p,q in candidates]
        score,p,q=min(choices)
        assert score==F(row['quality_squared_upper'])
        assert p==row['p'] and q==row['q'] and gcd(p,q)==1
        assert score<a/4
        # Completeness: every rational with objective<=score then satisfies
        # a*q^2*(p-theta*q)^2 < a/4, so |p/q-theta|<1/(2q^2).
        # Legendre's theorem forces a canonical convergent of rational theta;
        # theta itself is included as the final convergent.
        assert row['global_certificate_minimum_proved'] is True
        h=[p*x+q*y for x,y in zip(hp,hq)]
        vec=[p*x+q*y for x,y in zip(vp,vq)]
        assert vec[-1]==4*q and sum(vec[:-1])==-2**(4*m)*p
        assert q*q*dot(h,h,mu)==score
        L=lcm(*(v.denominator for v in vec))
        before=[int(v*L) for v in vec]
        assert all(v*L==k for v,k in zip(vec,before))
        content=gcd(*before)
        iv=[v//content for v in before]
        ih=[int(v*L/content) for v in h]
        assert all(v*L/content==k for v,k in zip(h,ih))
        assert L==row['clearing_denominator'] and content==row['removed_coordinate_content']
        assert iv==row['R_coefficients']+[row['alpha']]
        assert ih==row['H_coefficients'] and gcd(*iv)==1
        check_identity(iv,ih,m,c)
        rat=F(-4*sum(iv[:-1]),2**(4*m)*iv[-1])
        assert rat==F(p,q)
        g=gcd(abs(4*sum(iv[:-1])),abs(2**(4*m)*iv[-1]))
        assert abs(2**(4*m)*iv[-1])//g==q
        cert=16*F(q**4,iv[-1]**2)*dot(ih,ih,mu)
        assert cert==score
        qlo,qhi=quality_interval(p,q,pi_lo,pi_hi)
        assert qhi*qhi<score<F(row['strict_quality_upper'])**2
        sl,sh=sqrt_interval(score)
        summaries.append({'m':m,'forced':forced,'fraction':f'{p}/{q}',
            'convergents_checked':len(candidates),'certificate_global_minimum':True,
            'strict_quality_upper':row['strict_quality_upper'],
            'certificate_sqrt_lower':sl,'certificate_sqrt_upper':sh,
            'actual_quality_lower':str(qlo),'actual_quality_upper':str(qhi),
            'actual_raw_gcd':str(g),
            'primitive_coordinates':len(iv),'identity_and_endpoints_passed':True})
        print('global verified',m,forced,f'{p}/{q}',len(candidates),flush=True)
    assert len(summaries)==6
    result={'status':'passed','selections':summaries,
        'binary_form_reconstruction':'Direct independent moment integration of target-period polynomial lifts, without using stored Schur matrix or projection implementation.',
        'global_completeness':'All finite canonical convergents of rational theta checked. Incumbent<a/4 forces every equally good rational to satisfy the strict Legendre criterion, hence to be on that finite list.',
        'scope':'Global minimum of this particular rational majorant certificate within rationally optimized frozen period fibers; does not minimize true pi approximation quality or prove an infinite small-quality theorem.',
        'forced_scope':'Rational affine minimization is homogeneous and yields identical binary forms for forced/unforced bases. Prescribed integral quotient coefficients are not claimed.',
        'global_data_sha256':hashlib.sha256(selected_path.read_bytes()).hexdigest(),
        'verified_weighted_data_sha256':prior['data_sha256'],
        'elapsed_seconds':time.monotonic()-start}
    dest=Path(__file__).with_name('independent_global_verification.json')
    dest.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','selections':len(summaries),'elapsed_seconds':result['elapsed_seconds']},indent=2))


if __name__=='__main__':run()
