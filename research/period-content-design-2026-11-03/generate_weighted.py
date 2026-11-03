"""Exact rational lifts in the matched-weight period quotient.

Input integer bases are frozen to repository commit 4d84aa6.  Generation
uses no pi value.  A Machin interval is read only in the final evaluation.
"""
from pathlib import Path
from weighted_rational import F, moments, inner as weighted_inner
import sys
SOURCE_DIR=Path(__file__).resolve().parents[1]/'pi-period-quotient-design-2026-11-02'/'period-main-kernel'
sys.path.insert(0,str(SOURCE_DIR))
from quotient_lift import gauss_two, base
from math import gcd,lcm
import sympy as sp
import json,time


def encode(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [encode(v) for v in value]
    return value


def optimized_lifts(block,mu):
    C=[list(map(F,row)) for row in block['period_lift_H']]
    Z=[list(map(F,row)) for row in block['zero_period_H']]
    CB=[list(map(F,row)) for row in block['period_lift_coordinates']]
    ZB=[list(map(F,row)) for row in block['zero_period_lift_coordinates']]
    dot=lambda a,b:weighted_inner(a,b,mu)
    A=sp.Matrix([[sp.Rational(dot(a,b)) for b in Z] for a in Z])
    B=sp.Matrix([[sp.Rational(dot(a,b)) for b in C] for a in Z])
    correction=-(A.inv()*B)
    H=[];V=[]
    for j in range(2):
        H.append([C[j][k]+sum((F(correction[i,j])*Z[i][k] for i in range(len(Z))),F(0))
                  for k in range(len(C[j]))])
        V.append([CB[j][k]+sum((F(correction[i,j])*ZB[i][k] for i in range(len(Z))),F(0))
                  for k in range(len(CB[j]))])
    for h in H:
        assert all(dot(h,z)==0 for z in Z)
    S=[[dot(a,b) for b in H] for a in H]
    assert S[0][0]>0 and S[0][0]*S[1][1]>S[0][1]**2
    return H,V,S


def run():
    start=time.monotonic();blocks=[]
    for forced in [False,True]:
        source='hilbert_forced_period_lift.json' if forced else 'hilbert_quotient_lift.json'
        old=json.loads((SOURCE_DIR/source).read_text())
        for block in old['per_m']:
            m=block['m'];c=block['c'];degree=len(block['period_lift_H'][0])-1
            mu=moments(m,degree)
            H,V,S=optimized_lifts(block,mu)
            P=block['period_basis'];seen={}
            for e in [0,4*m,8*m,12*m,16*m,20*m,24*m]:
                w=F(1,2**e)
                cost=[[S[i][j]+w*w*sum(P[k][i]*P[k][j] for k in range(2))
                       for j in range(2)] for i in range(2)]
                a,b=gauss_two(cost)
                for v in [a,b,[a[0]+b[0],a[1]+b[1]],[a[0]-b[0],a[1]-b[1]]]:
                    h=[v[0]*H[0][j]+v[1]*H[1][j] for j in range(len(H[0]))]
                    vec=[v[0]*V[0][j]+v[1]*V[1][j] for j in range(len(V[0]))]
                    alpha=vec[-1];z=sum(vec[:-1]);assert alpha.denominator==z.denominator==1
                    if not alpha:continue
                    rat=-4*z/(2**(c-1)*alpha)
                    if rat<=0:continue
                    L=lcm(*(t.denominator for t in vec))
                    iv=[int(t*L) for t in vec];ih=[int(t*L) for t in h]
                    assert all(t*L==k for t,k in zip(h,ih))
                    content=gcd(*iv)
                    iv=[t//content for t in iv];ih=[t//content for t in ih]
                    ialpha=iv[-1];iz=sum(iv[:-1]);assert rat==F(-4*iz,2**(c-1)*ialpha)
                    q=rat.denominator
                    norm=weighted_inner(h,h,mu)
                    certificate=16*F(q**4,alpha**2)*norm
                    assert certificate==16*F(q**4,ialpha**2)*weighted_inner(ih,ih,mu)
                    seen[str(rat)]={'m':m,'c':c,'forced':forced,'weight_exponent':e,
                        'quotient_coordinates':v,'p':rat.numerator,'q':q,
                        'alpha':ialpha,'z':iz,'R_coefficients':iv[:-1],
                        'H_coefficients':ih,'clearing_denominator':L,
                        'removed_coordinate_content':content,
                        'raw_fraction_gcd':gcd(abs(4*iz),abs(2**(c-1)*ialpha)),
                        'quality_squared_upper':certificate,'rational_norm_before_clear':norm}
            rows=sorted(seen.values(),key=lambda r:r['quality_squared_upper'])
            blocks.append({'m':m,'forced':forced,'degree':degree,'majorant_truncation':4*c,
                           'moment_vector':mu,'schur_metric':S,'optimized_H':H,
                           'optimized_primitive':V,'candidates':rows})
            print('generated',m,forced,len(rows),flush=True)
    # Strict rational pi certification occurs only after all candidate generation.
    a,b=base.atan(5);u,v=base.atan(239);plo,phi=16*a-4*v,16*b-4*u
    for block in blocks:
        for row in block['candidates']:
            r=F(row['p'],row['q']);q=row['q']
            if r<plo:lo,hi=(plo-r)*q*q,(phi-r)*q*q
            elif r>phi:lo,hi=(r-phi)*q*q,(r-plo)*q*q
            else:raise AssertionError('need more precision')
            assert hi*hi<=row['quality_squared_upper']
            row['quality_lower']=lo;row['quality_upper']=hi
            row['displayed_quality']=float(hi)
        block['certificate_selected']=block['candidates'][0]
        block['postcheck_best']=min(block['candidates'],key=lambda r:r['quality_upper'])
    out={'source_commit':'4d84aa6c383e5eed6b0dc447946d689bd3bbe010',
         'rule':'Matched rational weight; exact rational zero-period lift; seven fixed dyadic period penalties',
         'generation_reads_pi':False,'elapsed_seconds':time.monotonic()-start,'blocks':blocks,
         'scope':'Finite m=1,2,3 control; no infinite small-quality theorem, no new approximation record.'}
    Path(__file__).with_name('weighted_results.json').write_text(json.dumps(encode(out),indent=2)+'\n')
    print(json.dumps([{'m':b['m'],'forced':b['forced'],'candidates':len(b['candidates']),
         'certificate_selected':str(F(b['certificate_selected']['p'],b['certificate_selected']['q'])),
         'selected_Q':b['certificate_selected']['displayed_quality'],
         'selected_bound':float(b['certificate_selected']['quality_squared_upper'])**.5,
         'postcheck_best':str(F(b['postcheck_best']['p'],b['postcheck_best']['q'])),
         'best_Q':b['postcheck_best']['displayed_quality']}
         for b in blocks],indent=2),flush=True)


if __name__=='__main__':run()
