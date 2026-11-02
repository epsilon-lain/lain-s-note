"""Finite rational rotation, prearranged norm cancellation, exact two-round audit.

Standard library only. No decimal pi, continued fraction, or target fraction is
used to select the candidates. Pi bounds come from rational integral remainders.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext
from pathlib import Path
from math import gcd
import json

HERE=Path(__file__).resolve().parent

def pack(x): return [str(x.numerator),str(x.denominator)]
def dec(x):
    with localcontext() as ctx:
        ctx.prec=50
        return str(Decimal(x.numerator)/Decimal(x.denominator))

def atan_bounds(base,n):
    s=sum((F((-1)**j,(2*j+1)*base**(2*j+1)) for j in range(n)),F(0))
    t=F((-1)**n,(2*n+1)*base**(2*n+1))
    return min(s,s+t),max(s,s+t)

def pi_bounds(n):
    n239=max(1,(n+2)//3)
    a,b=atan_bounds(5,n)
    c,d=atan_bounds(239,n239)
    return 16*a-4*d,16*b-4*c,n239

def abs_bounds(a,b):
    return (F(0) if a<=0<=b else min(abs(a),abs(b))),max(abs(a),abs(b))

def quality(p,q,L,U):
    a,b=abs_bounds(q*q*L-p*q,q*q*U-p*q)
    return a,b

def info_fraction(p,q,L,U):
    lower,upper=quality(p,q,L,U)
    e0,e1=abs_bounds(L-F(p,q),U-F(p,q))
    return {'p':str(p),'q':str(q),'gcd':str(gcd(p,q)),
            'quality_lower':pack(lower),'quality_upper':pack(upper),
            'quality_display':dec((lower+upper)/2),
            'absolute_error_lower':pack(e0),'absolute_error_upper':pack(e1),
            'absolute_error_display':dec((e0+e1)/2)}

def select(p,q):
    assert q>1 and gcd(p,q)==1
    # Canonical lower Farey mate, computed only from p,q.
    d=pow(p,-1,q)
    c=(p*d-1)//q
    assert 0<d<q and p*d-q*c==1
    N=p*p+q*q
    b0=(N*d-p)//q
    assert 0<b0<N and (b0*q+p)==N*d and (b0*p-q)==N*c
    attempts=[]
    for n in range(3,100):
        L,U,n239=pi_bounds(n)
        e0,e1=p-q*U,p-q*L
        h0,h1=c-d*U,c-d*L
        if not(e0>0 and h1<0):
            attempts.append({'atan5_terms':n,'atan239_terms':n239,'reason':'opposite sides not yet certified'})
            continue
        z0=(d*L-c)/(p-q*L)
        z1=(d*U-c)/(p-q*U)
        assert 0<z0<z1
        k0,k1=z0.numerator//z0.denominator,z1.numerator//z1.denominator
        if k0!=k1:
            attempts.append({'atan5_terms':n,'atan239_terms':n239,'reason':'floor(z) unresolved'})
            continue
        ks=sorted(set([0,k0,k0+1]))
        candidates=[]
        for k in ks:
            P,Q=p*k+c,q*k+d
            lower,upper=quality(P,Q,L,U)
            candidates.append((k,P,Q,lower,upper))
        good=[row for row in candidates if all(row is other or row[4]<other[3] for other in candidates)]
        if len(good)!=1:
            attempts.append({'atan5_terms':n,'atan239_terms':n239,'reason':'candidate quality intervals overlap'})
            continue
        chosen=good[0]
        break
    else: raise AssertionError('adaptive selection failed')
    k,P,Q=chosen[:3]
    b=b0+N*k
    raw_P=b*p-q
    raw_Q=b*q+p
    assert raw_P==N*P and raw_Q==N*Q
    assert gcd(raw_P,raw_Q)==N and gcd(P,Q)==1
    assert p*Q-q*P==1
    # Verify the tangent-addition identity and its norm ledger exactly.
    r,t=F(p,q),F(-1,b)
    assert (r+t)/(1-r*t)==F(P,Q)
    assert raw_P*raw_P+raw_Q*raw_Q==(b*b+1)*N
    assert (b*b+1)%N==0
    # Algebraic error identity for the entire family, checked symbolically in pi.
    for j in (0,1,k0,k0+1,k0+2,2*k0+1):
        PP,QQ=p*j+c,q*j+d
        assert p*QQ-q*PP==1
        assert F(PP,QQ)==F(p,q)-F(1,q*QQ)
        for alpha in (L,U):
            assert (QQ*abs(PP-alpha*QQ))==(q*j+d)*abs((p-alpha*q)*j+(c-alpha*d))
    # High-precision scoring is after candidate selection, independently refined.
    Lscore,Uscore,_=pi_bounds(80)
    inp=info_fraction(p,q,Lscore,Uscore)
    out=info_fraction(P,Q,Lscore,Uscore)
    a0,a1=quality(p,q,Lscore,Uscore)
    bq0,bq1=quality(P,Q,Lscore,Uscore)
    ratio0,ratio1=bq0/a1,bq1/a0
    # Angle error is an integral over rationally bracketed endpoints.
    angle0=(F(p,q)-U)/(1+F(p,q)**2)
    angle1=(F(p,q)-L)/(1+L*L)
    assert 0<angle0<angle1
    tangent0=(p-q*U)/(q+p*U)
    tangent1=(p-q*L)/(q+p*L)
    ideal_b0=(q+p*L)/(p-q*L)
    ideal_b1=(q+p*U)/(p-q*U)
    assert ideal_b0==b0+N*z0 and ideal_b1==b0+N*z1
    for alpha in (L,U,Lscore,Uscore):
        ee=p-alpha*q
        zz=(alpha*d-c)/ee
        eta=k-zz
        assert (d+q*zz)*ee==1
        assert Q*abs(P-alpha*Q)==abs(eta)*(1+q*ee*eta)
    eta0=k-(Uscore*d-c)/(p-Uscore*q)
    eta1=k-(Lscore*d-c)/(p-Lscore*q)
    return {
        'input':inp,'canonical_mate':{'c':str(c),'d':str(d),'determinant':1},
        'norm_N':str(N),'b0':str(b0),
        'selection':{'atan5_terms':n,'atan239_terms':n239,
                     'pi_lower':pack(L),'pi_upper':pack(U),
                     'input_linear_error_interval':[pack(e0),pack(e1)],
                     'mate_linear_error_interval':[pack(h0),pack(h1)],
                     'z_interval':[pack(z0),pack(z1)],
                     'certified_floor_z':k0,'critical_k':ks,'selected_k':k,
                     'candidate_intervals':[dict(k=j,**info_fraction(PP,QQ,L,U)) for j,PP,QQ,_,_ in candidates],
                     'discarded_coarser_attempts':attempts},
        'angle':{'integral_error_lower':pack(angle0),'integral_error_upper':pack(angle1),
                 'tangent_error_lower':pack(tangent0),'tangent_error_upper':pack(tangent1),
                 'ideal_b_interval':[pack(ideal_b0),pack(ideal_b1)],
                 'chosen_tangent':pack(t),'b':str(b)},
        'rotation_raw':{'P':str(raw_P),'Q':str(raw_Q),'gcd':str(N),
                        'norm_identity':True,'full_norm_factor_prearranged':True},
        'output':out,'output_denominator_grows':Q>q,
        'quality_ratio_lower':pack(ratio0),'quality_ratio_upper':pack(ratio1),
        'quality_ratio_display':dec((ratio0+ratio1)/2),
        'fractional_alignment':{'eta_lower':pack(eta0),'eta_upper':pack(eta1),
                                'eta_display':dec((eta0+eta1)/2),
                                'exact_quality_transfer':'gamma_new=abs(eta)*(1+gamma*eta), gamma=q*(p-pi*q), eta=k-z'},
        'quality_strictly_improves':bq1<a0,
        'quality_strictly_worsens':bq0>a1,
        'global_chain_certificate':'For k<=z quality is concave quadratic; for k>=z it is strictly increasing. Thus k=0,floor z,ceil z contain the minimum over every integer k>=0.'}

def main():
    first=select(22,7)
    second=select(int(first['output']['p']),int(first['output']['q']))
    assert first['quality_strictly_improves']
    assert second['quality_strictly_worsens']
    result={'rounds':[first,second],
            'selection_inputs':'Initial 22/7; canonical Bezout mate; rational integral bounds. No decimal pi or target output fraction.',
            'limits':'Known fractions recovered; second same-rule step worsens quality; no infinite improvement or new pi record.'}
    path=HERE/'angle_chain_results.json'
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    reread=json.loads(path.read_text())
    for row in reread['rounds']:
        p,q=int(row['input']['p']),int(row['input']['q'])
        c,d=int(row['canonical_mate']['c']),int(row['canonical_mate']['d'])
        k=row['selection']['selected_k']
        assert (p*k+c,q*k+d)==(int(row['output']['p']),int(row['output']['q']))
        assert gcd(int(row['rotation_raw']['P']),int(row['rotation_raw']['Q']))==int(row['norm_N'])
    print(json.dumps([{'input':(r['input']['p'],r['input']['q']),
                       'output':(r['output']['p'],r['output']['q']),
                       'selected_k':r['selection']['selected_k'],
                       'atan_terms':(r['selection']['atan5_terms'],r['selection']['atan239_terms']),
                       'quality':r['output']['quality_display'],
                       'quality_ratio':r['quality_ratio_display'],
                       'gcd':r['rotation_raw']['gcd']} for r in reread['rounds']],indent=2))

if __name__=='__main__': main()
