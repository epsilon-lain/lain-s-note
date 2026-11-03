"""Globally minimize each matched-weight certificate over primitive rationals.

Legendre's convergent criterion makes the finite rational CF list complete
once an incumbent square certificate is strictly below a/4.
"""
from fractions import Fraction as F
from pathlib import Path
from math import gcd,lcm,ceil
import json
from generate_weighted import encode, SOURCE_DIR


def convergents(x):
    pm2,pm1,qm2,qm1=0,1,1,0
    while True:
        t=x.numerator//x.denominator;p=t*pm1+pm2;q=t*qm1+qm2
        if q and p>0:yield p,q
        pm2,pm1,qm2,qm1=pm1,p,qm1,q
        rem=x-t
        if not rem:return
        x=1/rem


def period_form(block,original):
    S=[[F(x) for x in r] for r in block['schur_metric']]
    P=original['period_basis'];det=P[0][0]*P[1][1]-P[0][1]*P[1][0]
    Pinv=[[F(P[1][1],det),F(-P[0][1],det)],
          [F(-P[1][0],det),F(P[0][0],det)]]
    J=[[0,4],[-2**(4*block['m']),0]]
    T=[[sum(Pinv[i][k]*J[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    B=[[sum(T[k][i]*S[k][l]*T[l][j] for k in range(2) for l in range(2))
        for j in range(2)] for i in range(2)]
    return T,B


def run():
    path=Path(__file__).with_name('weighted_results.json')
    data=json.loads(path.read_text());out=[];summary=[];forms={}
    for block in data['blocks']:
        source='hilbert_forced_period_lift.json' if block['forced'] else 'hilbert_quotient_lift.json'
        original=json.loads((SOURCE_DIR/source).read_text())['per_m'][block['m']-1]
        T,B=period_form(block,original)
        a,b,d=B[0][0],B[0][1],B[1][1]
        theta=-b/a;epsilon=d-b*b/a
        assert a>0 and epsilon>0
        if block['m'] in forms:assert B==forms[block['m']]
        forms[block['m']]=B
        choices=[(q*q*(a*p*p+2*b*p*q+d*q*q),p,q) for p,q in convergents(theta)]
        score,p,q=min(choices)
        assert score<a/4 # completeness theorem applies
        v=[T[i][0]*p+T[i][1]*q for i in range(2)]
        V=[[F(x) for x in row] for row in block['optimized_primitive']]
        H=[[F(x) for x in row] for row in block['optimized_H']]
        vec=[v[0]*V[0][j]+v[1]*V[1][j] for j in range(len(V[0]))]
        h=[v[0]*H[0][j]+v[1]*H[1][j] for j in range(len(H[0]))]
        assert vec[-1]==4*q and sum(vec[:-1])==-2**(4*block['m'])*p
        L=lcm(*(x.denominator for x in vec));iv=[int(x*L) for x in vec];ih=[int(x*L) for x in h]
        assert all(x*L==j for x,j in zip(h,ih))
        content=gcd(*iv);iv=[x//content for x in iv];ih=[x//content for x in ih]
        assert F(-4*sum(iv[:-1]),2**(4*block['m'])*iv[-1])==F(p,q)
        upper=F(658,10000) if block['m']==1 else F(366 if block['m']==2 else 368,100000)
        assert score<upper*upper
        row={'m':block['m'],'forced':block['forced'],'a':a,'theta':theta,'epsilon':epsilon,
             'binary_form':B,'all_convergents_checked':len(choices),'p':p,'q':q,
             'quality_squared_upper':score,'strict_quality_upper':upper,
             'global_certificate_minimum_proved':True,
             'R_coefficients':iv[:-1],'H_coefficients':ih,'alpha':iv[-1],
             'clearing_denominator':L,'removed_coordinate_content':content}
        out.append(row)
        summary.append({k:row[k] for k in ['m','forced','all_convergents_checked','p','q',
            'strict_quality_upper','global_certificate_minimum_proved']})
    path.with_name('global_selected.json').write_text(json.dumps(encode(out),indent=2)+'\n')
    path.with_name('global_summary.json').write_text(json.dumps(encode(summary),indent=2)+'\n')
    print(json.dumps(encode(summary),indent=2))


if __name__=='__main__':run()
