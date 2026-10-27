"""Certified quadratic reference numbers converging to pi; stdlib only."""
from fractions import Fraction as F
from math import isqrt, gcd
from decimal import Decimal, localcontext
from pathlib import Path
import json
from verify import atan_bounds

def dec(x):
    with localcontext() as c:
        c.prec=20
        return str(Decimal(x.numerator)/Decimal(x.denominator))

def run(steps=500):
    al,au=atan_bounds(5,450);bl,bu=atan_bounds(239,140)
    pl,pu=16*al-4*bu,16*au-4*bl
    S=10**600
    l=F(pl.numerator*S//pl.denominator,S)
    u=F(-(-pu.numerator*S//pu.denominator),S)
    assert pl<=l+F(1,S) and l<=pl<=pu<=u
    pi_l,pi_u=l,u
    scale=10**650;root=isqrt(5*scale*scale)
    phi_l,phi_u=(1+F(root,scale))/2,(1+F(root+1,scale))/2
    p2,p1,q2,q1=0,1,1,0
    M=1;coefficients=[];rows=[];records=[]
    selected={1,3,4,20,40,80,120,200,300,400,500}
    for k in range(steps+1):
        a=l.__floor__()
        if a!=u.__floor__():
            raise ArithmeticError('Uncertain CF coefficient; increase precision')
        coefficients.append(a)
        p,q=a*p1+p2,a*q1+q2
        assert abs(p*q1-p1*q)==1
        A=q1*q1+q1*q-q*q
        B=-2*p1*q1-p1*q-p*q1+2*p*q
        C=p1*p1+p1*p-p*p
        assert B*B-4*A*C==5
        assert gcd(gcd(abs(A),abs(B)),abs(C))==1
        new_record=k>0 and a>M
        if k>0:M=max(M,a)
        if k in selected or new_record:
            end1=(p*phi_l+p1)/(q*phi_l+q1)
            end2=(p*phi_u+p1)/(q*phi_u+q1)
            alpha_l,alpha_u=min(end1,end2),max(end1,end2)
            eta=max(abs(pi_l-alpha_u),abs(pi_u-alpha_l))
            c=F(1,M+2)
            limit=c/(2*eta)
            Q=isqrt(limit.numerator//limit.denominator)
            assert Q*Q*eta<=c/2
            assert (Q+1)*(Q+1)*eta>c/2
            # This quadratic lies inside the common-prefix CF cylinder.
            left,right=sorted((F(p,q),F(p+p1,q+q1)))
            assert left<alpha_l<=alpha_u<right
            assert left<pi_l<=pi_u<right
            row={'k':k,'a_k':a,'max_partial_quotient':M,
                'convergent_numerator':str(p),'convergent_denominator':str(q),
                'reference_matrix':[str(p),str(p1),str(q),str(q1)],
                'minimal_polynomial':[str(A),str(B),str(C)],
                'eta_upper_approx':dec(eta),
                'eta_upper_exact':[str(eta.numerator),str(eta.denominator)],
                'reference_c':str(c),'guaranteed_pi_c':str(c/2),
                'guaranteed_denominator_limit':str(Q)}
            rows.append(row)
            if new_record:records.append({'k':k,'M':M,'guaranteed_pi_c':str(c/2)})
        l,u=1/(u-a),1/(l-a)
        p2,p1,q2,q1=p1,p,q1,q
    return {'pi_enclosure_places':600,'certified_cf_coefficients':coefficients,
            'record_changes':records,'rows':rows,
            'claim':'finite denominator certificates only; no infinite uniform bound'}

if __name__=='__main__':
    result=run()
    Path(__file__).with_name('quadratic_results.json').write_text(json.dumps(result,indent=2))
    print('Record changes:',result['record_changes'])
    for row in result['rows']:
        Q=int(row['guaranteed_denominator_limit'])
        print('k=',row['k'],'M=',row['max_partial_quotient'],
              'eta<=',row['eta_upper_approx'],'pi_c=',row['guaranteed_pi_c'],
              'Q=',str(Q) if Q<10**12 else str(Q)[:12]+'... ('+str(len(str(Q)))+' digits)')
