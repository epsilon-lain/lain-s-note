"""Pure square-root references and rigorous lower-envelope certificates."""
from fractions import Fraction as F
from math import isqrt,gcd
from decimal import Decimal,localcontext
from pathlib import Path
import json

ROOT=Path(__file__).parent
def dec(x):
    with localcontext() as c:
        c.prec=22
        return str(Decimal(x.numerator)/Decimal(x.denominator))

def run(max_m=20000):
    j=json.loads((ROOT/'pi_interval.json').read_text())
    L,U,S=map(int,(j['lower'],j['upper'],j['scale']))
    Plo,Phi=F(L,S),F(U,S)
    root_scale=10**80
    q_values=[1,2,5,10,20,50,100,200]
    best={q:(F(0),None) for q in q_values}
    record_q=0;records=[];changed_squares=0
    for m in range(2,max_m+1):
        n0=m*m*L*L//(S*S)
        assert n0==m*m*U*U//(S*S)
        r=isqrt(n0)
        n=n0-1 if r*r==n0 else n0
        changed_squares+=int(n!=n0)
        assert isqrt(n)**2!=n and 9*m*m<n<16*m*m
        root=isqrt(n*root_scale*root_scale)
        al,au=F(root,m*root_scale),F(root+1,m*root_scale)
        assert 3<al<=au<Plo
        eta_upper=Phi-al
        common=gcd(n,m*m)
        primitive_A=m*m//common
        primitive_B=n//common
        assert gcd(primitive_A,primitive_B)==1
        c_lower=1/(primitive_A*(4+au))
        ratio=c_lower/(2*eta_upper)
        q_half=isqrt(ratio.numerator//ratio.denominator)
        assert q_half*q_half*eta_upper<=c_lower/2
        if q_half>record_q or m in [2,10,100,1000,10000,max_m]:
            row={'m':m,'N':n,'raw_floor_was_square':n!=n0,
                 'polynomial_A':primitive_A,'polynomial_B':primitive_B,
                 'reference_lower_approx':dec(al),
                 'eta_upper_approx':dec(eta_upper),'c_lower_approx':dec(c_lower),
                 'eta_upper_exact':str(eta_upper),'c_lower_exact':str(c_lower),
                 'guaranteed_half_barrier_q_max':q_half}
            records.append(row)
            record_q=max(record_q,q_half)
        for q in q_values:
            h=max(F(0),c_lower/(q*q)-eta_upper)
            if h>best[q][0]:best[q]=(h,{'m':m,'N':n,'polynomial_A':primitive_A,'polynomial_B':primitive_B})
    return {'max_m':max_m,'certified_references':max_m-1,
            'perfect_square_floors_adjusted':changed_squares,
            'records':records,
            'selected_envelope':[
                {'q':q,'lower_bound_approx':dec(h),'lower_bound_exact':str(h),
                 'q_squared_lower_bound_approx':dec(q*q*h),'reference':ref}
                for q,(h,ref) in best.items()]}

if __name__=='__main__':
    result=run()
    (ROOT/'pure_root_results.json').write_text(json.dumps(result,indent=2))
    for row in result['records']:
        print('m',row['m'],'N',row['N'],'eta',row['eta_upper_approx'],
              'c',row['c_lower_approx'],'Q',row['guaranteed_half_barrier_q_max'])
    print('Adjusted square floors:',result['perfect_square_floors_adjusted'])
    for row in result['selected_envelope']:
        print('Envelope',row['q'],row['lower_bound_approx'],row['reference'])
