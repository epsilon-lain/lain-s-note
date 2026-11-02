"""Exact finite audits for the generic recurrence-window witness limitation.

Uses only standard-library integer arithmetic. Infinite-index and formal
polynomial statements are proved in COMMON_WINDOW_WITNESS_BARRIER.md.
"""
from math import gcd
from pathlib import Path
import json


def odd(x):
    x=abs(x)
    assert x
    return x//(x&-x)


def main():
    maximum=300
    max_h=8
    U=[1,2]
    for n in range(1,maximum+max_h):
        numerator=2*(2*n+1)*U[n]+4*n*U[n-1]
        assert numerator%(n+1)==0
        U.append(numerator//(n+1))
    # Independently check the nonzero residue direction behind the stronger
    # linear-witness theorem at all candidate poles in a larger formal window.
    pole_directions=0
    nonzero_poles=0
    for t in range(1,81):
        N=-t
        a=[1,2*(2*N+1)]
        b=[0,4*N]
        for j in range(1,100):
            a.append(2*(2*N+2*j+1)*a[j]+4*(N+j)**2*a[j-1])
            b.append(2*(2*N+2*j+1)*b[j]+4*(N+j)**2*b[j-1])
        assert (-1)**t*b[t]>0
        nonzero_poles+=1
        for j in range(t,101):
            assert a[j]*b[t]==b[j]*a[t]
            pole_directions+=1
    checks=0
    linear=quadratic=0
    samples=[]
    for n in range(1,maximum+1):
        # D_j V_j = a_j x+b_j y. The first step is handled separately.
        aa=[1,2*(2*n+1)]
        bb=[0,4*n]
        dd=[1,n+1]
        for j in range(1,max_h):
            aa.append(2*(2*n+2*j+1)*aa[j]+4*(n+j)**2*aa[j-1])
            bb.append(2*(2*n+2*j+1)*bb[j]+4*(n+j)**2*bb[j-1])
            dd.append((n+j+1)*dd[j])
        C=odd(U[n])
        for h in range(1,max_h+1):
            assert dd[h]*U[n+h]==aa[h]*U[n]+bb[h]*U[n-1]
            checks+=1
            F=dd[h]*U[n+h]-bb[h]*U[n-1]
            assert F==aa[h]*U[n] and F
            core=C//gcd(C,dd[h])
            assert odd(F)%core==0
            assert odd(F)*dd[h]>=C
            linear+=1
            # A quadratic witness made from two different linear eliminants.
            # Its generic restriction at x=0 vanishes as a polynomial in y.
            previous=dd[h-1]*U[n+h-1]-bb[h-1]*U[n-1]
            Q=F*U[n+h]-previous*U[n+h-1]
            assert Q==U[n]*(aa[h]*U[n+h]-aa[h-1]*U[n+h-1])
            assert Q
            core2=C//gcd(C,dd[h]**2)
            assert odd(Q)%core2==0
            assert odd(Q)*dd[h]**2>=C
            quadratic+=1
            if (n,h) in ((4,1),(21,2),(100,4),(300,8)):
                samples.append({'n':n,'h':h,'D_h':str(dd[h]),
                    'central_odd_core':str(C),
                    'retained_divisor_for_degree_one':str(core),
                    'retained_divisor_for_degree_two':str(core2)})
    result={'status':'passed','arithmetic':'Python standard-library exact integers',
        'maximum_central_index':maximum,'maximum_forward_window':max_h,
        'exact_formal_window_coefficient_checks':checks,
        'nonzero_y_residue_poles':nonzero_poles,
        'collinear_residue_direction_checks':pole_directions,
        'linear_witness_odd_divisibility_checks':linear,
        'quadratic_witness_odd_divisibility_checks':quadratic,
        'samples':samples,
        'scope':'Finite algebra audits. The written theorem proves the universal generic-window limitation. No S(n) upper bound, no no-go theorem for actual-initial-data arithmetic, and no pi approximation conclusion.'}
    Path(__file__).with_name('common_window_witness_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='samples'},indent=2))


if __name__=='__main__':
    main()
