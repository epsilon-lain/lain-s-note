"""Independent audit of the prime-power correction theorem; standard library."""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import json

def primes(n):
    return [p for p in range(3,n+1,2) if all(p%d for d in range(2,isqrt(p)+1))]

def atan_interval(d,n):
    value=F()
    for j in range(n):
        value+=F((-1)**j,(2*j+1)*d**(2*j+1))
    other=value+F((-1)**n,(2*n+1)*d**(2*n+1))
    return min(value,other),max(value,other)

def main():
    limit=2001
    powers={}
    for p in primes(limit):
        m=p;r=1
        while m<=limit:
            powers[m]=(p,r)
            m*=p;r+=1
    x,y=atan_interval(5,90),atan_interval(239,30)
    lo,hi=16*x[0]-4*y[1],16*x[1]-4*y[0]
    assert lo<hi
    S=F();count=0;exceptions=[];samples=[]
    for m in range(1,limit+1):
        S+=F((-1)**(m-1),2*m-1)
        if m not in powers:continue
        p,r=powers[m]
        U=4*S;V=U-F(1,m)
        qu,qv=U.denominator,V.denominator
        assert qu%m==0 and qu%qv==0
        special=(p==3 and r%2==0)
        assert qu//qv==(3 if special else 1)
        if special:exceptions.append({'m':m,'p':p,'r':r,'denominator_loss':3})
        E0,E1=U-hi,U-lo
        d0,d1=lo-V,hi-V
        assert E0>0 and d0>0
        assert d0>=F(1,4*m**3)-F(5,16*m**5)
        assert d1<=F(1,4*m**3)
        actual_ratio_upper=F(qv*qv,qu*qu)*d1/E0
        assert actual_ratio_upper<=F(1,4*m*m-1)
        if m in [3,5,7,9,27,81,729,2001]:
            samples.append({'m':m,'p':p,'r':r,
                'q_U_bits':qu.bit_length(),'q_V_bits':qv.bit_length(),
                'denominator_ratio_U_over_V':qu//qv,
                'relative_quality_upper_bound':str(F(1,4*m*m-1))})
        count+=1
    report={'status':'PASS','limit':limit,'odd_prime_powers_audited':count,
        'exact_denominator_classification':True,
        'rational_error_bounds_and_actual_quality_ratios':count,
        'exceptions':exceptions,'selected_samples':samples,
        'scope':'Infinite relative gain theorem audited; no claim that corrected absolute quality tends to zero.'}
    Path(__file__).with_name('audit_common_denominator.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
