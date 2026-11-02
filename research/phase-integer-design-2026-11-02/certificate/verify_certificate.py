"""Independent exact positive-integral error certificates; standard library."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
from decimal import Decimal, localcontext
import json

def beta3(s):
    return F(2,s*(s+1)*(s+2))

def bounds(m,K=20):
    assert m>=1 and m%2==1 and K>=0
    center=F(4,3)*sum((sum((comb(k,j)*(-2)**j*beta3(2*m+2*j)
                    for j in range(k+1)),F())/3**k for k in range(K+1)),F())
    radius=2*beta3(2*m)/3**(K+1)
    return center-radius,center+radius

def independent_center(m,K):
    coefficients=[1]
    result=F()
    for k in range(K+1):
        # Integrate x^(2m-1) (1-x)^2 (1-2x^2)^k directly.
        moment=sum((F(a,2*m+j)-2*F(a,2*m+j+1)+F(a,2*m+j+2)
                    for j,a in enumerate(coefficients)),F())
        result+=F(4,3**(k+1))*moment
        nxt=[0]*(len(coefficients)+2)
        for j,a in enumerate(coefficients):
            nxt[j]+=a
            nxt[j+2]-=2*a
        coefficients=nxt
    return result

def atan_bounds(k,n):
    total=sum((F((-1)**j,(2*j+1)*k**(2*j+1)) for j in range(n)),F())
    nxt=total+F((-1)**n,(2*n+1)*k**(2*n+1))
    return min(total,nxt),max(total,nxt)

def machin():
    a,b=atan_bounds(5,70),atan_bounds(239,20)
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]

def display(x):
    with localcontext() as ctx:
        ctx.prec=14
        return str(Decimal(x.numerator)/Decimal(x.denominator))

def source(m):
    return 4*sum((F((-1)**j,2*j+1) for j in range(m)),F())

def main():
    K=20
    pi_lo,pi_hi=machin()
    tests=[]
    for m in range(1,100,2):
        lo,hi=bounds(m,K)
        assert (lo+hi)/2==independent_center(m,K)
        U=source(m);V=U-F(1,m)
        assert lo>0 and hi<F(1,m)
        assert lo<pi_lo-V<pi_hi-V<hi
        assert hi-lo==4*beta3(2*m)/3**(K+1)
        tests.append(m)
    rows=[]
    for p,r in [(3,1),(3,2),(3,3),(3,4),(3,5),(5,1),(5,2),(5,3),
                (7,1),(7,2),(7,3),(11,1),(11,2)]:
        m=p**r
        U=source(m);V=U-F(1,m)
        q,qq=U.denominator,V.denominator
        assert q%qq==0
        assert q//qq==(3 if p==3 and r%2==0 else 1)
        lo,hi=bounds(m,K)
        assert lo<pi_lo-V<pi_hi-V<hi
        oldlo,oldhi=q*q*(F(1,m)-hi),q*q*(F(1,m)-lo)
        newlo,newhi=qq*qq*lo,qq*qq*hi
        # Direct error bounds also imply this general upper bound; certify
        # it by an independently sharper interval in each actual finite row.
        assert newhi/oldlo<F(1,4*m*m-1)
        rows.append({'p':p,'r':r,'m':m,'source':str(U),'corrected':str(V),
                     'q_source':str(q),'q_corrected':str(qq),
                     'q_source_bits':q.bit_length(),'q_corrected_bits':qq.bit_length(),
                     'delta_interval':[str(lo),str(hi)],
                     'quality_source_interval':[str(oldlo),str(oldhi)],
                     'quality_corrected_interval':[str(newlo),str(newhi)],
                     'quality_source_display':display((oldlo+oldhi)/2),
                     'quality_corrected_display':display((newlo+newhi)/2),
                     'relative_quality_upper':str(newhi/oldlo),
                     'relative_quality_upper_display':display(newhi/oldlo)})
    report={'status':'PASS','K':K,'general_odd_indices_checked':len(tests),
            'primepower_rows':len(rows),'certificate_uses_pi_digits':False,
            'crosscheck_uses_independent_machin_interval':True,'rows':rows,
            'scope':'Loss-free relative acceleration; does not establish absolute quality tending to zero.'}
    Path(__file__).with_name('checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'},indent=2))
    for row in rows[:5]:
        print({k:row[k] for k in ['m','q_source_bits','q_corrected_bits','quality_corrected_display','relative_quality_upper_display']})

if __name__=='__main__':main()
