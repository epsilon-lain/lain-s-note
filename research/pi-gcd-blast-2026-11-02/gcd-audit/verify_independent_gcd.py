"""Independent finite checks supporting the all-index gcd/valuation audit.

Standard library only. The default run is intentionally small; it does not
rerun the supplied N=20000 sieve or infer an all-index bound from finite data.
"""
from fractions import Fraction as F
from math import comb, factorial, gcd, isqrt
from pathlib import Path
import json


def vp(x,p):
    if not x:
        return float('inf')
    if isinstance(x,F):
        return vp(x.numerator,p)-vp(x.denominator,p)
    out=0
    while x%p==0:
        x//=p
        out+=1
    return out


def flog(n,p):
    out=0
    while n>=p:
        n//=p
        out+=1
    return out


def fmod(x,p):
    x=F(x)
    assert x.denominator%p
    return x.numerator*pow(x.denominator,-1,p)%p


def prime_list(n):
    return [p for p in range(2,n+1) if all(p%d for d in range(2,isqrt(p)+1))]


def mm(a,b,mod):
    return [[sum(a[i][k]*b[k][j] for k in range(2))%mod
             for j in range(2)] for i in range(2)]


def mpow(a,n,mod):
    out=[[1,0],[0,1]]
    while n:
        if n&1:
            out=mm(out,a,mod)
        a=mm(a,a,mod)
        n>>=1
    return out


def poly_mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j]+=x*y
    return out


def main():
    N=300
    aa=[1,1]; bb=[0,1]; ff=[1,1]
    for n in range(1,N+1):
        aa.append((2*n+1)*aa[-1]+n*n*aa[-2])
        bb.append((2*n+1)*bb[-1]+n*n*bb[-2])
        ff.append(ff[-1]*(n+1))
    ell=[F(a,f) for a,f in zip(aa,ff)]
    ww=[F(b,f) for b,f in zip(bb,ff)]
    counts={}
    for n in range(61):
        explicit=sum(F(factorial(n),2**k*factorial(k)**2*factorial(n-2*k))
                     for k in range(n//2+1))
        assert ell[n]==explicit
    counts['independent_binomial_ell_checks']=61
    for n in range(1,81):
        assert ww[n]==sum(ell[j-1]*ell[n-j]/j for j in range(1,n+1))
    counts['convolution_checks']=80
    for n in range(1,N+1):
        assert aa[n]*bb[n-1]-aa[n-1]*bb[n]==(-1)**n*ff[n-1]**2
    counts['integer_wronskian_checks']=N
    lucas=pole=cap=strong=0
    primes=[p for p in prime_list(47) if p>2]
    for p in primes:
        chi=(-1)**((p-1)//2)
        for n in range(1,201):
            b,r=divmod(n,p)
            assert fmod(ell[n]-ell[b]*ell[r],p)==0
            lucas+=1
            h=flog(n,p)
            b,r=divmod(n,p**h)
            residue=fmod(chi**h*ell[r]*ww[b],p)
            assert fmod(p**h*ww[n],p)==residue
            if residue:
                assert vp(ww[n],p)==-h
            pole+=1
            if n<2:
                continue
            e=min(vp(ell[n],p),vp(ww[n],p))
            c=min(flog(n-1,p)-vp(n,p),flog(n+1,p)-vp(n+1,p))
            assert e<=c
            cap+=1
            if e>0:
                u=vp(ell[n],p)
                assert u<=c
                assert vp(ww[n-1],p)==-vp(n,p)-u
                assert vp(ww[n+1],p)==-vp(n+1,p)-u
                strong+=1
        block=[[1,0],[0,1]]
        for n in range(1,p+1):
            block=mm([[2*n+1,n*n],[1,0]],block,p)
        want=[[(-chi)%p,chi%p],[(-chi)%p,chi%p]]
        assert block==want and mm(block,block,p)==[[0,0],[0,0]]
    counts.update(lucas_checks=lucas,leading_pole_checks=pole,
                  two_sided_cap_checks=cap,positive_excess_neighbor_checks=strong,
                  prime_nilpotent_blocks=len(primes))
    block_checks=0
    for modulus in (3,9,25,49,11,27):
        prefixes=[[[1,0],[0,1]]]
        for n in range(1,modulus+1):
            prefixes.append(mm([[2*n+1,n*n],[1,0]],prefixes[-1],modulus))
        for n in (1,2,modulus,modulus+1,2*modulus+4):
            k,s=divmod(n-1,modulus)
            product=mm(prefixes[s],mpow(prefixes[-1],k,modulus),modulus)
            assert (product[0][0]+product[0][1])%modulus==aa[n]%modulus
            assert (product[1][0]+product[1][1])%modulus==aa[n-1]%modulus
            assert product[0][0]==bb[n]%modulus
            assert product[1][0]==bb[n-1]%modulus
            block_checks+=1
    counts['prime_power_transfer_checks']=block_checks
    assert (vp(aa[72],17),vp(bb[72],17),vp(ff[72],17))==(5,5,4)
    first_odd=None
    for n in range(1,73):
        gain=F(gcd(aa[n],4*bb[n]),ff[n])
        odd_numerator=gain.numerator//2**vp(gain.numerator,2)
        if odd_numerator>1 and first_odd is None:
            first_odd=(n,odd_numerator)
    assert first_odd==(72,17)
    repeated=[]
    for end in (4,72,1228,20880):
        baseline=sum(end//17**k for k in range(1,flog(end,17)+1))
        modulus=17**(baseline+4)
        ap,a,bp,b=1,1,0,1
        for n in range(1,end):
            ap,a=a,((2*n+1)*a+n*n*ap)%modulus
            bp,b=b,((2*n+1)*b+n*n*bp)%modulus
        va,vb=vp(a,17)-baseline,vp(b,17)-baseline
        repeated.append([end,va,vb,min(va,vb)])
    assert repeated==[[4,1,0,0],[72,1,1,1],[1228,2,-1,-1],[20880,2,0,0]]

    # Independently derive the U-only upper gain from exact A,B for a short range.
    up,u=1,2
    lcm=1
    survivors=[]
    exact_rows=[]
    for n in range(1,N+1):
        if n>1:
            up,u=u,(2*(2*n-1)*u+4*(n-1)*up)//n
            lcm=lcm*(n-1)//gcd(lcm,n-1)
        assert F(u,2**n)==ell[n]
        E=lcm//gcd(lcm,n)
        odd_u=u//2**vp(u,2)
        odd_E=E//2**vp(E,2)
        G=gcd(odd_u,odd_E)
        gain=F(gcd(aa[n],4*bb[n]),ff[n])
        upper_gain=F(2)**(vp(u,2)-n)*G
        assert gain<=upper_gain
        floor=F(2,2*n+1)/upper_gain**2
        if floor<F(1,5):
            survivors.append(n)
            truefloor=F(2,2*n+1)/gain**2
            exact_rows.append([n,str(truefloor),str(2*truefloor)])
    assert survivors==[2,3,6,7,15]
    assert [n for n,lo,hi in exact_rows if F(lo)<F(1,5)]==[2]
    assert exact_rows[0]==[2,'1/10','1/5']
    counts['sieve_bound_vs_actual_gcd_checks']=N

    # Formal inverse-binomial ODE and duplication identity, exactly in Q[z].
    order=40
    c=[F(factorial(j)**2,factorial(2*j+1)) for j in range(order+1)]
    for j in range(1,order+1):
        assert (4*j+2)*c[j]==j*c[j-1]
    compose=[F(0)]*(order+1)
    for j in range(order+1):
        for k in range(min(j,order-j)+1):
            compose[j+k]+=c[j]*comb(j,k)*4**(j-k)*(-1)**k
    for k in range(order+1):
        assert 2*compose[k]-(compose[k-1] if k else 0)==2*c[k]
    counts['formal_ode_coefficients']=order
    counts['formal_duplication_coefficients']=order+1
    H=[F(0)]*10
    P=[F(1)]
    for j in range(1,6):
        P=poly_mul(P,[F(1),F(-1,j)])
        square=poly_mul(P,P)
        square[0]-=1
        assert square[0]==0
        t=F(2**j*factorial(j)**2,factorial(2*j+1))
        for k,value in enumerate(square[1:]):
            H[k]+=t*value
    H=[value*2**k for k,value in enumerate(H)]
    H[0]-=2
    W=[value*F(155925,4) for value in H]
    want=[-133890,106091,-112150,95560,-57235,25853,-8200,1800,-240,16]
    assert W==list(map(F,want))
    for n in range(N+1):
        m,r=divmod(n,4)
        expected=(2*m,2*m,2*m+2,2*m+3+vp(m+1,2))[r]
        assert vp(aa[n],2)==expected
    counts['two_adic_recurrence_checks']=N+1
    for j in range(6,501):
        assert vp(factorial(j),2)-flog(j,2)>=2
    counts['two_adic_tail_bound_checks']=495
    result={'status':'passed','small_range_limit':N,'counts':counts,
            'first_odd_excess_through_72':list(first_odd),
            'n72_p17_A_B_factorial_valuations':[5,5,4],
            'n72_p17_first_units':[str((aa[72]//17**5)%17),str((bb[72]//17**5)%17)],
            'repeated_base17_digit_checks':repeated,
            'small_range_sieve_survivors':survivors,'exact_survivor_floors_and_ceilings':exact_rows,
            'H5_polynomial_ascending':want,
            'scope':'Independent finite checks plus proof audit; full N=20000 sieve was code-reviewed but not redundantly rerun. No all-n conclusion is inferred from the checked range.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
