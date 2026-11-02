"""Exact multi-frequency endpoint and root-branch checks; standard library."""
from fractions import Fraction as F
from math import comb, factorial, gcd, lcm, prod
from pathlib import Path
import json


def det(matrix):
    a=[list(map(F,row)) for row in matrix]
    out=F(1)
    for i in range(len(a)):
        p=next((j for j in range(i,len(a)) if a[j][i]),None)
        if p is None:return F(0)
        if p!=i:a[i],a[p]=a[p],a[i];out=-out
        pivot=a[i][i];out*=pivot
        for j in range(i+1,len(a)):
            ratio=a[j][i]/pivot
            for k in range(i+1,len(a)):a[j][k]-=ratio*a[i][k]
    return out


def A(n):return prod(factorial(k) for k in range(n))


def cols(n,profile):
    return [(0,a) for a in range(n)]+[(h,a) for h,r in profile for a in range(r)]


def raw(n,profile,z,ys):
    return det([[ys.get(h,F(1))**j*sum(F(comb(a,k))*(j*z)**(a-k)*
                    F(h**(s-k),factorial(s-k)) for k in range(min(s,a)+1))
                for h,a in cols(n,profile)] for j in (0,1) for s in range(n)])


def moment(k,w,y):
    v=(y-1)/w
    for j in range(1,k+1):v=y/w-j*v/w
    return v


def moment_formula(n,profile,z,ys):
    if not z:return F(0)
    exponent=n*(n+1)//2+sum(r*(r-1)//2 for h,r in profile)
    factor=prod(h**(n*r) for h,r in profile)
    m=[[moment(i+a,h*z,ys[h]) for h,r in profile for a in range(r)] for i in range(n)]
    return F(factor,A(n)**2)*z**exponent*det(m)


def evaluate(poly,x):
    value=F(0)
    for c in reversed(poly):value=value*x+c
    return value


def interpolate(values):
    out=[F(0)]*len(values);basis=[F(1)];differences=list(values)
    for k in range(len(values)):
        for j,c in enumerate(basis):out[j]+=c*differences[0]
        differences=[b-a for a,b in zip(differences,differences[1:])]
        nxt=[F(0)]*(len(basis)+1)
        for j,c in enumerate(basis):nxt[j]-=c*F(k,k+1);nxt[j+1]+=c/F(k+1)
        basis=nxt
    while len(out)>1 and not out[-1]:out.pop()
    return out


def trim(p):
    while len(p)>1 and not p[-1]:p.pop()
    return p


def remainder(p,q):
    p=p[:]
    while len(p)>=len(q) and any(p):
        k=len(p)-len(q);c=p[-1]/q[-1]
        for j,v in enumerate(q):p[j+k]-=c*v
        trim(p)
    return p


def sturm(poly):
    seq=[trim(poly[:]),trim([k*c for k,c in enumerate(poly)][1:])]
    while any(seq[-1]):
        r=[-c for c in remainder(seq[-2],seq[-1])]
        if not any(r):break
        seq.append(r)
    return seq


def variations(seq,x):
    signs=[1 if v>0 else -1 for p in seq if (v:=evaluate(p,x))]
    return sum(a!=b for a,b in zip(signs,signs[1:]))


def count_roots(poly,lo,hi):
    assert evaluate(poly,lo) and evaluate(poly,hi)
    seq=sturm(poly)
    return variations(seq,lo)-variations(seq,hi)


def atan_interval(q,N=65):
    v=sum(F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(N))
    w=v+F((-1)**N,(2*N+1)*q**(2*N+1))
    return min(v,w),max(v,w)


def period_interval():
    a,b=atan_interval(5),atan_interval(239)
    return 32*a[0]-8*b[1],32*a[1]-8*b[0]


def interval_eval(poly,x):
    y=(F(0),F(0))
    for c in reversed(poly):
        products=[a*b for a in y for b in x]
        y=min(products)+c,max(products)+c
    return y


def absolute_interval(poly,x):
    lo,hi=interval_eval(poly,x)
    assert lo*hi>0
    return (lo,hi) if lo>0 else (-hi,-lo)


def integer_enclosure(interval):
    lo,hi=interval
    return [str(lo.numerator//lo.denominator),str(-(-hi.numerator//hi.denominator))]


def main():
    profiles=[(2,[(1,1),(2,1)]),(3,[(1,2),(2,1)]),
              (4,[(1,2),(2,2)]),(4,[(1,2),(2,1),(3,1)]),
              (5,[(1,3),(2,2)]),(6,[(1,3),(2,3)]),(8,[(1,4),(2,4)])]
    identity_checks=0;packets=[]
    for n,profile in profiles:
        assert sum(r for h,r in profile)==n
        H=max(2*n+max(h for h,r in profile),5)
        assert F(H,n)>1 and F(4,H)<1 and F(4*n,H*H)<1
        assert [(s,beta) for s in range(n+1) for beta in range(2)
                if F(H,n)*s+H*beta<H]==[(s,0) for s in range(n)]
        assert len(cols(n,profile))==2*n
        assert all(h+2*a<=H for h,a in cols(n,profile))
        for z in (F(1),F(-2),F(3,2)):
            for option in range(3):
                ys={h: (F(1) if option==0 else F(h+1,2) if option==1 else F((-1)**h*(h+1))) for h,r in profile}
                assert raw(n,profile,z,ys)==moment_formula(n,profile,z,ys)
                identity_checks+=1
        packets.append({'n':n,'profile':profile,'height':H,'source_rows':2*n})

    fixture=json.loads(Path(__file__).with_name('old_period_branch_fixture.json').read_text())
    old={c['n']:c for c in fixture['cases']}
    period=period_interval();assert 6<period[0]<period[1]<F(44,7)
    brackets={4:(F('3.87298'),F('3.87299')),
              6:(F('6.26801'),F('6.26803')),
              8:(F('6.28317604'),F('6.28317605'))}
    records=[]
    for n in (4,6,8):
        profile=[(1,n//2),(2,n//2)];ys={1:F(1),2:F(1)}
        assert all(h+2*a<=2*n for h,a in cols(n,profile))
        bound=n*(n-1)//2+2*(n//2)*(n//2-1)//2
        polynomial=interpolate([moment_formula(n,profile,F(z),ys) for z in range(bound+1)])
        for z in (F(-1),F(1,3),F(bound+2)):
            assert evaluate(polynomial,z)==raw(n,profile,z,ys)
        assert all(not c for k,c in enumerate(polynomial) if k<n or k%2)
        denominator=lcm(*(c.denominator for c in polynomial))
        primitive=[int(c*denominator) for c in polynomial]
        assert gcd(*primitive)==1
        imag=[c*(-1)**(k//2) if k%2==0 else F(0) for k,c in enumerate(polynomial)]
        reduced_imag=imag[n:]
        lo,hi=brackets[n]
        assert count_roots(reduced_imag,lo,hi)==1
        assert evaluate(reduced_imag,lo)*evaluate(reduced_imag,hi)<0
        for _ in range(90):
            mid=(lo+hi)/2
            if evaluate(reduced_imag,lo)*evaluate(reduced_imag,mid)<0:hi=mid
            else:lo=mid
        oldpoly={int(k):F(v) for k,v in old[n]['polynomial'].items()}
        oi=[oldpoly.get(k,F(0))*(-1)**(k//2) if k%2==0 else F(0) for k in range(max(oldpoly)+1)]
        oldlo,oldhi=map(F,old[n]['root_imaginary_part_interval'])
        if oldlo==oldhi:assert evaluate(oi,oldlo)==0
        else:assert count_roots(oi[n:],oldlo,oldhi)==1
        assert hi<oldlo<=oldhi<period[0]
        newabs=absolute_interval(imag,period)
        exponent=n*(n+1)//2+sum(r*(r-1)//2 for h,r in profile)
        factor=prod(h**(n*r) for h,r in profile)
        gram_product=F(A(n)**4,A(2*n))
        for h,r in profile:gram_product*=F(A(r)**4,A(2*r))
        bound_squared=F(factor*factor,A(n)**4)*F(44,7)**(2*exponent)*gram_product
        assert newabs[1]**2<bound_squared
        oldabs=absolute_interval(oi,period)
        oldden=lcm(*(c.denominator for c in oldpoly.values()))
        primitive_abs=tuple(v*denominator for v in newabs)
        oldprimitive_abs=tuple(v*oldden for v in oldabs)
        comparison='smaller' if primitive_abs[1]<oldprimitive_abs[0] else 'larger' if primitive_abs[0]>oldprimitive_abs[1] else 'overlap'
        records.append({'n':n,'polynomial':{str(k):str(c) for k,c in enumerate(polynomial) if c},
                        'degree':len(polynomial)-1,'least_coefficient_denominator':str(denominator),
                        'primitive_leading_coefficient':str(primitive[-1]),
                        'same_original_height_2n_budget_certified':True,
                        'gram_hadamard_majorant_certified':True,
                        'positive_imaginary_root_interval':[str(lo),str(hi)],
                        'branch_farther_from_first_period_than_old_same_row_family':True,
                        'primitive_actual_value_absolute_integer_enclosure':integer_enclosure(primitive_abs),
                        'old_primitive_actual_value_absolute_integer_enclosure':integer_enclosure(oldprimitive_abs),
                        'primitive_actual_value_comparison':comparison})
    assert records[0]['polynomial']=={'4':'-5/3','6':'-1/9'}
    assert moment_formula(2,[(1,1),(2,1)],F(3),{1:F(1),2:F(1)})==0
    output={'status':'passed','complete_packets':packets,'independent_two_parameter_identity_checks':identity_checks,
            'balanced_comparisons':records,
            'scope':'General multiweight moment and Gram identities are proved in the note. Finite branches and primitive values are exact certificates only. No general convergence, best rational approximation, or pi bad-approximability claim.'}
    Path(__file__).with_name('multiweight_verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'status':'passed','complete_packets':len(packets),'moment_identity_checks':identity_checks,'balanced_cases':len(records),'comparisons':[{k:c[k] for k in ('n','least_coefficient_denominator','primitive_actual_value_absolute_integer_enclosure','primitive_actual_value_comparison')} for c in records]},indent=2))


if __name__=='__main__':main()
