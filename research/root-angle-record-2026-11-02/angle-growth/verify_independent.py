"""Independent integer/Fraction verification; standard library only."""
from pathlib import Path
from fractions import Fraction as F
from math import gcd
import json

HERE=Path(__file__).resolve().parent

def atan_bounds(base,n):
    s=F(0)
    for j in range(n): s+=F((-1)**j,(2*j+1)*base**(2*j+1))
    t=F((-1)**n,(2*n+1)*base**(2*n+1))
    return min(s,s+t),max(s,s+t)

def pi_bounds():
    a,b=atan_bounds(5,220);c,d=atan_bounds(239,75)
    return 16*a-4*d,16*b-4*c

def unpack(v):return F(int(v[0]),int(v[1]))
def pack(v):return [str(v.numerator),str(v.denominator)]
def score(v,L,U):
    a,b=v-L,v-U
    assert a*b>0
    return v.denominator**2*min(abs(a),abs(b)),v.denominator**2*max(abs(a),abs(b))

def add(a,b):
    out=[F(0)]*max(len(a),len(b))
    for i,v in enumerate(a):out[i]+=v
    for i,v in enumerate(b):out[i]+=v
    return trim(out)
def scale(a,c):return trim([v*c for v in a])
def trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a
def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return trim(out)
def power(a,n):
    out=[F(1)]
    for _ in range(n):out=mul(out,a)
    return out
def derivative(a):return [F(0)] if len(a)==1 else [j*a[j] for j in range(1,len(a))]
def value(a,x):
    out=F(0)
    for c in reversed(a):out=out*x+c
    return out

def verify_polynomials():
    B=[F(0),F(1),F(-1)]
    P=mul([F(1),F(-2)],add([F(1)],scale(B,-1)))
    Q=mul(B,add([F(1)],scale(B,-2)))
    D=add(power(P,2),power(Q,2))
    N0=add(mul(P,derivative(Q)),scale(mul(Q,derivative(P)),-1))
    NJ=add(add([F(5)],scale(B,-30)),add(scale(power(B,2),40),scale(power(B,3),-20)))
    N=scale(power(P,2),4)
    assert add(N0,scale(NJ,F(3,5)))==N
    # A sum-of-squares form supplies a rigorous denominator lower bound 1/64.
    V=[F(1,4),F(-1),F(1)]
    positive=[F(1,64)]
    for exponent,coefficient in [(1,F(9,4)),(2,F(11,2)),(3,F(4)),(4,F(4))]:
        positive=add(positive,scale(power(V,exponent),coefficient))
    assert positive==D
    A=N; endpoint_differences={}
    for j in range(8):
        if j in (1,3,5,7):endpoint_differences[j]=value(A,F(1))-value(A,F(0))
        if j==6:
            L1=sum(abs(c) for c in A)
            C6=F(64**7)*L1/F(15360)
        A=add(mul(D,derivative(A)),scale(mul(derivative(D),A),-(j+1)))
    assert endpoint_differences=={1:F(0),3:F(0),5:F(5760),7:F(645120)}
    assert C6==1077181949922975036148285440
    return {'polynomial_identities':True,'denominator_lower_bound':pack(F(1,64)),
            'odd_endpoint_differences':{str(k):pack(v) for k,v in endpoint_differences.items()},
            'sixth_derivative_numerator_l1':pack(L1),
            'all_node_count_error_bound_C6':pack(C6),
            'bound':'For all N>=1, |T_N(K*)-pi| <= C6/N^6; conservative derivative bound, not a quality-zero result.'}

def verify_rotations(L,U):
    count=0; full=0
    for q in range(1,26):
        for p in range(3*q,4*q+1):
            if gcd(p,q)!=1:continue
            N=p*p+q*q
            for a in range(-3,4):
                for b in range(1,9):
                    if gcd(a,b)!=1:continue
                    P,Q=b*p+a*q,b*q-a*p
                    G=gcd(P,Q)
                    assert gcd(N,a*a+b*b)%G==0
                    count+=1
            b0=(-p*pow(q,-1,N))%N
            for k in [0,1,2,16,293]:
                b=b0+N*k
                if b==0:continue
                P,Q=b*p-q,b*q+p
                assert gcd(P,Q)==N
                assert p*(Q//N)-q*(P//N)==1
                full+=1
    stored=json.loads((HERE/'rotation/angle_chain_results.json').read_text())
    for r in stored['rounds']:
        p,q=int(r['input']['p']),int(r['input']['q'])
        c,d=int(r['canonical_mate']['c']),int(r['canonical_mate']['d'])
        k=r['selection']['selected_k'];P,Q=p*k+c,q*k+d
        assert (str(P),str(Q))==(r['output']['p'],r['output']['q'])
        e0,e1=p-q*U,p-q*L
        assert e0>0 and c-d*L<0
        z0=(L*d-c)/(p-L*q);z1=(U*d-c)/(p-U*q)
        assert z0.numerator//z0.denominator==z1.numerator//z1.denominator==r['selection']['certified_floor_z']
        candidates=[0,z0.numerator//z0.denominator,z0.numerator//z0.denominator+1]
        bounds={j:score(F(p*j+c,q*j+d),L,U) for j in candidates}
        assert all(j==k or bounds[k][1]<bounds[j][0] for j in candidates)
        for alpha in (L,U,(L+U)/2):
            e=p-alpha*q;z=(alpha*d-c)/e;eta=k-z
            assert Q*abs(P-alpha*Q)==abs(eta)*(1+q*e*eta)
        v0=score(F(p,q),L,U);v1=score(F(P,Q),L,U)
        assert (v1[1]<v0[0])==r['quality_strictly_improves']
        assert (v1[0]>v0[1])==r['quality_strictly_worsens']
    # Periodic polynomial multiplier: ordinary rational clock, two real coordinates.
    for n in range(-3,4):
        for t in (F(n),F(2*n+1,2)):
            u=t-(t.numerator//t.denominator)
            H=16*u*u*(1-u)*(1-u)
            assert H==(0 if t.denominator==1 else 1)
            real,imag=F(7)+F(22)*H/8601,F(22)-F(7)*H/8601
            if H==0:assert (real,imag)==(7,22)
            else:assert imag/real==F(355,113)
    return {'general_gcd_divisibility_cases':count,'full_norm_cancellation_cases':full,
            'two_round_optimal_selection_independently_rechecked':True,
            'quality_transfer_exactly_checked':True,'periodic_growth_multiplier_checked':True}

def verify_kernel_candidates(L,U):
    data=json.loads((HERE/'angle_kernel/trace_angle_results.json').read_text())
    def trace(u):
        B=u*(1-u);P=(1-2*u)*(1-B);Q=B*(1-2*B)
        # Derive W directly from spinor coordinates instead of stored expanded numerator.
        bp=1-2*u
        pp=-2*(1-B)-(1-2*u)*bp
        qp=bp*(1-4*B)
        return (P*qp-Q*pp)/(P*P+Q*Q)
    def star(u):
        B=u*(1-u);P=(1-2*u)*(1-B);Q=B*(1-2*B)
        return 4*P*P/(P*P+Q*Q)
    def base(u):return 4/(1+u*u)
    cache={}
    def T(fun,n):
        key=(fun.__name__,n)
        if key not in cache:cache[key]=((fun(F(0))+fun(F(1)))/2+sum((fun(F(j,n)) for j in range(1,n)),F(0)))/n
        return cache[key]
    for row in data['rows']:
        n=row['n'];m=row['method']
        if m=='trace_trapezoid':v=T(trace,n)
        elif m=='flat6_trapezoid':v=T(star,n)
        elif m=='standard_trapezoid':v=T(base,n)
        elif m=='standard_EM1':v=T(base,n)+F(1,6*n*n)
        elif m=='trace_Simpson':v=(4*T(trace,n)-T(trace,n//2))/3
        elif m=='standard_Simpson':v=(4*T(base,n)-T(base,n//2))/3
        elif m=='trace_Richardson4':v=(16*T(trace,n)-T(trace,n//2))/15
        elif m=='flat6_Richardson6':v=(64*T(star,n)-T(star,n//2))/63
        else:raise AssertionError(m)
        assert v==unpack(row['fraction'])
        a,b=score(v,L,U)
        assert unpack(row['quality_lower'])<=a<=b<=unpack(row['quality_upper'])
    ratios={}
    for n in [21,26]:
        a=T(star,n).denominator;b=T(trace,n).denominator
        ratios[str(n)]=pack(F(a,b))
    assert ratios=={'21':pack(F(7)),'26':pack(F(13,2))}
    return {'all_candidates_rebuilt_from_spinor':len(data['rows']),
            'all_quality_bounds_independently_rechecked':True,
            'actual_denominator_growth_counterexamples':ratios}

def main():
    L,U=pi_bounds()
    out={'polynomials':verify_polynomials(),'rotations':verify_rotations(L,U),
         'finite_kernel_candidates':verify_kernel_candidates(L,U),
         'scope':'Finite exact audits plus proved algebraic identities; no new pi record, no infinite quality-zero construction.'}
    (HERE/'independent_verification.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
