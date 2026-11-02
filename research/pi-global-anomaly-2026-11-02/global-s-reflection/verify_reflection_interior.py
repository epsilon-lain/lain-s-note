"""Exact reflection/folding/edge-capacity checks, Python standard library only."""
from pathlib import Path
from math import gcd
import json


def primes_up_to(n):
    flags=bytearray(b'\x01')*(n+1)
    flags[:2]=b'\x00\x00'
    for p in range(2,n+1):
        if flags[p]:
            for k in range(p*p,n+1,p):flags[k]=0
    return [p for p in range(3,n+1) if flags[p]]


def v2(n):
    assert n>0
    return (n&-n).bit_length()-1


def quadratic_product(x,y):
    a,b=x;c,d=y
    return a*c+2*b*d,a*d+b*c


def quadratic_power(n):
    result=(1,0);base=(2,1)
    while n:
        if n&1:result=quadratic_product(result,base)
        base=quadratic_product(base,base);n//=2
    return result


def bounded_by_quadratic_integer(value,pair):
    a,b=pair
    assert b>=0
    return value<=a or (value-a)**2<=2*b*b


def main():
    prime_limit=2000;index_limit=1000
    primes=primes_up_to(prime_limit)
    events={n:[] for n in range(1,index_limit+1)}
    residue_checks=folding_checks=0
    examples=[]
    for p in primes:
        ell=[1,1]
        for r in range(1,p-1):
            ell.append(((2*r+1)*ell[-1]+r*ell[-2])*pow(r+1,-1,p)%p)
        for r in range(p):
            assert ell[p-1-r]%p==((-1)**((p-1)//2+r)*ell[r])%p
            residue_checks+=1
        for n in range(p+1,2*p):
            r=n-p;k=min(r,p-1-r)
            assert 0<=3*k<=n-1
            assert (ell[r]%p==0)==(ell[k]%p==0)
            assert p in (n-k,(n+1+k)//2)
            folding_checks+=1
            if n<=index_limit and ell[r]%p==0:events[n].append((p,k))
        if p in (17,23,37,41,107):
            zeros=[r for r in range(p) if ell[r]%p==0]
            examples.append({'p':p,'zero_degrees':zeros})
    U=[1,2]
    for k in range(1,index_limit):
        numerator=2*(2*k+1)*U[-1]+4*k*U[-2]
        assert numerator%(k+1)==0
        U.append(numerator//(k+1))
    C=[u//2**v2(u) for u in U]
    height_checks=0
    for k,(u,c) in enumerate(zip(U,C)):
        assert 2*v2(u)>=k
        assert bounded_by_quadratic_integer(c,quadratic_power(k))
        height_checks+=1
    edge_checks=anomalies=0
    for n,entries in events.items():
        K=0
        while (K+1)**3<=n:K+=1
        capacity=1
        for c in C[:K+1]:capacity*=c
        edge_product=1
        for p,k in entries:
            anomalies+=1
            assert C[k]%p==0
            assert bounded_by_quadratic_integer(p,quadratic_power(k))
            if k<=K:edge_product*=p
        assert capacity%edge_product==0
        assert bounded_by_quadratic_integer(edge_product,quadratic_power(K*(K+1)//2))
        edge_checks+=1
    byp={x['p']:x['zero_degrees'] for x in examples}
    assert all(r in byp[p] for p,r in ((17,4),(107,11),(37,5),(23,7)))
    assert 6 in byp[41] and 13 in byp[41]
    result={'status':'passed','prime_limit':prime_limit,'index_limit':index_limit,
            'full_reflection_residue_checks':residue_checks,
            'folding_equivalence_checks':folding_checks,
            'odd_core_height_checks':height_checks,
            'edge_product_divisibility_checks':edge_checks,
            'anomaly_events_in_finite_index_range':anomalies,
            'residue_class_and_parity_examples':examples,
            'scope':'Finite exact checks of the displayed structural identities. No uniform moving-diagonal S(n) bound is inferred.'}
    Path(__file__).with_name('reflection_interior_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
