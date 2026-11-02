"""Exact modular search for simultaneous odd-prime factorial surplus.

For each p, a fixed modulus p**(v_p(N!)+3) preserves every tested condition
p**(v_p(n!)+1) | gcd(A_n,B_n), n<=N. No floating point or heuristic inference.
"""
from math import gcd, factorial, log
from pathlib import Path
from time import perf_counter
import json

N = 50000
PRIME_BOUND = 500


def primes_upto(n):
    return [p for p in range(3,n+1,2)
            if all(p%d for d in range(3,int(p**0.5)+1,2))]


def vp_fact(n,p):
    v = 0
    while n:
        n //= p
        v += n
    return v


def vp(n,p):
    v = 0
    while n%p == 0:
        n //= p
        v += 1
    return v


def digit_zero_set(p):
    a0=a1=1
    roots=[]
    for n in range(1,p-1):
        a0,a1=a1,((2*n+1)*a1+n*n*a0)%p
        if a1==0:
            roots.append(n+1)
    return roots


def scan(p):
    started=perf_counter()
    modulus=p**(vp_fact(N,p)+3)
    threshold=p
    a0=a1=1
    b0=0
    b1=1
    hits=[]
    for n in range(1,N+1):
        if n>=2:
            k=n-1
            a0,a1=a1,((2*k+1)*a1+k*k*a0)%modulus
            b0,b1=b1,((2*k+1)*b1+k*k*b0)%modulus
        k=n
        while k%p==0:
            threshold*=p
            k//=p
        if a1%threshold==0 and b1%threshold==0:
            hits.append(n)
    return {'p':p,'hits':hits,'seconds':round(perf_counter()-started,3)}


def certify(indices):
    if not indices:
        return []
    targets=set(indices)
    out=[]
    a0=a1=1
    b0=0
    b1=1
    fac=1
    for n in range(1,max(indices)+1):
        if n>=2:
            k=n-1
            a0,a1=a1,(2*k+1)*a1+k*k*a0
            b0,b1=b1,(2*k+1)*b1+k*k*b0
        fac*=n
        if n in targets:
            g=gcd(a1,4*b1)
            excess=g//gcd(g,fac)
            odd=excess
            while odd%2==0:
                odd//=2
            out.append({'n':n,'odd_excess_factor':str(odd),
                        'A_bits':a1.bit_length(),'g_bits':g.bit_length(),
                        'log_h':log(g)+0.5*log(n)-log(fac),
                        'surplus_valuations':[{'p':p,'v_g':vp(g,p),'v_factorial':vp_fact(n,p)}
                                              for p in primes_upto(PRIME_BOUND)
                                              if vp(g,p)>vp_fact(n,p)]})
    return out


if __name__=='__main__':
    started=perf_counter()
    eligible=[]
    skipped=[]
    for p in primes_upto(PRIME_BOUND):
        roots=digit_zero_set(p)
        if roots:
            eligible.append((p,roots))
        else:
            skipped.append(p)
    results=[]
    by_index={}
    for p,roots in eligible:
        result=scan(p)
        result['digit_roots']=roots
        results.append(result)
        for n in result['hits']:
            by_index.setdefault(n,[]).append(p)
    combined=[{'n':n,'primes':ps} for n,ps in sorted(by_index.items()) if len(ps)>=2]
    # Certify the first three simultaneous events with full unreduced integers.
    certificates=certify([r['n'] for r in combined[:3]])
    output={'status':'PASS','max_index':N,'prime_bound':PRIME_BOUND,
            'eligible_prime_count':len(eligible),
            'primes_with_no_digit_zero_skipped':skipped,
            'combined_events':combined,'full_integer_certificates':certificates,
            'prime_scans':results,'seconds':round(perf_counter()-started,3),
            'scope':'Finite exact search only; no assertion about unbounded surplus or h_n.'}
    Path(__file__).with_name('combined_surplus_checks.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='prime_scans'},indent=2))
