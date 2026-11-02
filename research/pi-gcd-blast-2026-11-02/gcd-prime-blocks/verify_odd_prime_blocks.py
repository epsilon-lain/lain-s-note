"""Exact finite audits of ODD_PRIME_BLOCKS_AND_EXCESS.md.

The all-n claims are proved in the note. Finite checks certify explicit
counterexamples and audit the coefficient and modular formulas.
"""
from fractions import Fraction as F
from functools import cache
from math import gcd, isqrt
from pathlib import Path
import json


def primes_to(limit):
    return [p for p in range(3,limit+1,2)
            if all(p%d for d in range(2,isqrt(p)+1))]


def valuation(x,p):
    assert x
    result=0
    while x%p==0:
        result+=1
        x//=p
    return result


def factorial_valuation(n,p):
    result=0
    while n:
        n//=p
        result+=n
    return result


def floor_log(n,p):
    assert n>=1
    h=0
    while n>=p:
        n//=p
        h+=1
    return h


def residue(x,p):
    x=F(x)
    assert x.denominator%p
    return x.numerator*pow(x.denominator,-1,p)%p


def sequences(last):
    A=[1,1]
    B=[0,1]
    fs=[1,1]
    for n in range(1,last):
        A.append((2*n+1)*A[n]+n*n*A[n-1])
        B.append((2*n+1)*B[n]+n*n*B[n-1])
        fs.append(fs[n]*(n+1))
    return A[:last+1],B[:last+1],fs[:last+1]


def matmul(A,B,m):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2))%m
                       for j in range(2)) for i in range(2))


def matpower(A,k,m):
    result=((1,0),(0,1))
    while k:
        if k&1: result=matmul(result,A,m)
        A=matmul(A,A,m)
        k//=2
    return result


def transfer(n,m):
    return ((2*n+1)%m,n*n%m),(1,0)


@cache
def blocks(m):
    prefixes=[((1,0),(0,1))]
    for n in range(1,m+1):
        prefixes.append(matmul(transfer(n,m),prefixes[-1],m))
    return prefixes,prefixes[-1]


def blocked_pair(n,m,initial):
    assert n>=1
    k,s=divmod(n-1,m)
    prefixes,S=blocks(m)
    product=matmul(prefixes[s],matpower(S,k,m),m)
    return tuple(sum(product[i][j]*initial[j] for j in range(2))%m
                 for i in range(2))


def direct_modular(n,m,initial):
    pair=initial
    for k in range(1,n):
        T=transfer(k,m)
        pair=tuple(sum(T[i][j]*pair[j] for j in range(2))%m
                   for i in range(2))
    return pair


def poly_mul(A,B,p):
    out=[0]*(len(A)+len(B)-1)
    for i,x in enumerate(A):
        for j,y in enumerate(B): out[i+j]=(out[i+j]+x*y)%p
    return out


def main():
    A,B,fs=sequences(800)
    ell=[F(a,f) for a,f in zip(A,fs)]
    w=[F(b,f) for b,f in zip(B,fs)]
    ps=primes_to(47)
    counts={'w_convolution':0,'lucas_digits':0,'leading_pole':0,
            'generic_exact_deficit':0,'backward_forward_excess_bound':0,
            'frobenius_polynomial':0,'nilpotent_prime_block':0,
            'prime_power_block_pairs':0,'counterexample_modular_certificates':0,
            'positive_excess_neighbor_valuation_identities':0}
    for n in range(1,51):
        assert w[n]==sum(ell[j-1]*ell[n-j]/j for j in range(1,n+1))
        counts['w_convolution']+=1
    for p in ps:
        poly=[1]
        for _ in range((p-1)//2): poly=poly_mul(poly,[1,-2,-1],p)
        assert len(poly)==p
        assert all(x==residue(ell[k],p) for k,x in enumerate(poly))
        counts['frobenius_polynomial']+=1
        chi=1 if ((p-1)//2)%2==0 else -1
        expected=tuple(tuple((-chi*x)%p for x in row)
                       for row in ((1,-1),(1,-1)))
        _,S=blocks(p)
        assert S==expected
        assert matmul(S,S,p)==((0,0),(0,0))
        counts['nilpotent_prime_block']+=1
        for n in range(1,501):
            digits=[]
            number=n
            while number:
                number,d=divmod(number,p)
                digits.append(d)
            product=1
            for d in digits: product=product*residue(ell[d],p)%p
            assert residue(ell[n],p)==product
            counts['lucas_digits']+=1
            h=floor_log(n,p)
            power=p**h
            b,r=divmod(n,power)
            left=residue(power*w[n],p)
            right=(pow(chi,h,p)*residue(ell[r],p)*residue(w[b],p))%p
            assert left==right
            counts['leading_pole']+=1
            e=valuation(gcd(A[n],4*B[n]),p)-factorial_valuation(n,p)
            if right:
                assert valuation(w[n].numerator,p)-valuation(w[n].denominator,p)==-h
                assert e==-h
                counts['generic_exact_deficit']+=1
            if n>=2:
                upper=min(floor_log(n-1,p)-valuation(n,p),
                          floor_log(n+1,p)-valuation(n+1,p))
                assert e<=upper
                counts['backward_forward_excess_bound']+=1
        if p<=17:
            for exponent in (1,2,3):
                m=p**exponent
                if m>350: continue
                for n in sorted({1,2,p,p+1,m,m+1,2*m,2*m+17}):
                    if n>800:continue
                    for initial,X in (((1,1),A),((1,0),B)):
                        expected_pair=(X[n]%m,X[n-1]%m)
                        assert blocked_pair(n,m,initial)==expected_pair
                        assert direct_modular(n,m,initial)==expected_pair
                        counts['prime_power_block_pairs']+=1
    # Complete exact search through n=500, all potentially relevant odd primes.
    excess=[]
    all_ps=primes_to(500)
    for n in range(1,501):
        g=gcd(A[n],4*B[n])
        for p in all_ps:
            if p>n:break
            e=valuation(g,p)-factorial_valuation(n,p)
            if e>0:
                excess.append((n,p,e))
                u=valuation(A[n],p)-factorial_valuation(n,p)
                vprev=valuation(B[n-1],p)-factorial_valuation(n-1,p)
                vnext=valuation(B[n+1],p)-factorial_valuation(n+1,p)
                assert vprev==-valuation(n,p)-u
                assert vnext==-valuation(n+1,p)-u
                cap=min(floor_log(n-1,p)-valuation(n,p),
                        floor_log(n+1,p)-valuation(n+1,p))
                assert u<=cap
                counts['positive_excess_neighbor_valuation_identities']+=1
    assert excess[0]==(72,17,1)
    certificate_pairs=((72,17),(107,73),(166,131),(170,41),(176,23),(480,31))
    certificates=[]
    for n,p in certificate_pairs:
        va,vb=valuation(A[n],p),valuation(B[n],p)
        vf=factorial_valuation(n,p)
        g=gcd(A[n],4*B[n])
        assert valuation(g,p)==min(va,vb)>vf
        # Independent modular recurrence suffices to certify both exact powers.
        for initial,X,v in (((1,1),A,va),((1,0),B,vb)):
            m=p**(v+1)
            r=direct_modular(n,m,initial)[0]
            assert r==X[n]%m and r%(p**v)==0 and r!=0
            counts['counterexample_modular_certificates']+=1
        certificates.append({'n':n,'p':p,'v_p_A':va,'v_p_B':vb,'v_p_factorial':vf,
                             'v_p_g':min(va,vb),'excess':min(va,vb)-vf,
                             'first_nonzero_A_residue':(A[n]//p**va)%p,
                             'first_nonzero_B_residue':(B[n]//p**vb)%p})
    # A specific repeated-zero-digit extrapolation fails at the next lengths.
    wanted={4,72,1228,20880}
    aa0,aa=1,1
    bb0,bb=0,1
    branch=[]
    for n in range(1,max(wanted)+1):
        if n>1:
            aa0,aa=aa,(2*n-1)*aa+(n-1)**2*aa0
            bb0,bb=bb,(2*n-1)*bb+(n-1)**2*bb0
        if n in wanted:
            vf=factorial_valuation(n,17)
            va,vb=valuation(aa,17)-vf,valuation(bb,17)-vf
            branch.append({'n':n,'v_17_ell':va,'v_17_w':vb,'excess':min(va,vb)})
    assert [(r['v_17_ell'],r['v_17_w'],r['excess']) for r in branch]==[
        (1,0,0),(1,1,1),(2,-1,-1),(2,0,0)]
    result={'status':'passed','arithmetic':'exact integers and fractions; Python standard library',
            'scope':'All-n bounds and congruences are proved analytically. Finite data certify counterexamples; no sparse-subsequence or pi bad-approximability conclusion.',
            'counts':counts,
            'odd_excesses_through_n_500':[{'n':n,'p':p,'excess':e} for n,p,e in excess],
            'selected_exact_valuation_certificates':certificates,
            'repeated_base17_digit4_branch':branch}
    Path(__file__).with_name('odd_prime_blocks_verification.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
