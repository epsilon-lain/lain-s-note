"""Exact finite audits for PRIME_BLOCKS_AND_LIFTING.md.

The all-index statements are proved in the note. These checks use no pi,
floating point, third-party library, or numerical inference.
"""
from fractions import Fraction as F
from math import comb, factorial, gcd
from pathlib import Path
import json


def primes_upto(N):
    flags = [True]*(N+1)
    flags[:2] = [False,False]
    for p in range(2,int(N**0.5)+1):
        if flags[p]:
            flags[p*p:N+1:p] = [False]*((N-p*p)//p+1)
    return [p for p in range(2,N+1) if flags[p]]


def valuation(x,p):
    assert x != 0
    n = 0
    while x%p == 0:
        x //= p
        n += 1
    return n


def residue(x,p):
    x = F(x)
    assert x.denominator%p != 0
    return x.numerator*pow(x.denominator,-1,p)%p


def matmul(a,b,p):
    return [[sum(a[i][k]*b[k][j] for k in range(2))%p
             for j in range(2)] for i in range(2)]


def main():
    N = 600
    A = [1,1]
    B = [0,1]
    D = [0,2]
    for n in range(1,N):
        A.append((2*n+1)*A[-1]+n*n*A[-2])
        B.append((2*n+1)*B[-1]+n*n*B[-2])
        D.append((2*n+1)*D[-1]+n*n*D[-2]+2*A[-2]+2*n*A[-3])
    fac = [factorial(n) for n in range(N+1)]
    a = [F(A[n],fac[n]) for n in range(N+1)]
    b = [F(B[n],fac[n]) for n in range(N+1)]
    g = [gcd(A[n],4*B[n]) for n in range(N+1)]
    counts = {'constant_term':0,'convolution':0,'Lucas':0,
              'reflection':0,'leading_denominator_residue':0,
              'exact_generic_valuation':0,'exceptional_valuation_floor':0,
              'transfer_blocks':0,'top_half_prime_square':0,
              'annihilation_tail':0}
    for n in range(51):
        ct = sum(F(comb(n,2*j)*comb(2*j,j),2**j) for j in range(n//2+1))
        assert a[n] == ct
        counts['constant_term'] += 1
    for n in range(1,101):
        assert b[n] == sum(a[n-j]*a[j-1]/j for j in range(1,n+1))
        counts['convolution'] += 1
    for p in primes_upto(199):
        if p == 2:
            continue
        s = (-1)**((p-1)//2)
        R = [[1,0],[0,1]]
        for k in range(1,p):
            R = matmul([[2*k+1,k*k],[1,0]],R,p)
        assert R == [[-s%p,s%p],[0,-s%p]]
        assert matmul(R,[[1,0],[1,0]],p) == [[0,0],[-s%p,0]]
        counts['transfer_blocks'] += 2
        for r in range(p):
            assert residue(a[p-1-r],p) == s*(-1)**r*residue(a[r],p)%p
            counts['reflection'] += 1
        for n in range(N+1):
            m,r = divmod(n,p)
            assert residue(a[n],p) == residue(a[m],p)*residue(a[r],p)%p
            counts['Lucas'] += 1
        for n in range(p,N+1):
            e = 0
            pe = 1
            while pe*p <= n:
                e += 1
                pe *= p
            top,r = divmod(n,pe)
            lhs = residue(pe*b[n],p)
            rhs = pow(s,e,p)*residue(a[r],p)*residue(b[top],p)%p
            assert lhs == rhs
            counts['leading_denominator_residue'] += 1
            actual = valuation(g[n],p)
            baseline = valuation(fac[n],p)-e
            if rhs:
                assert actual == baseline
                counts['exact_generic_valuation'] += 1
            else:
                assert actual >= baseline+1
                counts['exceptional_valuation_floor'] += 1
            if n >= 2*p:
                assert actual >= 1
                counts['annihilation_tail'] += 1
        for r in range(1,p):
            n = p+r
            if A[r]%p:
                predicted = 0
            else:
                predicted = 2 if (A[r]//p+D[r])%p == 0 else 1
            assert valuation(g[n],p) == predicted
            counts['top_half_prime_square'] += 1
    excess = []
    for n in range(1,501):
        z = g[n]//gcd(g[n],fac[n])
        while z%2 == 0:
            z //= 2
        if z > 1:
            excess.append({'n':n,'odd_excess_factor':str(z)})
    examples = []
    for p,r in [(17,4),(73,34),(131,35)]:
        n = p+r
        examples.append({'p':p,'r':r,'n':n,
                         'A_r_over_p_residue':A[r]//p%p,
                         'D_r_residue':D[r]%p,
                         'g_valuation':valuation(g[n],p),
                         'factorial_valuation':valuation(fac[n],p)})
    result = {'status':'PASS','max_sequence_index':N,'counts':counts,
              'prime_square_examples':examples,
              'odd_excess_examples_through_500':excess,
              'scope':'Exact finite audits only; the all-index identities require the proofs in the note. No aggregate gcd growth conclusion.'}
    Path(__file__).with_name('prime_block_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
