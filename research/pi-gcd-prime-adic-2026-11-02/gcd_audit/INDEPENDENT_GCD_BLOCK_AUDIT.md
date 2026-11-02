# Independent audit: prime blocks and prime-square lifting

2026-11-02. This note audits `../gcd_blocks/PRIME_BLOCKS_AND_LIFTING.md`. The claims are about the particular recurrence

    A_0=A_1=1, B_0=0, B_1=1,
    X_(n+1)=(2n+1)X_n+n^2 X_(n-1),
    g_n=gcd(A_n,4B_n).

The audit independently checked the all-prime arguments. No mathematical gap was found. This is an independent agent review, not external peer review.

## 1. Leading denominator residue

Put `a_n=A_n/n!`, `b_n=B_n/n!`. Their generating functions give

    a(z)=(1-2z-z^2)^(-1/2),
    b_n=sum_(j=1)^n a_(n-j)*a_(j-1)/j.

The constant-term expression for `a_n` ensures its denominator is a power of two. For every odd prime p, Frobenius in the Laurent-polynomial expression proves `a_(mp+r)=a_m*a_r (mod p)` for `0<=r<p`. In particular this Lucas property applies to all base-p digits, not just a single block.

Let `e=floor(log_p n)` and `n=N*p^e+r`, where `1<=N<p` and `0<=r<p^e`. Multiplication by `p^e` makes every term in the convolution p-integral. All terms except `j=l*p^e`, `1<=l<=N`, vanish modulo p. The surviving factors are

    a_(n-l*p^e)=a_(N-l)*a_r (mod p),
    a_(l*p^e-1)=a_(l-1)*a_(p-1)^e (mod p).

Since `a_(p-1)=s_p=(-1)^((p-1)/2) (mod p)`, summing the remaining N terms proves exactly

    p^e*b_n=s_p^e*a_r*b_N (mod p).

If the residue is nonzero, `v_p(B_n)=v_p(n!)-e`, while `v_p(A_n)>=v_p(n!)`; therefore the claimed exact gcd valuation follows. If it vanishes, the lower bound is still valid because `e>=1` gives `v_p(n!)-e+1<=v_p(n!)`. No unjustified exact valuation is asserted in that case.

## 2. Transfer matrices and singular indices

The initial fundamental matrix at index one is

    [[A_1,B_1],[A_0,B_0]]=[[1,1],[1,0]].

Wilson's theorem, the Lucas property, the leading b_p residue, and paired terms in b_(p-1) give the final matrix modulo p

    [[A_p,B_p],[A_(p-1),B_(p-1)]]=[[0,-s_p],[-s_p,0]].

Multiplying this by the inverse of the initial matrix yields exactly

    M_(p-1)...M_1 = s_p*[[-1,1],[0,-1]] (mod p).

Right multiplication by `M_0=[[1,0],[1,0]]` gives

    M_(p-1)...M_0 = [[0,0],[-s_p,0]] (mod p).

This singular matrix must not be used to advance the given B initial condition from index zero, since that would assume `B_1=B_0`. The draft correctly applies periodicity from index p. It yields `X_(2p)=0 (mod p)` for both actual solutions, and the next singular step gives `X_(2p+1)=0`; all later indices then vanish modulo p. The argument also handles p=3 without any change of the sign character.

## 3. Exact second-power lifting for p<n<2p

At the first shifted step the recurrence reads

    X_(p+1)=(2p+1)X_p+p^2 X_(p-1).

The second term vanishes modulo p^2. Subsequent recurrence steps therefore prove, by induction, `X_(p+r)=C_r(p)*X_p (mod p^2)`. This is valid for both A and B. Taylor expansion of the integer polynomial gives

    C_r(p)=A_r+p*D_r (mod p^2).

The differentiated recurrence for D_r has both forcing terms `2*A_r` and `2r*A_(r-1)`; neither is missing.

The exact Wronskian

    A_n B_(n-1)-A_(n-1) B_n=(-1)^n*((n-1)!)^2

shows `v_p(g_n)<=2` when `p<n<2p`. Also `B_p=-s_p (mod p)` is a unit. Consequently:

- If `A_r` is a unit modulo p, then `B_(p+r)` is a unit, so the gcd valuation is zero.
- If `p|A_r` and `A_r/p+D_r` is a unit modulo p, then `C_r(p)` has valuation exactly one. The B value has valuation one, while the A value is zero modulo p^2, so the gcd valuation is exactly one.
- If both residue tests vanish, then both values are zero modulo p^2. The Wronskian upper bound then forces the gcd valuation to be exactly two.

This establishes the three-way rule without requiring the individual A or B valuation beyond modulus p^2. The restriction `0<r<p` is essential for the Wronskian bound used here. The index r=0 is separately covered by `p` not dividing `g_p`.

## 4. Scope of the arithmetic advance

The local rule permits `v_p(g_n)>v_p(n!)`, as the stated n=107, p=73 example illustrates. Thus an assumed universal divisibility `g_n|n!` would be false. The proof does not bound the combined contribution from exceptional residues at different primes, and does not establish unboundedness or boundedness of `g_n*sqrt(n)/n!`. Higher lifts and a simultaneous-prime argument remain genuinely separate work.

## 5. Supplement: universal odd-surplus cap

The later Section 5 in `../gcd_arithmetic/ODD_PRIME_DIGIT_CRITERION.md` was also checked. Let `n>=2`, let p be odd, and set

    s=min(v_p(a_n),v_p(b_n))=v_p(g_n)-v_p(n!),
    e=floor(log_p(n-1)).

Here e is zero when `p>n-1`. Division of the exact integer Wronskian by `n!*(n-1)!` gives

    a_n*b_(n-1)-a_(n-1)*b_n=(-1)^n/n.

The first product on the left has valuation at least `s-e`, because `v_p(b_(n-1))>=-e`; the second has valuation at least s, because `a_(n-1)` is p-integral. Their difference therefore has valuation at least `s-e`, regardless of cancellation. The exact right side has valuation `-v_p(n)`, proving

    s<=e-v_p(n).

Thus the positive surplus exponent `max(s,0)` is at most `max(e-v_p(n),0)`. Since the latter is exactly the p-valuation of `L/gcd(L,n)`, for `L=lcm(1,...,n-1)`, the claimed divisibility of the odd positive-surplus product E_n follows prime by prime.

The stated corollaries are valid: each surplus prime is less than n; primes larger than `sqrt(n-1)` have surplus exponent at most one; and at `n=p^k` the bound is at most -1 for that p, so positive p-surplus is excluded. This universal finite cap is genuine additional information. Its total logarithm is only bounded on the linear scale O(n), and it does not establish an o(n) bound or resolve the global unboundedness question.
