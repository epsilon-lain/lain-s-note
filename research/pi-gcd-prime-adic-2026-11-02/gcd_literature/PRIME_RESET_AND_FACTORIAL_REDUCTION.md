# Prime-reset theorem and literature for the fixed arctangent continuants

2026-11-02. Context: epsilon-lain/lain-s-note at dff048f73da83e63978e1c01576d38ea71af1546, research/pi-joint-arithmetic-2026-11-02. This note adds a completely proved first-level odd-prime criterion and a partial 2-adic theorem. It does not decide whether g_n sqrt(n)/n! is unbounded.

Define A0=A1=1, B0=0, B1=1, and for either sequence X,

    X_(n+1)=(2n+1)X_n+n^2 X_(n−1), n>=1.

Put g_n=gcd(A_n,4B_n).

## 1. The exact first-level odd-prime theorem

For an odd prime p set chi_p=(-1)^((p+1)/2). Then:

(a) p divides A_n for every n>=p.
(b) A_(p−1) ≡ B_p ≡ chi_p (mod p).
(c) B_(p+r) ≡ chi_p A_r (mod p), 0<=r<=p.
(d) Consequently,

    p ∤ g_n                       if n<p;
    p | g_n iff p | A_(n−p)       if p<=n<2p;
    p | g_n                       if n>=2p.                          (1)

This identifies every first-level exceptional prime in n/2<p<=n via a small-index Legendre residue. It gives no upper bound on higher valuations.

Proof. Let L_n denote the standard Legendre polynomial. The recurrence gives

    A_n=n! i^(−n)L_n(i).

The coefficients of L_n have only powers of 2 in their denominators, and i^(−n)L_n(i) is rational (indeed in Z[1/2]). Thus A_n/n! is p-integral for every odd p, proving (a).

The identity

    L_n(x)=sum_(k=0)^n binom(n,k)binom(n+k,k)((x−1)/2)^k

shows L_(p−1)(x)≡1 mod p: all k>=1 coefficients contain p, while denominators 2^k are units. Wilson's theorem now yields

    A_(p−1) ≡ −i^(−(p−1)) = chi_p (mod p).

At n=p, the known Wronskian gives

    A_p B_(p−1)−A_(p−1)B_p = −((p−1)!)^2 ≡−1 (mod p),

hence B_p≡chi_p. At the reset index n=p the recurrence gives B_(p+1)≡B_p, and thereafter its coefficients repeat modulo p. The block B_(p+r) therefore has the same recurrence and initial data as chi_p A_r, proving (c), including r=p. Hence B_(2p)≡B_(2p+1)≡0, and both sequences remain divisible by p for all later indices. The n<p clause follows because a common factor p would divide the nonzero-unit Wronskian ((n−1)!)^2; n=0 can be handled directly. This proves (1).

For example A_7=463680 is divisible by 23, so p=23 enters g_30, although 23>30/2. This is a real exceptional prime, not an error in factorial normalization. In contrast p=3,5,7,11,13 have no residues r in 0<=r<p with A_r≡0 mod p; their first entry into g_n occurs at n=2p.

## 2. Provable 2-adic information

For every k>=0 one has

    v2(A_(4k))=2k,
    v2(A_(4k+1))=2k,
    v2(A_(4k+2))=2k+2,
    v2(A_(4k+3))>=2k+3.                                           (2)

Here v2(A0)=0 by convention for A0=1. Proof: inductively put a=A_(4k)/2^(2k), b=A_(4k+1)/2^(2k). Both are odd, and b≡a mod 8. This is true at k=0. The next recurrence gives

    A_(4k+2)/2^(2k) = (8k+3)b+(4k+1)^2 a ≡4a (mod 8),

so its valuation is exactly 2. The two terms producing A_(4k+3) each have valuation 2k+2 and odd normalized coefficients, so their sum has valuation at least 2k+3. The next step producing A_(4k+4) adds a term of valuation >=2k+3 to a term of exact valuation 2k+2, proving the first assertion for k+1. Finally

    A_(4k+5)−A_(4k+4)
      =8(k+1)A_(4k+4)+16(k+1)^2 A_(4k+3),

which proves equality of the normalized odd residues modulo 8 for k+1. This closes the induction.

Since v2(g_n)<=v2(A_n), (2) already exhibits an exponential 2-adic loss relative to n! on three residue classes. It cannot by itself rule out the target because odd-prime excess factors may compensate.

The companion `../gcd_arithmetic/TWO_ADIC_VALUATION.md` now proves the stronger full formula

    v2(A_n)=floor(n/2)+[n≡2 mod4]
                             +[n≡3 mod4]v2(n+1).

Its independent audit is `../gcd_audit2/INDEPENDENT_TWO_ADIC_AUDIT.md`. The elementary argument above provides an independent check of three residue classes. The separate observation v2(g_n)=v2(A_n) except n=3 remains conjectural here and is not used in any global bound.

## 3. An attractive false shortcut, with exact counterexamples

The odd part of g_n need not divide n!. Exact computation through n=2000 finds the first excess cases

    n=72:  v17(g_n)=v17(n!)+1;
    n=107: v73(g_n)=v73(n!)+1;
    n=166: v131(g_n)=v131(n!)+1;
    n=170: v41(g_n)=v41(n!)+1.

Thus a proof that combines the 2-adic loss with the unproved blanket inequality g_(n,odd)|n! would be invalid. Any global estimate must bound the collective contribution of these higher-valuation exceptions, not discard them.

## 4. Closest primary literature and the precise level of applicability

### Ben David et al. (2024), On the Connection Between Irrationality Measures and Polynomial Continued Fractions

https://link.springer.com/article/10.1007/s40598-024-00250-z
Published 18 June 2024, Arnold Mathematical Journal 10, 529–566.

This is the most directly useful research reference. It analyzes balanced-degree polynomial CFs, defines factorial reduction through the exponential growth of gcd_n/(n!)^deg a, and separates proved criteria from conjectured exact gcd forms. Appendix D Theorem 5 proves a factorial/lcm divisibility bound under reflection symmetries a_n=−a_(−1−n), b_n=b_(−n), b0=0. Our coefficients a_n=2n+1,b_n=n^2 satisfy these conditions.

Indexing matters: in the paper's standard CF normalization, its convergents are p_j=A_(j+1), q_j=B_(j+1). Therefore the quoted theorem applies with j=n−1. It implies

    (n−1)!/lcm(1,...,n−1) divides 2^(n−1) gcd(A_n,B_n),

hence also 2^(n−1)g_n. The Legendre explicit formulas may give somewhat stronger denominators, but none of these factorial-reduction lower bounds establishes the target.

The paper's FR criterion permits g_n/n!=exp(−c n+o(n)), c>0. Such a sequence has FR while g_n sqrt(n)/n! tends to zero. The current critical target is thus much stronger than mere FR. If one proves a negative exponential upper rate for this exact sequence, the route fails; if the rate is exactly zero, the polynomial factor sqrt(n) remains decisive. Section 2.5's more precise gcd structures include a still-unformulated 'Zebra' pattern and are explicitly conjectural.

The first-level reset p|g_n for n>=2p is already part of the general proof behind their Theorem 5. Formula (1)'s exact intermediate-block criterion specializes the reset further for the present initial data. No literature-priority claim is made for that specialization.

### Chao Wang (2026), A Family of Continued Fraction Identities for Arctangent Values

https://arxiv.org/abs/2601.11892
Submitted 17 January 2026; v3 28 March 2026.

This is exactly the Gauss fraction arctan(p/q)=p/(q+p^2/(3q+(2p)^2/(5q+...))). At p=q=1 its continuants are precisely the present recurrence (up to index and inversion conventions). It proves geometric convergence, but no factorial gcd asymptotic sufficient for the current critical target is claimed in its abstract. It should be cited to identify the classical CF mechanism, not as a solution of the arithmetic bottleneck.

### Alekseyev–Amdeberhan–Shallit–Vukusic (2026), On the p-adic valuations of values of Legendre polynomials

https://link.springer.com/article/10.1007/s40993-026-00743-0
https://www.math.tulane.edu/~tamdeberhan/Legendre_polynomials.pdf
Published 11 May 2026; author manuscript dated 4 July 2025.

This recent paper proves valuation formulae and p-regularity for P_n(p), and certain rational arguments r with v_p(r)>=1. It is a relevant source of proof methods (digit recurrences, valuation isolation, p-regular sequences), but its hypotheses DO NOT cover P_n(i): i is a p-adic unit, and the argument is algebraic rather than one of their rational inputs. It cannot directly certify the observed formula in Section 2 or the gcd's higher-prime behavior.

### Tanguy Rivoal (2007), Convergents and irrationality measures of logarithms

https://doi.org/10.4171/RMI/519
https://ems.press/content/serial-article-files/38218

The paper studies restricted denominator shapes via detailed arithmetic of the upper logarithmic Padé table. This is useful if the project moves from a single diagonal family to off-diagonal/multiple logarithm Padé with a controlled common denominator. Its stated real positive rational-logarithm results do not automatically transfer to pi as a complex logarithm at a root of unity.

## 5. Concrete next arithmetic work

1. Lift (1) from modulo p to modulo p^e while retaining the complete two-component state, rather than approximating all valuations by v_p(n!). The first excess n=72,p=17 is a mandatory regression case.
2. Use the companion's completed full 2-adic theorem to account for the exponential loss at every index. Any further exact formula for v2(g_n), rather than v2(A_n), would refine its polynomial factor but is not needed for the proved upper budget.
3. Bound the SUM of exceptional odd-prime valuation gains. Controlling each fixed p separately is insufficient when p ranges with n. If that sum is o(n), the proved 2-adic deficit would force g_n sqrt(n)/n! ->0 and close this particular route negatively.
4. Any use of the general FR literature must retain the zero exponential-rate boundary: a theorem that proves only FR does not advance the current success condition.
