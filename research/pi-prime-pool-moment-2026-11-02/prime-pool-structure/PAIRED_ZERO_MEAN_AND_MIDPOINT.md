# Paired zeros, the midpoint exception, and a sufficient first-moment target

2026-11-02. This note refines the active-prime capacity problem for the fixed sequence

\[
U_0=1,\quad U_1=2,\qquad
(n+1)U_{n+1}=2(2n+1)U_n+4nU_{n-1}.
\]

It does not prove the requested uniform active-prime capacity bound. The new deduction is that the exceptional unpaired midpoint zeros have logarithmic prime mass o(x), so a quantitative first-moment bound on all digit zeros would suffice without any midpoint main term.

## 1. Exact pairing identity

For each odd prime p, put

\[
z_p=\#\{0\le r<p:U_r\equiv0\pmod p\},\qquad
m_p=\mathbf1_{U_{(p-1)/2}\equiv0\pmod p}.
\]

The established Holt reflection is

\[
U_{p-1-r}\equiv(-4)^{(p-1)/2-r}U_r\pmod p.
\]

The scalar factor is nonzero modulo p, so the zero set is invariant under r -> p-1-r. Its unique possible fixed point is h=(p-1)/2. Thus

\[
z_p=2k_p+m_p,\qquad k_p\in\mathbb Z_{\ge0}.
\]

Since p belongs to the active pool Z exactly when z_p>0, this gives the exact pointwise inequality

\[
\boxed{\mathbf1_{p\in\mathcal Z}\le\frac{z_p+m_p}{2}.}
\tag{1}
\]

No probabilistic independence is used. The midpoint correction cannot simply be discarded at the pointwise level: for example p=31 has z_p=m_p=1.

## 2. The midpoint curve and a short proof that it has no CM

Consider

\[
E:\ y^2=x(x^2+2x+2).
\]

Its discriminant is -256 and its j-invariant is 128. For p>=5 it has good reduction. Write a_p=p+1-#E(F_p). The already established coefficient/Hasse argument gives

\[
a_p\equiv U_{(p-1)/2}\pmod p,
\qquad |a_p|\le2\sqrt p<p,
\]

and consequently

\[
m_p=1\quad\Longleftrightarrow\quad a_p=0
\qquad(p\ge5).
\tag{2}
\]

For completeness, the coefficient identity follows from Euler's criterion applied to f(x)=x(x^2+2x+2): the only positive power divisible by p-1 in f(x)^h is x^(p-1). Its coefficient is [x^h](x^2+2x+2)^h=U_h, the central coefficient being unchanged by reversal. Summing powers over F_p gives the stated trace congruence.

Here one can certify that E is non-CM without invoking a table of rational CM j-invariants. Direct point counts give

| p | #E(F_p) | a_p | a_p^2-4p | Ordinary Frobenius field |
|---:|---:|---:|---:|---|
| 3 | 2 | 2 | -8 | Q(sqrt(-2)) |
| 5 | 8 | -2 | -16 | Q(i) |

Both reductions are ordinary because p does not divide a_p. If E had CM by an imaginary quadratic field K, specialization at any prime of good ordinary reduction would embed K into the geometric endomorphism algebra of that reduction. The latter is the quadratic field Q(sqrt(a_p^2-4p)), so it must equal K. The two displayed, nonisomorphic fields rule this out. Possible changes of the endomorphism order at primes dividing its conductor do not affect this argument on the rational endomorphism algebra.

The specialization principle and the ordinary Frobenius-field test are standard CM theory; see Charles (2004), Section 5.1, especially the discussion following Theorem 5.8. The finite fields and point counts here are exact certificates for this particular curve.

## 3. The midpoint correction has logarithmic mass o(x)

Serre's non-CM trace-zero theorem (1981, Section 8.2, Theorem 20 and Corollary 2) implies

\[
\#\{p\le x:a_p(E)=0\}=o(x/\ln x).
\]

Using (2), and absorbing the finite small-prime exceptions, it follows that

\[
\boxed{\sum_{p\le x\atop p\text{ odd}}m_p\ln p=o(x).}
\tag{3}
\]

Indeed each logarithm is at most ln x. This application is unconditional; no Lang-Trotter asymptotic or GRH is required. It is a statement about this one fixed elliptic curve, not a density theorem for all degrees r.

## 4. A rigorous first-moment sufficient condition

Define the weighted zero first moment

\[
T(x)=\sum_{p\le x\atop p\text{ odd}}z_p\ln p.
\]

Combining (1) and (3) proves

\[
\boxed{\vartheta_{\mathcal Z}(x)\le\tfrac12T(x)+o(x).}
\tag{4}
\]

Thus a sufficient target is: for some eta>0 and all sufficiently large x,

\[
\boxed{T(x)\le\left(1+\frac{\ln2}{2}-\eta\right)x.}
\tag{5}
\]

The threshold coefficient is 1.346573590279... . If (5) holds, (4) gives a strict capacity gap below tau=(2+ln2)/4. More directly, the established multi-band bound and W(n)<=theta_Z(n) yield

\[
\ln Q_n\ge(2+\ln2)n-2T(n)+o(n)
\ge2\eta n+o(n)\longrightarrow+\infty.
\]

This would exclude every bounded-quality subsequence in the fixed family. The first-moment estimate (5) is not proved here. Its finite numerical plausibility does not replace a uniform estimate.

The factor one half in (4) is sharp using only reflection and the information m_p=0: an abstract reflected zero set with exactly one non-midpoint pair has z_p=2 and active indicator one. Any stronger coefficient must use extra arithmetic information about additional pairs or their occurrence across primes.

## 5. A representation-theoretic route inspected, and its limit

Anandavardhanan (2025), Theorem 1.1, relates Legendre values to p-adic reductions of torus correlation coefficients. Our U_r is a nonzero scalar multiple of P_r(-i), so for p=3 modulo 4 this gives another exact interpretation using alpha=-1 and u=1.

Theorem 1.2, however, concerns existence of a suitable parameter over F_(p^m), m>=2. It is proved by counting the roots of a product of the first half of the Legendre polynomials; the available degree bound is quadratic in p. It neither controls the fixed parameter -i as p varies nor gives (5). Remark 5 explicitly leaves the analogous prime-field existence problem unresolved. This is useful scope information, not a new obstruction theorem for the present pool.

## Sources

- Jean-Pierre Serre, Quelques applications du theoreme de densite de Chebotarev, Publications Mathematiques de l'IHES 54 (1981), 123-201, Section 8.2, Theorem 20 / Corollary 2, printed p.189. Primary scan: https://www.numdam.org/item/PMIHES_1981__54__123_0.pdf . Only the non-CM trace-zero density consequence is used.
- Denis Xavier Charles, Complex Multiplication Tests for Elliptic Curves, arXiv:math/0409501 (2004), Section 5.1 and Corollary 6.2. https://arxiv.org/pdf/math/0409501 . Used for the ordinary-reduction CM-field principle; the two small point counts are computed here.
- U.K. Anandavardhanan, Simultaneous nonvanishing of the correlation constant, arXiv:2507.14616 (2025), Theorems 1.1, 1.2, Sections 5-6 and Remark 5. https://arxiv.org/pdf/2507.14616 . Its parameter-existence conclusion is not substituted for the required fixed-parameter prime-density bound.

The established reflection and midpoint trace identification also appear in the prior repository checkpoint's EDGE_CAPACITY_AND_MIDPOINT.md. This note's additional step is the non-CM certification and unconditional removal of the midpoint main term from the first-moment criterion.
