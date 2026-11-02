# Reflection folds the exceptional-prime problem, but does not bound its diagonal mass

2026-11-02. This note studies the fixed sequence

\[
A_0=A_1=1,\qquad A_{r+1}=(2r+1)A_r+r^2A_{r-1}.
\]

Write ell_r=A_r/r!, lambda=1+sqrt(2), and let p always denote an odd prime. The intended target is the uniform upper bound

\[
S(n)=\sum_{\substack{n/2<p<n\\p\mid A_{n-p}}}\ln p\le n/6
\]

for all sufficiently large integer n. This target is not proved below. We give an exact reflection reformulation and remove a sublinear-degree edge contribution. Neither is a pointwise distribution theorem for the remaining interior.

## 1. An exact reflected coefficient identity

The generating function is F(z)=(1−2z−z²)^(-1/2). In F_p[[z]], Frobenius and the common constant term 1 give

\[
F(z)=(1-2z-z^2)^{(p-1)/2}F(z^p).
\]

Thus, for 0<=r<p, ell_r is the coefficient of z^r in the polynomial

\[
H_p(z)=(1-2z-z^2)^{(p-1)/2}.
\]

Put h=(p−1)/2 and chi_p=(−1)^h. Reversing that polynomial gives

\[
z^{p-1}H_p(1/z)
=(z^2-2z-1)^h
=\chi_p H_p(-z).
\]

Coefficient comparison proves the complete residue identity

\[
\boxed{\ell_{p-1-r}\equiv\chi_p(-1)^r\ell_r\pmod p.}
\tag{1}
\]

In particular, its zero set is invariant under r↦p−1−r. This is the specialization of the classical Holt reflection to this scaled Legendre family, not a new general congruence.

Because r! is a p-unit for r<p, vanishing of ell_r and of A_r are equivalent. Therefore, for p in the open interval (n/2,n), put

\[
k(n,p)=\min\{n-p,\ 2p-n-1\}.
\]

Both arguments are integers between 0 and p−1. Equation (1) yields

\[
\boxed{p\mid A_{n-p}\iff p\mid A_{k(n,p)},
\qquad 0\le k(n,p)\le(n-1)/3.}
\tag{2}
\]

The last inequality follows by separating p<=(2n+1)/3 and p>=(2n+1)/3. In the first range use k<=2p−n−1; in the second use k<=n−p. The transition value can lie at a half-integer; no rounding ambiguity affects the displayed real bound.

For a fixed integer k, the only possible primes assigned to this folded degree are

\[
p=n-k\quad\hbox{or}\quad p=(n+1+k)/2.
\tag{3}
\]

The latter value is used only when it is an integer and the prime lies in the required interval. At the transition, the two values can coincide and must be counted only once.

## 2. A smaller integer for the divisibility condition

Define the integer generalized central trinomial coefficient

\[
U_k=2^k\ell_k=[x^k](1+2x+2x^2)^k
\]

and its positive odd part

\[
C_k=U_k/2^{v_2(U_k)}.
\]

For k<p, all omitted factorials and powers of two are p-units, so the anomalous condition in (2) is precisely p|C_k.

The previous checkpoint proved ell_k<=lambda^k and the exact binary valuation of A_k. In particular

\[
v_2(U_k)=k+v_2(A_k)-v_2(k!)\ge k/2.
\]

Hence a useful uniform height bound is

\[
\boxed{1\le C_k\le(\sqrt2\lambda)^k=(2+\sqrt2)^k.}
\tag{4}
\]

This uses the binary valuation theorem for A_k, not the additional saturation theorem for gcd(A_k,4B_k). It is also valid at k=0, where C_0=1.

An anomalous prime must consequently satisfy

\[
\boxed{k(n,p)\ge\frac{\ln p}{\ln(2+\sqrt2)}
>\frac{\ln(n/2)}{\ln(2+\sqrt2)}.}
\tag{5}
\]

Thus genuinely empty logarithmic-width strips occur near both ends p=n/2 and p=n. This is an exact exclusion with a proved constant, not a statistical assertion about zeros. It does not remove any fixed positive fraction of the interval.

For reference, by the Schur congruence ell_n≡ell_(n−p) modulo p in this interval, S(n) is also the logarithm of the squarefree product of the prime divisors of C_n lying in (n/2,n). The elementary height bound for C_n is much too large to imply n/6; the prime-number-theorem bound S(n)<=n/2+o(n) is stronger than that direct height estimate.

## 3. All folded degrees up to o(sqrt(n)) cost only o(n)

Fix an integer K>=0 and define

\[
S_{\mathrm{edge}}(n,K)=
\sum_{\substack{n/2<p<n\\p\mid A_{n-p}\\k(n,p)\le K}}\ln p.
\]

Assign each distinct prime in this sum to its unique folded degree k(n,p). For each k<=K, the product of the assigned primes divides C_k. Even if both distinct primes in (3) occur, their product still divides C_k. Multiplying over k gives the exact integer divisibility

\[
\prod_{\substack{n/2<p<n\\p\mid A_{n-p}\\k(n,p)\le K}}p
\ \Bigm|\ \prod_{k=0}^{K}C_k.
\tag{6}
\]

Using (4),

\[
\boxed{S_{\mathrm{edge}}(n,K)
\le\sum_{k=0}^K\ln C_k
\le\frac{\ln(2+\sqrt2)}2K(K+1).}
\tag{7}
\]

Therefore K=o(sqrt(n)) gives S_edge(n,K)=o(n), uniformly in n. One convenient choice is K=floor(n^(1/3)). It leaves only

\[
n^{1/3}<k(n,p)\le(n-1)/3
\]

as the possible source of a linear-size exceptional mass. This proves a genuine reduction in the global lemma: boundary degrees do not affect its asymptotic constant. It does not prove the necessary bound for this interior range.

## 4. Why the remaining diagonal problem is not a fixed-prime zero count

For each fixed p, reflection pairs zero degrees and the recurrence rules out adjacent zeros. But S(n) samples degree n−p at a different modulus for every summand. A bound on the number of zero degrees for one prime does not bound how many different primes hit that particular moving diagonal.

Nor does an elementary residue-class exclusion close this gap. There are anomalies in all four odd residue classes modulo 8. Concrete certified examples are p=17 with r=4, p=107 with r=11, p=37 with r=5, and p=23 with r=7. Prime p=41 has zeros of both degree parities, including r=6 and r=13. These disprove simple rules excluding one class modulo 8 or one degree parity uniformly. They do not disprove every possible more refined arithmetic restriction.

The exact reflected identity is useful for search and for separating edges from the interior. It supplies no independence across moduli. The estimate S(n)<=n/6 for all large n remains unproved.

## Sources, verification, and scope

Tony D. Noe, *On the Divisibility of Generalized Central Trinomial Coefficients*, Journal of Integer Sequences 9 (2006), Article 06.2.7, Sections 6–8, especially equations (12)–(14) and Theorem 8.6: https://cs.uwaterloo.ca/journals/JIS/VOL9/Noe/noe35.pdf. The publisher-hosted primary paper was read. It establishes the general Schur/Holt framework; the folded-degree and edge-budget consequences here are direct deductions for the specified fixed sequence. No uniform moving-diagonal estimate is attributed to that paper.

Proof inputs from previous checkpoints are ell_k<=lambda^k and the all-index binary valuation formula for A_k. The companion standard-library script checks (1), the folding equivalence, exact odd cores, the height bound in Q(sqrt(2)), divisibility (6), and the small counterexamples to class/parity shortcuts. Finite checks do not establish a uniform S(n) bound or any bad-approximability conclusion.
