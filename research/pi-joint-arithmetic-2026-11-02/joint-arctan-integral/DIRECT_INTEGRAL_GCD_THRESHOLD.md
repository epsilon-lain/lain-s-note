# A direct positive π integral and its exact reduced-denominator threshold

2026-11-02. This note identifies a natural direct rational-integral route with the Legendre–arctangent Padé sequence in the companion arithmetic study. The identification prevents counting a change of integral representation as a new approximating family. It also gives a sharp joint error/gcd criterion. No assertion is made that the criterion holds along an infinite subsequence.

## 1. A positive integral gives the same recurrence

For each integer n≥0 define

\[
I_n=\int_0^1\frac{x^n(1-x)^n}{(1+x^2)^{n+1}}\,dx.
\]

In particular I_n>0. With x=tan t, 0≤t≤π/4, put

\[
h(t)=\frac{\sin(2t)+\cos(2t)-1}{2}.
\]

Then I_n=∫_0^(π/4) h(t)^n dt, h(0)=h(π/4)=0, and

\[
(h')^2=1-4h-4h^2,\qquad h''=-4h-2.
\]

For n≥1 the boundary values of h^n h' vanish. Integrating its derivative gives

\[
0=nI_{n-1}-(4n+2)I_n-(4n+4)I_{n+1}.
\tag{1}
\]

The initial values are I_0=π/4 and I_1=1/2−π/8. Therefore J_n=(−2)^n I_n obeys

\[
(n+1)J_{n+1}=(2n+1)J_n+nJ_{n-1},
\quad J_0=\pi/4,\quad J_1=\pi/4-1.
\tag{2}
\]

Define rational sequences ℓ_n,w_n by this same recurrence and by

\[
\ell_0=\ell_1=1,\qquad w_0=0,\quad w_1=1.
\]

Equality of the initial values and the recurrence proves the all-n identity

\[
\boxed{(-2)^n I_n=\ell_n\frac\pi4-w_n.}
\tag{3}
\]

The recurrence makes ℓ_n>0. Thus r_n=4w_n/ℓ_n approximates π with exact error

\[
\boxed{\pi-r_n=\frac{4(-2)^n I_n}{\ell_n},
\qquad |\pi-r_n|=\frac{4\,2^n I_n}{\ell_n}.}
\tag{4}
\]

This is exactly the a=1 Legendre–arctangent family of the companion note. In particular, it is not a second independent family with separate approximation opportunities.

## 2. The analytic part of the joint cost has an explicit constant

Set λ=1+√2 and ρ=(√2−1)/2, so 2ρλ=1. The kernel h has a unique maximum at t=π/8, with

\[
h(\pi/8)=\rho,\qquad h''(\pi/8)=-2\sqrt2.
\]

An elementary Laplace estimate therefore gives

\[
I_n=\rho^n\sqrt{\frac{\pi\rho}{n\sqrt2}}
\bigl(1+O(n^{-1})\bigr).
\tag{5}
\]

For completeness, ℓ_n has a real integral representation

\[
\ell_n=\frac1\pi\int_0^\pi(1+\sqrt2\cos\theta)^n\,d\theta.
\tag{6}
\]

This does not require importing a Legendre asymptotic. The right side has initial values 1,1. Integrating the derivative of sin θ(1+√2 cos θ)^n shows that it obeys (2); hence (6) follows. Its unique largest absolute base is λ at θ=0. The negative-base portion has absolute value at most √2−1<λ and is exponentially negligible. Since

\[
1+\sqrt2\cos\theta=\lambda-\frac{\sqrt2}{2}\theta^2+O(\theta^4),
\]

the boundary version of the same estimate gives

\[
\ell_n=\lambda^n\sqrt{\frac{\lambda}{2\pi n\sqrt2}}
\bigl(1+O(n^{-1})\bigr).
\tag{7}
\]

One direct justification of the error terms in (5) and (7) is to split off a fixed neighborhood of the unique maximum. The complement is exponentially smaller. On a shrinking window of width n^(−2/5), the logarithm of the base has its quadratic term plus O(distance^4); after Gaussian rescaling, the integrated remainder is O(1/n). The tails outside that window are bounded by exp(−c n^(1/5)). Both expansions are even about the relevant maximum, so no cubic term occurs.

Equations (4)–(7) imply the precise analytic error scale

\[
\boxed{|\pi-r_n|=4\pi\lambda^{-2n-1}
\bigl(1+O(n^{-1})\bigr),}
\tag{8}
\]

and, crucially for the denominator comparison,

\[
\boxed{\ell_n^2|\pi-r_n|
=\frac{\sqrt2}{n}\bigl(1+O(n^{-1})\bigr).}
\tag{9}
\]

The square in (9) involves the rational coefficient ℓ_n, not the reduced integer denominator of r_n.

## 3. The entire remaining cost is an explicit integer gcd

Define integers a_n=n!ℓ_n and b_n=n!w_n. Their initial values and recurrence are

\[
a_0=a_1=1,\qquad b_0=0,\ b_1=1,
\]
\[
\boxed{x_{n+1}=(2n+1)x_n+n^2x_{n-1}\quad(n\ge1),}
\tag{10}
\]

for either sequence x=a or x=b. This proves integrality without any separate coefficient clearing argument. Let

\[
g_n=\gcd(a_n,4b_n),\qquad
p_n=\frac{4b_n}{g_n},\quad q_n=\frac{a_n}{g_n}.
\tag{11}
\]

Then q_n>0, gcd(p_n,q_n)=1, and p_n/q_n=r_n. The ratio of the genuine denominator to the rational coefficient is exactly

\[
\boxed{\frac{q_n}{\ell_n}=\frac{n!}{g_n}.}
\tag{12}
\]

Writing Q_n=q_n²|π−p_n/q_n|, combination with (9) now proves

\[
\boxed{Q_n=\frac{\sqrt2}{n}
\left(\frac{n!}{g_n}\right)^2
\bigl(1+O(n^{-1})\bigr).}
\tag{13}
\]

This is a joint error and *reduced* denominator statement. For any sequence of indices n_j→∞, it gives the exact equivalence

\[
\boxed{Q_{n_j}\longrightarrow0
\quad\Longleftrightarrow\quad
\frac{n_j!}{g_{n_j}\sqrt{n_j}}\longrightarrow0.}
\tag{14}
\]

Thus the arithmetic question for this fixed family is whether g_n can dominate n!/√n along an infinite subsequence. No such infinite subsequence is proved here. The full consecutive sequence cannot satisfy (14), by the separate adjacent-quality obstruction in the Legendre companion. This does not exclude sparse subsequences.

The reduced denominators q_n tend to infinity: (8) makes r_n→π, and a bounded-denominator sequence of rational numbers in a bounded interval has only finitely many possible values. Irrationality of π excludes convergence to π through such a finite set. This argument does not assert that q_n is monotone.

## 4. Adjacent gcd gains have an exact divisibility constraint

The recurrence (10), with the initial values, proves the integer Wronskian

\[
a_n b_{n-1}-a_{n-1}b_n
=(-1)^n[(n-1)!]^2\qquad(n\ge1).
\tag{15}
\]

Indeed the determinant at index n+1 is −n² times the determinant at index n; the index-1 determinant is −1. Because g_n divides both a_n and 4b_n, the two products in four times (15) are each divisible by g_ng_(n−1). Hence

\[
\boxed{g_n g_{n-1}\mid4[(n-1)!]^2.}
\tag{16}
\]

In particular, define h_n=g_n√n/n! for n≥1. Then for n≥2,

\[
\boxed{h_nh_{n-1}\le4\sqrt{\frac{n-1}{n}}.}
\tag{17}
\]

Equation (13) becomes Q_n=(√2/h_n²)(1+O(1/n)); therefore

\[
\boxed{\liminf_{n\to\infty}Q_nQ_{n-1}\ge\frac18.}
\tag{18}
\]

This is an asymptotic statement; it is not a uniform bound Q_nQ_(n−1)≥1/8 at every small n. It strengthens the limiting consequence of the companion's finite-n adjacent-quality bound. Large gcd gains sufficient for Q_n→0 cannot occur at every index or at two consecutive indices along a diverging sequence. An isolated sparse subsequence remains logically possible.

## 5. Every fixed index gap has a positive limiting quality obstruction

The sharp signed error (4), (8) supplies a broader statement with no gcd estimates. Fix an integer k≥1. For n≥k put e_n=π−r_n. Then

\[
\frac{e_n}{e_{n-k}}=(-1)^k\lambda^{-2k}(1+o(1)).
\]

In particular r_n and r_(n−k) are distinct for every sufficiently large n. Their reduced rational representations give the nonzero integer

\[
p_nq_{n-k}-p_{n-k}q_n,
\]

so q_nq_(n−k)|r_n−r_(n−k)|≥1. Consequently

\[
Q_nQ_{n-k}\ge
\frac{|e_ne_{n-k}|}{|e_{n-k}-e_n|^2}.
\]

Taking the limit inferior proves, for each fixed k,

\[
\boxed{\liminf_{n\to\infty}Q_nQ_{n-k}
\ge\frac{\lambda^{-2k}}
{[1-(-1)^k\lambda^{-2k}]^2}>0.}
\tag{19}
\]

For k=1 the right side is 1/8, agreeing with (18); for k=2 it is 1/32. In particular, neither the complete even-index sequence nor the complete odd-index sequence can have Q_n→0.

If an increasing subsequence n_j has Q_(n_j)→0, its successive gaps must satisfy n_(j+1)−n_j→∞. Otherwise an infinite subsequence of its adjacent index gaps would equal one fixed k, contradicting (19). This gives a necessary sparsity condition, not an exclusion of all possible infinite good subsequences.

## Verification and limits

`verify_direct_integral.py` independently integrates the rational function for small n by an exact Hermite-reduction coefficient solve; it does not assume recurrence (1) when doing that check. It then checks the integer recurrence, reduction, exact error signs, and joint-cost identities through n=100 using a rigorous rational Machin interval for π. The numerical values are finite certificates, not evidence that an infinite gcd condition holds.

The associated error/gcd criterion isolates a concrete target. It neither proves nor disproves π's bad approximability, and it does not establish a new rational approximation record.
