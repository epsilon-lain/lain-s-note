# From the actual reduced denominator to a fixed error exponent

2026-11-02. This note studies the same fixed Legendre–arctangent sequence. It gives an exact test and a sharp asymptotic threshold for an error bound C/q^c, with arbitrary fixed real constants C>0 and c>0. The constant C may in particular be Euler's number e or any prescribed positive error budget. It is independent of n.

The criteria below concern this constructed sequence. They do not establish an upper bound on the approximation quality of all rational numbers, nor a global irrationality exponent for π.

## 1. Objects and the two established asymptotics

For integers n≥0 define a_n,b_n by a_0=a_1=1, b_0=0,b_1=1 and

\[
x_{n+1}=(2n+1)x_n+n^2x_{n-1}\quad(n\ge1).
\]

Put

\[
g_n=\gcd(a_n,4b_n),\quad
p_n=4b_n/g_n,\quad q_n=a_n/g_n,
\]
\[
\ell_n=a_n/n!,\quad R_n=g_n/n!,
\quad \delta_n=|\pi-p_n/q_n|>0.
\]

The denominator q_n is positive and reduced; q_n=ℓ_n/R_n exactly. The direct positive-integral proof establishes, with λ=1+√2,

\[
\ell_n=L\,n^{-1/2}\lambda^n(1+O(n^{-1})),
\qquad L=\sqrt{\frac{\lambda}{2\pi\sqrt2}},
\tag{1}
\]
\[
\delta_n=D\lambda^{-2n}(1+O(n^{-1})),
\qquad D=\frac{4\pi}{\lambda}.
\tag{2}
\]

Also q_n→∞, although monotonicity is neither assumed nor needed. These statements are proved for the actual sequence, not conditional growth models for the gcd.

## 2. Exact test and the sharp asymptotic threshold

For fixed C,c>0 the following equivalence is exact at every index:

\[
\boxed{\delta_n\le Cq_n^{-c}
\quad\Longleftrightarrow\quad
R_n\ge T_n(C,c),
\qquad
T_n(C,c)=\ell_n(\delta_n/C)^{1/c}.}
\tag{3}
\]

All factors are positive real numbers, so raising to the positive power c preserves the inequality. Equations (1)–(2) compute the threshold rather than estimating the numerator and denominator separately:

\[
\boxed{T_n(C,c)=K(C,c)\,n^{-1/2}
\lambda^{(1-2/c)n}(1+O(n^{-1})),}
\tag{4}
\]
\[
\boxed{K(C,c)=
\frac{\lambda^{1/2-1/c}(4\pi/C)^{1/c}}
{\sqrt{2\pi\sqrt2}}.}
\tag{5}
\]

The O-term is for fixed C,c. It is not claimed uniform as c approaches zero or infinity.

A form that avoids any assumption that a growth rate exists is useful. Define α=ln λ and

\[
S_n(C,c)=\ln R_n-(1-2/c)\alpha n
+\tfrac12\ln n-\ln K(C,c).
\]

Then the actual normalized error satisfies

\[
\boxed{\ln\frac{\delta_n q_n^c}{C}
=-cS_n(C,c)+O(n^{-1}).}
\tag{6}
\]

Thus S_n→+∞ gives δ_nq_n^c→0, whereas S_n→−∞ gives δ_nq_n^c→∞. If S_n→s is finite, then δ_nq_n^c/C→exp(−cs). At s=0 this determines a limit equal to one; it does not decide the eventual direction of the non-strict inequality in (3). The exact threshold (3) retains that boundary information.

## 3. Exponential and polynomial parts of the criterion

The exponential threshold is

\[
\tau_c=(1-2/c)\alpha.
\]

If liminf (ln R_n)/n>τ_c, the bound δ_n≤Cq_n^(−c) holds eventually with an exponentially vanishing normalized error. If limsup (ln R_n)/n<τ_c, it fails eventually with an exponentially growing normalized error. Equality of these exponential rates leaves the polynomial and constant factors decisive.

For example, suppose that along any specified subsequence of indices,

\[
R_n=K_0 n^{-a}\exp(\tau n)(1+o(1)),
\quad K_0>0,\ a,\tau\in\mathbb R.
\]

This is a conditional description of the actual gcd ratio, not a claim that arbitrary such profiles are realizable. On that subsequence,

\[
\boxed{\delta_nq_n^c
=D(L/K_0)^c\,n^{c(a-1/2)}
\exp([c(\alpha-\tau)-2\alpha]n)(1+o(1)).}
\tag{7}
\]

At the critical rate τ=τ_c, a<1/2 gives a vanishing normalized error, a>1/2 gives divergence, and a=1/2 gives the finite limit D(L/K_0)^c. Comparing that constant with C settles strict-margin cases; exact equality again needs finer terms.

## 4. The exponent two retains a separate constant problem

At c=2,

\[
\boxed{K(C,2)=2^{1/4}/\sqrt C,
\qquad T_n(C,2)=\frac{2^{1/4}}{\sqrt{Cn}}
(1+O(n^{-1})).}
\tag{8}
\]

In particular, with h_n=R_n√n,

\[
\boxed{q_n^2\delta_n=\frac{\sqrt2}{h_n^2}
(1+O(n^{-1})).}
\tag{9}
\]

Therefore h_n→∞ is exactly the condition for q_n²δ_n→0 along a diverging sequence of indices. A positive finite limit h gives a positive limiting quality √2/h²; h_n→0 gives divergent quality. Merely proving an exponent tending to two does not distinguish these outcomes.

For instance, if R_n=K_0n^(−a)(1+o(1)), all such profiles have constructed exponents tending to two, but their quadratic qualities are

\[
q_n^2\delta_n=(\sqrt2/K_0^2)n^{2a-1}(1+o(1)).
\]

They tend respectively to zero, a positive constant, or infinity as a is below, equal to, or above 1/2.

## 5. Limsup and liminf exponents without a rate-limit assumption

For the eventually valid q_n>1 define the observed exponent

\[
c_n=\frac{-\ln\delta_n}{\ln q_n}.
\]

The two established expansions give

\[
\boxed{c_n=
\frac{2\alpha n+\alpha-\ln(4\pi)+O(n^{-1})}
{\alpha n-\ln R_n-\frac12\ln n
+\frac12\alpha-\frac12\ln(2\pi\sqrt2)+O(n^{-1})}.}
\tag{10}
\]

Let τ⁺=limsup (ln R_n)/n and τ⁻=liminf (ln R_n)/n, allowing extended values. Since q_n→∞, τ⁺≤α. Then

\[
\boxed{\limsup c_n=\frac{2\alpha}{\alpha-\tau^+},
\qquad
\liminf c_n=\frac{2\alpha}{\alpha-\tau^-}.}
\tag{11}
\]

Use 2α/0=+∞ and 2α/(+∞)=0. The proof is simply that (ln q_n)/n=α−(ln R_n)/n+o(1) is positive, while (−ln δ_n)/n→2α. No existence of a limit for (ln R_n)/n is assumed. The same formulas apply after restricting to any diverging subsequence.

These are exponents of this constructed sequence. Their value below two would only say that this sequence misses the critical scale; it would not bound the exponent of π over all rational approximants.

## 6. Positive odd excess and odd deficit enter with opposite signs

Write the reduced positive rational number R_n uniquely as

\[
R_n=2^{d_n}\frac{E_n}{D_n},
\quad d_n\in\mathbb Z,
\quad E_n,D_n\in\mathbb Z_{>0}\text{ odd},
\quad\gcd(E_n,D_n)=1.
\]

The odd numerator E_n is positive excess beyond n!, and D_n is the remaining odd deficit. They must both be retained:

\[
\ln R_n=d_n\ln2+\ln E_n-\ln D_n.
\tag{12}
\]

The exact criterion (3) is therefore equivalent to

\[
\boxed{\ln(E_n/D_n)
\ge\ln T_n(C,c)-d_n\ln2.}
\tag{13}
\]

The right side equals

\[
\tau_c n-\tfrac12\ln n+\ln K(C,c)-d_n\ln2+O(n^{-1}).
\]

There is now an exact 2-adic gcd statement, proved in the [saturation note](../exponent-audit/TWO_ADIC_GCD_SATURATION.md). Here v_2 is the exponent of 2 in a nonzero rational number, and s_2(n) counts the ones in the binary expansion of n. The preceding valuation theorem gives f(n)=v_2(a_n) explicitly:

\[
f(4m)=f(4m+1)=2m,\quad f(4m+2)=2m+2,\quad
f(4m+3)=2m+3+v_2(m+1).
\]

Then

\[
d_n=f(n)-[n-s_2(n)]-\mathbf1_{\{n=3\}}.
\tag{14}
\]

The single exceptional index n=3 has no asymptotic effect. In particular, define the explicit correction β_n=d_n+n/2=O(ln n), and put

\[
\sigma_c=\tfrac12\ln2+(1-2/c)\alpha.
\]

The threshold in (13) is now an exact net-prime-budget test with asymptotic value

\[
\boxed{\frac{E_n}{D_n}\ge2^{-d_n}T_n(C,c),
\qquad
2^{-d_n}T_n(C,c)=
2^{-\beta_n}K(C,c)n^{-1/2}
\exp(\sigma_c n)(1+O(n^{-1})).}
\tag{15}
\]

Thus the power-of-two correction is retained explicitly rather than hidden inside an unspecified polynomial loss. In particular success at exponent c requires

\[
\boxed{\ln(E_n/D_n)
\ge[\tfrac12\ln2+(1-2/c)\alpha]n-O_{C,c}(\ln n).}
\tag{16}
\]

Equation (16) discards the explicit logarithmic correction and is only a necessary exponential-rate condition. The exact test is (13), or equivalently (15). A large E_n alone does not establish it if D_n is also large.

Writing σ⁺=limsup ln(E_n/D_n)/n and σ⁻=liminf ln(E_n/D_n)/n, one may also substitute τ⁺=−ln2/2+σ⁺ and τ⁻=−ln2/2+σ⁻ in (11). This yields

\[
\boxed{\limsup c_n=
\frac{2\alpha}{\alpha+\frac12\ln2-\sigma^+},
\qquad
\liminf c_n=
\frac{2\alpha}{\alpha+\frac12\ln2-\sigma^-},}
\tag{17}
\]

with the same extended-value conventions. Since q_n→∞, the first denominator is nonnegative.

For example, if limsup ln(E_n/D_n)/n≤0, the 2-adic loss implies the conditional upper bound

\[
\boxed{\limsup c_n\le
\frac{2\ln(1+\sqrt2)}{\ln(1+\sqrt2)+\frac12\ln2}<2.}
\tag{18}
\]

Proving the required net odd-prime upper bound is a separate arithmetic problem. Neither a finite list of positive excess primes nor a bound on each isolated prime establishes it.

## 7. An unconditional sufficient positive exponent for this sequence

There is also a weak positive exponent that this family actually achieves; this is a sufficient bound, not an optimality or novelty claim. The exact convolution identity is

\[
w_n=\sum_{j=1}^n\frac{\ell_{j-1}\ell_{n-j}}j,
\qquad w_n=b_n/n!.
\tag{19}
\]

For a short proof, let F(z)=∑ℓ_nz^n=(1−2z−z²)^(−1/2). Then G(z)=F(z)∫_0^zF(t)dt has G(0)=0,G'(0)=1 and satisfies

\[
(1-2z-z^2)G'(z)-(1+z)G(z)=1.
\]

Its coefficient recurrence is exactly that of w_n; taking coefficients gives (19). Since every ℓ_j has only powers of two in its denominator, (19) shows that the remaining odd deficit D_n divides the odd part of

\[
\mathcal L_n=\operatorname{lcm}(1,2,\ldots,n).
\]

Indeed, for each odd prime p, v_p(w_n)≥−⌊log_p n⌋ while v_p(ℓ_n)≥0; therefore v_p(g_n/n!)≥−⌊log_p n⌋. The exact four-case formula for d_n also gives d_n≥−n/2 for every n≥1. Hence

\[
\boxed{R_n\ge\frac{2^{-n/2}}{\mathcal L_n}.}
\tag{20}
\]

Now invoke the standard prime number theorem in its equivalent Chebyshev-function form

\[
\ln\mathcal L_n=\psi(n)=n+o(n).
\tag{21}
\]

Together with (1), (20) gives ln q_n≤(α+ln2/2+1)n+o(n). Equations (2) and this bound prove

\[
\boxed{\liminf c_n\ge c_{\mathrm{PNT}},
\qquad
c_{\mathrm{PNT}}=
\frac{2\ln(1+\sqrt2)}
{\ln(1+\sqrt2)+\frac12\ln2+1}
=0.7911979206\ldots.}
\tag{22}
\]

More explicitly, for **every fixed** C>0 and 0<c<c_PNT,

\[
\boxed{\delta_nq_n^c\longrightarrow0,
\qquad \delta_n\le Cq_n^{-c}\text{ for all sufficiently large }n.}
\tag{23}
\]

This includes C=e and c=3/4. For an elementary strict comparison, λ>12/5, e<11/4 and √2<10/7 give λ^5>e³2^(3/2), which is equivalent to 3/4<c_PNT. It is an existence of an eventual threshold, not an effective minimal starting index. The equality case c=c_PNT is not decided by the bound. The exponent in (22) is much smaller than two and does not establish critical quadratic approximation.

An elementary version avoids (21): the divisibilities
\(\mathcal L_{2m}/\mathcal L_m\mid\binom{2m}m\) and
\(\mathcal L_{2m+1}/\mathcal L_m\mid(2m+1)\binom{2m}m\)
give ln \(\mathcal L_n\)≤n ln4+O((ln n)²) by repeated halving. The prime-power factors in each quotient occur at most once and are present in the indicated binomial coefficient or the extra factor 2m+1. Substituting this bound proves (23) for

\[
0<c<\frac{2\alpha}{\alpha+\frac12\ln2+\ln4}
=0.6742862693\ldots
\]

without the prime number theorem. Both sufficient bounds use the full reduced-denominator model, and neither asserts that the family cannot do better.

## Verification and dependencies

Published proof paths, the prior checkpoint, and the PNT reference are listed in [SOURCES.md](../SOURCES.md).

The exact finite verifier checks actual rational approximants and rigorous π intervals for several rational c,C. Rational exponents c=k/l are certified by comparing δ_n^l q_n^k with C^l, avoiding floating-point powers. Finite classifications do not prove the asymptotic thresholds.

The proof inputs are the direct positive-integral identities and asymptotics, the exact reduced denominator formula, and the all-n 2-adic numerator valuation theorem from the preceding checkpoint. No claim about π's bad approximability follows here.
