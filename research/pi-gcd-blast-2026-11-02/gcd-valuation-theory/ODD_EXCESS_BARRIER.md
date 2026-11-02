# The gcd-growth hypothesis: a false shortcut and a necessary odd-excess barrier

2026-11-02. This note studies only the fixed Legendre–arctangent π approximation family. It attacks its reduced-denominator gain directly. A tempting divisibility shortcut is false. An all-n 2-adic valuation theorem nevertheless forces any successful gcd gain to contain exponentially large odd excess. No bound excluding that odd excess is proved here.

## 1. Fixed integer model and the false divisibility shortcut

For n≥0 define integers a_n,b_n by a_0=a_1=1, b_0=0, b_1=1 and

\[
x_{n+1}=(2n+1)x_n+n^2x_{n-1}\quad(n\ge1).
\]

Let

\[
g_n=\gcd(a_n,4b_n),\quad
p_n=4b_n/g_n,\quad q_n=a_n/g_n,
\quad Q_n=q_n^2|\pi-p_n/q_n|.
\]

The previously proved joint asymptotic is

\[
Q_n=\frac{\sqrt2}{n}\left(\frac{n!}{g_n}\right)^2
\bigl(1+O(n^{-1})\bigr).
\tag{1}
\]

The assertion that the odd part of g_n always divides the odd part of n! is **false**. Already at n=72,

\[
v_{17}(g_{72})=5,
\qquad v_{17}(72!)=4.
\tag{2}
\]

The numerator of the reduced rational number g_72/72! has odd part exactly 17. The companion exact certificate reconstructs the integer recurrence, gcd, factorial and prime valuations. The later examples n=107,p=73 and n=176,p=23 are additional finite counterexamples, not an inferred general pattern.

## 2. Exact powers of two in a_n

Write v_2(k) for the exponent of 2 dividing a nonzero integer k. The [companion valuation proof](two-adic-proof/TWO_ADIC_NUMERATOR_VALUATION.md) establishes, for every integer m≥0,

\[
\boxed{
\begin{aligned}
v_2(a_{4m})&=2m,\\
v_2(a_{4m+1})&=2m,\\
v_2(a_{4m+2})&=2m+2,\\
v_2(a_{4m+3})&=2m+3+v_2(m+1).
\end{aligned}}
\tag{3}
\]

This is an all-n theorem. Its odd-index case uses a convergent 2-adic inverse-binomial series and an exact formal duplication identity; it is not extrapolated from the recurrence table. It implies the deliberately loose uniform bound

\[
v_2(a_n)\le n/2+\log_2 n+1\qquad(n\ge1).
\tag{4}
\]

Since g_n divides a_n and v_2(n!)=n-s_2(n), where s_2(n) counts binary ones, we obtain

\[
v_2(g_n/n!)\le -n/2+2\log_2 n+3.
\tag{5}
\]

Here v_2 on rational numbers means numerator valuation minus denominator valuation. The estimate uses s_2(n)≤⌊log_2 n⌋+1; the extra constant in (5) is generous.

## 3. Odd excess must offset an exponential loss

Write the reduced positive rational number g_n/n! uniquely as

\[
\frac{g_n}{n!}=2^{d_n}\frac{u_n}{v_n},
\qquad u_n,v_n\in\mathbb Z_{>0}\text{ odd},\quad\gcd(u_n,v_n)=1.
\]

Thus u_n is the odd part of its reduced numerator. It measures *positive* odd-prime excess beyond n!, not the total odd part of g_n. Equations (5) and v_n≥1 give

\[
\boxed{\frac{g_n}{n!}\le8n^2\,2^{-n/2}u_n.}
\tag{6}
\]

The analytic multiplier in (1), after multiplying by n, tends to √2>1. Therefore for all sufficiently large n,

\[
\boxed{Q_n\ge\frac{2^n}{64n^5u_n^2}.}
\tag{7}
\]

Consequences for this fixed family are now precise:

* A polynomial bound u_n≤Cn^d for all sufficiently large n would imply Q_n→∞, and hence would completely exclude a Q→0 subsequence within this family.
* More generally, any bound u_n≤exp((ln 2/2−ε)n) with a fixed ε>0 eventually would still imply Q_n→∞.
* If Q_n≤M for some fixed positive M along arbitrarily large indices, then those indices must satisfy

\[
\boxed{u_n\ge\frac{2^{n/2}}{8\sqrt M\,n^{5/2}}.}
\tag{8}
\]

In particular, a Q→0 subsequence would require exponentially large odd excess at least at the critical rate (ln 2)/2, up to polynomial factors. The finite counterexample (2) shows that setting u_n=1 cannot prove this obstruction. It does not show that the exponentially large excess required by (8) occurs.

There is also an all-index version. The [Legendre error identity](../../pi-joint-arithmetic-2026-11-02/joint-log-pade/LEGENDRE_ARCTAN_JOINT_COST.md), together with ∫_(−1)^1 P_n(t)² dt=2/(2n+1) and 1/(1+t²)>1/2 almost everywhere, proves

\[
Q_n>\frac2{2n+1}\left(\frac{n!}{g_n}\right)^2.
\]

Consequently, without an unspecified asymptotic threshold, for every n≥1,

\[
\boxed{Q_n>\frac{2^n}{32(2n+1)n^4u_n^2}
\ge\frac{2^n}{96n^5u_n^2}.}
\tag{9}
\]

Thus Q_n≤M at any such index forces u_n>2^(n/2)/(√(96M)n^(5/2)). The asymptotic estimate (8) has the slightly better constant; (9) is uniform. The symbol P_n in the cited error identity denotes the classical Legendre polynomial, not the earlier auxiliary determinant polynomial.

This isolates the remaining arithmetic task: prove a subcritical odd-excess bound, or construct and certify a genuinely supercritical subsequence. Neither outcome is established by this note. Excluding this particular family would not prove π badly approximable, because that requires a bound for every reduced rational approximation.

## Verification

`verify_odd_excess_counterexamples.py` gives exact integer certificates for the counterexamples, including the first failure through n=72. The companion 2-adic verifier checks the finite algebraic identities used by its analytic proof and compares the valuation formula with a finite recurrence table. All-n statements above rely on the proofs, not on the extent of those tables.
