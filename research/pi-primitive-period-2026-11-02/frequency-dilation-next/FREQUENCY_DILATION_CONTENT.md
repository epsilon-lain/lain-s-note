# Integer frequency dilation: exact content savings and their limits

2026-11-02. This note isolates the effect of replacing raw frequencies 0,1 by 0,b, where b is a positive integer. It gives an exact scaling law, exact primitive coefficient content, and a separately sourced central-limit corollary showing that each fixed b>1 is eventually worse than b=1 after primitive normalization. It also gives a fixed-dimension obstruction to making primitive values small by sending b to infinity. It does not determine the best joint choice of b and dimension or prove a rational approximation result.

## 1. Keep the complete source and budget explicit

For integers n>=2 and b>=1, take two centers and one spatial variable. The 2n source rows have labels (j,s,0), j=0,1 and 0<=s<n. The 2n raw columns are Y^h X^a with h=0,b and 0<=a<n, ordered by frequency then a.

A legal complete-packet realization is

\[
H_b=\max\{2n,2n+b-2,5\},\quad w_0=1,\quad w_1=2,
\quad v_0=H_b/n,\quad v_1=H_b,\quad\theta=2/H_b.
\]

The strict source inequality v_0s+v_1 beta<H_b gives exactly beta=0 and s<n. Every column satisfies h+2a<=H_b. The choice also gives v_0>w_0 and, with K=2, both K theta<1 and K(w_0/v_0)theta<1. Thus increasing b does have a height cost; it is not a free change inside the original height 2n when b>2. For n>=3 and b<=2, the height remains exactly 2n; the harmless lower bound 5 covers n=2 without violating the strict theta condition.

Define

\[
P_{n,b}(Z)=\det\left[[z^s]e^{hz}(jZ+z)^a\right],
\qquad P_n=P_{n,1}.
\]

These are formal-period polynomials. At omega=2pi i, integer b makes exp(b omega)=1, so the displayed formal rows agree with the actual exponential endpoint rows.

## 2. Exact scaling and the lowest coefficient

Changing local variable x=bz in each entry gives

\[
[z^s]e^{b\varepsilon z}(jZ+z)^a
=b^{s-a}[x^s]e^{\varepsilon x}(jbZ+x)^a,
\quad\varepsilon\in\{0,1\}.
\]

The sums of row time indices and column spatial indices both equal n(n-1). The determinant factors therefore cancel exactly:

\[
\boxed{P_{n,b}(Z)=P_n(bZ).}
\tag{1}
\]

In particular, this is not a rescaling with an extra hidden power of b.

The known finite commutator formula is P_n(Z)=det(E exp(ZD)-exp(ZD)E), with E truncated multiplication by exp(x) and D differentiation on polynomials of degree below n. Its boundary identity is DE-ED-E=-nuv^t, where u is the last coordinate vector and v_a=1/(n-a)!.

The first-order commutator is ED-DE=-E+nuv^t. Since det(E)=1, E^(-1)u=u and v^tu=1, the determinant lemma gives

\[
\boxed{P_n(Z)=(-1)^n(1-n)Z^n+O(Z^{n+1}).}
\tag{2}
\]

Hence P_n is nonzero for every n>=2 and its lowest degree is exactly n. Write delta_n=deg(P_n). Coefficient extraction already gives delta_n<=n(n-1). The companion [endpoint theorem](../denominator-reflection-route/ENDPOINT_CONTENT_AND_DENOMINATORS.md) proves the sharper exact values delta_n=n^2/2 for even n and delta_n=(n^2-3)/2 for odd n. In particular delta_n<=n^2/2.

## 3. Integer denominator and primitive multiplier are different

Let d_n be the least positive integer clearing all rational coefficients of P_n. Let h_n be the gcd of the integer coefficients of d_nP_n. Define

\[
R_n(Z)=\frac{d_n}{h_n}P_n(Z)=\sum_k a_kZ^k\in\mathbb Z[Z],
\qquad\gcd_k a_k=1,
\]

with multiplier c_n=d_n/h_n>0. This definition does not presume h_n=1 beyond verified cases.

For integer b put

\[
g_n(b)=\gcd_k(a_k b^k).
\]

Then the primitive integer polynomial associated with the dilated determinant is exactly

\[
\boxed{R_{n,b}(Z)=\frac{R_n(bZ)}{g_n(b)}
=\frac{c_n}{g_n(b)}P_{n,b}(Z).}
\tag{3}
\]

The positive multiplier c_n/g_n(b) can be rational. It is not necessarily the least positive integer denominator multiplier. If h_n=1, the latter is

\[
d_{n,b}=\frac{d_n}{\gcd(d_n,g_n(b))},
\]

and d_{n,b}P_{n,b} may still have a nontrivial coefficient gcd. For general h_n one can compute d_{n,b} directly as the lcm of the reduced coefficient denominators of P_n(bZ).

Because the lowest nonzero degree is n, b^n divides every coefficient of R_n(bZ). Conversely a Bezout combination of the primitive coefficients a_k, after multiplication by b^(delta_n-k), shows that g_n(b) divides b^delta_n. Thus

\[
\boxed{b^n\mid g_n(b)\mid b^{\delta_n},
\qquad\frac{c_n}{b^{\delta_n}}\le\frac{c_n}{g_n(b)}\le\frac{c_n}{b^n}.}
\tag{4}
\]

All primes in g_n(b) divide b. More precisely,

\[
v_p(g_n(b))=\min_{k:a_k\ne0}\bigl(v_p(a_k)+k v_p(b)\bigr).
\tag{5}
\]

Equations (3)--(5) retain the exact arithmetic gain instead of treating a dilation as denominator-free.

## 4. Actual-period upper bounds and growing-frequency decay

Let U_n=(2pi)^(n^2) F_n^2/F_(2n), where F_n=product_(a=0)^(n-1) a!. The preceding Fourier--Hankel identity implies for every real t

\[
|P_n(2\pi i b)|
=b^{n^2}U_n\,|\chi_n(2\pi b)|,
\qquad |\chi_n(t)|\le1,
\tag{6}
\]

where chi_n(t)=E cos(tX_n) for the centered sum under the Vandermonde-squared probability density on [0,1]^n. The identity uses b integer to identify the formal and actual endpoint determinants. The previously proved lower bound chi_n(2pi)>1/6 does not automatically extend to b>1.

Combining (3), (4) and (6) gives the unconditional upper estimate

\[
\boxed{|R_{n,b}(\omega)|\le c_nU_n b^{n^2-n}.}
\tag{7}
\]

For comparison with b=1, the exact ratio is

\[
\frac{|R_{n,b}(\omega)|}{|R_n(\omega)|}
=\frac{b^{n^2}}{g_n(b)}
\frac{|\chi_n(2\pi b)|}{\chi_n(2\pi)}.
\tag{8}
\]

Without an additional lower bound on the numerator characteristic function, (8) alone does not compare the primitive values. The next section supplies such a bound asymptotically for each fixed b by applying an existing central limit theorem. The finite-n positivity estimate at b=1 cannot be silently reused for b>1.

A consequence of the companion endpoint denominator theorem makes the joint scales explicit. Put r=floor(n/2), s=ceil(n/2) and D'_n=F_n^2/(F_r^2 F_s^2). The theorem proves d_n<=D'_n in every dimension. Since c_n=d_n/h_n<=d_n,

\[
\ln(c_nU_n)\le\ln(D'_nU_n)
=-\tfrac12n^2\ln n+O(n^2).
\]

Hence every fixed integer b preserves decay of R_(n,b)(omega). More generally, for any fixed 0<epsilon<1/2 and any positive integer sequence b_n<=n^(1/2-epsilon), equation (7) gives

\[
\boxed{\log|R_{n,b_n}(\omega)|
\le-\epsilon n^2\log n+O(n^2)\longrightarrow-\infty.}
\]

This holds along all integer dimensions, and does not require an exact formula for the odd-dimensional primitive content. The constant in the O(n^2) bound comes from the denominator theorem and is independent of the chosen b_n in the displayed range. The legal height H_(b_n) from section 1 must still be used.

This is an upper estimate for a transcendental-period value. A value tending to zero here is not an integer contradiction or a rational approximation result.

## 5. Fixed integer frequency b>1 eventually worsens the primitive value

The probability law used in (6) is the orthogonal-polynomial ensemble for the uniform probability measure on [0,1]: its orthogonal-polynomial determinant is a constant times the Vandermonde, so the joint density is the same Vandermonde squared. Apply Theorem 2.6 of Breuer and Duits, *Central Limit Theorems for Biorthogonal Ensembles and Asymptotics of Recurrence Coefficients*, [arXiv:1309.6224v2](https://arxiv.org/pdf/1309.6224v2), with f(u)=u. Its assumptions hold: the measure is compactly supported, its essential support is the single interval [0,1], the density is positive almost everywhere there, and f is real-valued C^1 on the real line. Their variance uses the positive Fourier indices of f((1+cos(theta))/2), so only f-hat_1=1/4 contributes, giving variance 1/16. Reflection u_j -> 1-u_j gives exact mean n/2. Consequently

\[
X_n=\sum_{j=1}^n u_j-n/2\;\Longrightarrow\;N(0,1/16),
\qquad\chi_n(t)\longrightarrow e^{-t^2/32}
\tag{9}
\]

for each fixed real t. The characteristic-function convergence follows directly from convergence in distribution. These are an application of the cited theorem, not a new central limit theorem.

There is an independent check of the variance convention. For orthonormal shifted Legendre polynomials, multiplication by u has diagonal 1/2 and off-diagonal coefficient a_k=k/(2 sqrt(4k^2-1)). The determinantal projection variance of the sum is Tr(P_n J^2 P_n-(P_n J P_n)^2)=a_n^2, because only the recurrence edge crossing the n-dimensional projection survives. Thus its exact variance is n^2/[4(4n^2-1)], tending to 1/16. This agrees with the separately proved second-moment formula in the preceding Hankel note and rules out a factor-two convention error.

Fix an integer b>=2. Equation (8) therefore gives the exact asymptotic relation

\[
\frac{|R_{n,b}(\omega)|}{|R_n(\omega)|}
=\frac{b^{n^2}}{g_n(b)}
\left(e^{-\pi^2(b^2-1)/8}+o_b(1)\right).
\tag{10}
\]

Thus for all sufficiently large n, with the threshold allowed to depend on b,

\[
\boxed{\frac{|R_{n,b}(\omega)|}{|R_n(\omega)|}
\ge\tfrac12e^{-\pi^2(b^2-1)/8}b^{n^2-\delta_n}
\ge\tfrac12e^{-\pi^2(b^2-1)/8}b^{n^2/2}\longrightarrow\infty.}
\tag{11}
\]

This compares each fixed b>1 to b=1 in growing dimension, after the best primitive coefficient normalization. Both values tend to zero by section 4, while their ratio tends to infinity. The comparison is not uniform in b and says nothing by itself about a changing frequency b_n. The theorem supplies existence of the fixed-b threshold, not a numerical value for it. Even the coarser degree bound delta_n<=n(n-1) would suffice for a weaker lower bound proportional to b^n.

## 6. Sending frequency to infinity at fixed dimension cannot make primitive values small

Fix n. From (3)--(4),

\[
|R_{n,b}(\omega)|\ge\frac{|R_n(b\omega)|}{b^{\delta_n}}.
\]

The right side tends to |a_(delta_n)| (2pi)^delta_n as b tends to infinity. Therefore

\[
\boxed{\liminf_{b\to\infty}|R_{n,b}(\omega)|
\ge |a_{\delta_n}|(2\pi)^{\delta_n}>1.}
\tag{12}
\]

This holds with the best primitive normalization, not merely a chosen row clearing. It is a fixed-n statement and does not control b growing together with n.

There is an explicit finite threshold. Put L_n=sum_(k<delta_n)|a_k|. If b>=max(1,ceil(L_n/(3|a_delta|))), then 2pi b>6b>=2L_n/|a_delta|, so the leading-term triangle estimate yields

\[
|R_{n,b}(\omega)|>\tfrac12|a_{\delta_n}|(2\pi)^{\delta_n}>1.
\]

This sufficient threshold need not be sharp.

## 7. Rational centers and scope

At a rational imaginary center Z=ip/q, q^delta_n R_(n,b)(Z) is a Gaussian integer. Since Z^n is a known factor, one may instead apply the degree delta_n-n clearing to R_(n,b)(Z)/Z^n when p!=0. Any comparison with approximation error must retain the relevant q power and the chosen column-height budget H_b.

If gcd(b,q)>1, the argument b p/q has a reduced denominator; the same exact polynomial relation (3) still governs coefficient content. This reduction does not by itself control the value at b omega or its difference from a rational-center value.

The results separate three issues: frequency dilation, primitive coefficient normalization, and actual-period magnitude. They do not establish the optimal joint frequency/dimension choice, an infinite near-root theorem, or a new statement about pi's bad approximability.

## Verification

`verify_frequency_dilation.py` reconstructs P_n from the exact Laurent-Hankel identity, checks independent direct raw matrices at several centers and integer frequencies, and verifies primitive content, p-adic minima, integer denominators, and the explicit fixed-dimension threshold on small instances. It uses only the Python standard library. The infinite statements are proved above, with the fixed-b asymptotic explicitly relying on the cited Breuer--Duits theorem; finite checks are not used to extrapolate them.
