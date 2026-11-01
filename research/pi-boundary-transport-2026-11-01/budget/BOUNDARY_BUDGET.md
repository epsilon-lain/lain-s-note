# Boundary coefficient budget and the remaining determinant estimate

2026-11-01. Companion to the exact boundary transport identity. This is an algebraic reformulation and an explicit coefficient bound, not a proof of bad approximability of pi. All asymptotics fix dimension, weights, centers and radius before H tends to infinity.

## 1. Setup and exact coefficients

Use complete source rows rho=(j,s,beta), with cost

\[
\ell_\rho=v_0s+\sum_i v_i\beta_i<H,\qquad
v_i=w_i/\theta,\qquad c_i=v_i-v_0>0.
\]

Let omega be nonzero, let every r_i be nonzero, and put

\[
b_i=r_i/\omega,\quad \zeta_i=1-b_i^{-1},\quad
 g_\rho=H-\ell_\rho>0.
\]

For raw monomial columns Y^h X^alpha, rescale columns by b^alpha and rows by b^{-beta}. Change complete logarithmic jets to ordinary local jets. Internal triangular transport can then be removed with determinant one. The surviving coefficient replacing row (j,s,beta) by (j,s-|gamma|,beta+gamma) outside the packet is

\[
W_{\rho,\gamma}
=\binom{\beta+\gamma}{\beta}(-\zeta)^\gamma
\sum_{\substack{0\le\kappa\le\gamma\\c\cdot\kappa<g_\rho}}
(-1)^{|\kappa|}\binom\gamma\kappa,
\qquad |\gamma|\le s,\quad c\cdot\gamma\ge g_\rho.       \tag{1}
\]

This coefficient formula includes the cancellations from inverting the internal triangular block. It is not the original termwise translation expansion. The companion transport note proves the exact identity. Types with w(beta+gamma)>H act as zero on all allowed columns and can be discarded.

## 2. Explicit suppression at the first budget exit

Assume 0<|zeta_i|<1 and set

\[
\tau=\min_i\frac{-\log|\zeta_i|}{c_i},\qquad
w_*=\min_iw_i,\qquad
N_A=\#\{a\ge0:wa\le H\}.
\]

For every surviving boundary coefficient,

\[
\boxed{|W_{\rho,\gamma}|
\le 3^{H/w_*}\exp(-\tau g_\rho).}                 \tag{2}
\]

Indeed, the absolute value of the alternating sum is at most 2^{|gamma|}. For each coordinate, the term binom(beta_i+gamma_i,beta_i)2^{gamma_i} is one nonnegative term in (1+2)^{beta_i+gamma_i}. Their product is at most 3^{|beta+gamma|}. The column support gives |beta+gamma|<=H/w_*. Finally,

\[
|\zeta^\gamma|
\le e^{-\tau c\cdot\gamma}\le e^{-\tau g_\rho}.
\]

There are at most N_A relevant gamma for one source row, so the sum of its absolute boundary coefficients is at most N_A times the right side of (2). This controls coefficients only. It does not bound the sizes of the outside derivative rows or their mixed determinants.

At the true period omega=2 pi i and r_i=2 i p_i/q_i,

\[
|\zeta_i|=\frac{|\pi-p_i/q_i|}{|p_i/q_i|}.
\]

Writing the exact approximation quality as
\(\Lambda_i=-\log(q_i^2|\pi-p_i/q_i|)\), one obtains

\[
-\log|\zeta_i|=2\log q_i+\Lambda_i+\log|p_i/q_i|.   \tag{3}
\]

Thus the endpoint quality enters explicitly. When w_i=ceil(log q_i), v_0/w_* tends to zero, and Lambda_i/w_i tends uniformly to zero, tau tends to 2 theta. This is a statement under these additional limiting assumptions; non-bad-approximability alone does not prescribe the rate of Lambda_i relative to log q_i. If the errors are smaller, formula (3), rather than this limit, is the appropriate budget.

## 3. Exact mixed-minor target

Let B be the ordinary-jet matrix at the common period, on the same raw columns. Let B^{S,gamma} be obtained from its source row list by replacing precisely the rows rho in S with their chosen outside rows, preserving row order. Put

\[
D_i=\sum_{\rm columns}\alpha_i-\sum_{\rm source\ rows}\beta_i.
\]

The transport identity and determinant multilinearity give the exact finite expansion

\[
\boxed{\det A(r)=b^D
\sum_{S\subseteq\mathcal R_H}
\sum_{\text{allowed boundary }\gamma_\rho,\ \rho\in S}
\left(\prod_{\rho\in S}W_{\rho,\gamma_\rho}\right)
\det B^{S,\gamma}.}                               \tag{4}
\]

The empty subset contributes the original period determinant. Repeated replacement row labels force zero; other algebraic dependencies may force further zeros. Raw-column arithmetic scaling can be carried through (4); arbitrary mixed columns require expanding or carrying their full diagonal scaling, and must not be assigned a fictitious single alpha.

For a nonzero raw minor, D_i>=0 for every i: a nonzero determinant term must match each row beta to a column alpha>=beta. Also sum_i w_i D_i<=MH. Consequently

\[
\frac{\log|b^D|}{MH}
\le \max_i\frac{(\log|b_i|)_+}{w_i}.              \tag{5}
\]

This diagonal rescaling is a small, explicitly payable cost for genuinely close centers. It should not be confused with the potentially large boundary-minor cost.

For fixed parameters, the number of expanded terms is at most (1+N_A)^M. From (2)--(4), with the logarithm of zero interpreted as minus infinity,

\[
\frac{\log|\det A(r)|}{MH}
\le \max_i\frac{(\log|b_i|)_+}{w_i}
+\frac{\log3}{w_*}+\frac{\log(1+N_A)}H
+\max_{S,\gamma}
\left\{\frac{\log|\det B^{S,\gamma}|}{MH}
-\frac{\tau}{MH}\sum_{\rho\in S}g_\rho\right\}.    \tag{6}
\]

Here log(1+N_A)/H tends to zero. Formula (6) supplies an explicit sufficient target, not a verified favorable bound. One may retain exact W in (4) for a sharper estimate, including inter-term cancellation.

## 4. Why the coefficient saving alone is insufficient

Weighted simplex counting gives

\[
\frac1{MH}\sum_{\rho\in\mathcal R_H}g_\rho
\longrightarrow\frac1{m+2}.                     \tag{7}
\]

There are m+1 nonnegative weighted coordinates, and their simplex average cost is (m+1)H/(m+2). The K center copies cancel in the normalized mean.

For the all-boundary part of (4), bound (2) supplies the suppression exponent tau/(m+2), before coefficient combinatorics and the outside determinants are counted. At endpoint quality tau approximately 2 theta, this is only 2 theta/(m+2). A mixed term receives only the gaps of its replaced rows; the empty term receives none. This calculation does not bound the actual determinant from below or prove that boundary cancellation cannot help. It shows that a proof cannot stop at the coefficient powers in (2).

The remaining concrete task is to bound the coupled mixed minors in (4), with their actual pattern (s,beta)->(s-|gamma|,beta+gamma), sharply enough to beat the arithmetic lower bound on the same raw minor. Applying the zero-error full consecutive-jet Cauchy product unchanged to these altered row lists would be unjustified. Arbitrary independent Hardy components would also discard relations inherent in the polynomial columns.

## 5. A raw-minor contradiction criterion with exact arithmetic cost

Keep the earlier tail-cancellation assumptions and minimal logarithmic truncations
\(T_i=\lceil w_i/(\theta v_0)\rceil\),
\(L_i=\operatorname{lcm}(1,\ldots,T_i-1)\), with empty lcm equal to one. For a nonzero raw rational-center minor, multiplying column alpha by \(\prod_i(q_iL_i)^{\alpha_i}\) and row beta by \(\prod_i(q_iL_i)^{-\beta_i}\) gives a nonzero Gaussian-integer determinant. Tail cancellation transfers its determinant back to the full-logarithm matrix. Consequently

\[
|\det A(r)|\ge\prod_i(q_iL_i)^{-D_i}.             \tag{8}
\]

Combining (8) with (4) gives the necessary lower bound

\[
\left|\sum_{S,\gamma}
\Bigl(\prod_{\rho\in S}W_{\rho,\gamma_\rho}\Bigr)
\det B^{S,\gamma}\right|
\ge\prod_i(q_iL_i|b_i|)^{-D_i}.                  \tag{9}
\]

Define the explicit raw-minor cost

\[
\mathfrak A_H=\frac1{MH}\sum_iD_i\log(q_iL_i|b_i|).
\]

A sufficient route to a contradiction is a strict eventual inequality

\[
\max_{S,\gamma}
\left\{\frac{\log|\det B^{S,\gamma}|}{MH}
-\frac\tau{MH}\sum_{\rho\in S}g_\rho\right\}
+\frac{\log3}{w_*}+\frac{\log(1+N_A)}H
< -\mathfrak A_H.                               \tag{10}
\]

This is a newly explicit proof obligation, not a proved inequality. It uses the same nonzero raw minor on both sides, requires no common mixed-column selection, and charges the arithmetic denominator through the actual degree differences D_i rather than a uniform maximum. Any analytic bound for the period mixed minors must still include their genuine column growth; (10) does not make it free. For the source's critical parameters the existence of the raw nonzero minor retains the manuscript's interpolation hypothesis. No independent verification of that hypothesis is claimed here.
