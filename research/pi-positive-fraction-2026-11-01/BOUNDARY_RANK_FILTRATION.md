# Exact low-degree boundary compression and growth of the full transport rank

2026-11-01. This note studies the complete boundary remainder rather than an isolated degree-two term. The statements concern formal degree and matrix rank; they supply no norm estimate for an inverse period matrix and no Diophantine conclusion.

## 1. Every fixed error degree lies in a thin source strip

Let the source rows be rho=(j,s,beta), with

\[
\ell_\rho=v_0s+\sum_i v_i\beta_i<H,
\qquad c_i=v_i-v_0\ge0,
\qquad g_\rho=H-\ell_\rho>0.
\]

Assume c_max=max_i c_i>0. The exact boundary coefficient of degree gamma can occur only if c dot gamma>=g_rho. Consequently its total error degree obeys

\[
|\gamma|\ge d_\rho:=\left\lceil\frac{g_\rho}{c_{\max}}\right\rceil.
\tag{1}
\]

This includes coefficients that subsequently vanish by the alternating sum; such cancellations can only raise the order. If every c_i=0, the full remainder vanishes and no strip argument is needed.

For integer D>=0, retain only error monomials of total degree at most D. Write E_{<=D}(zeta) for the corresponding remainder matrix on any chosen column family. Its nonzero rows lie in

\[
S_D(H)=\{\rho:0<g_\rho\le c_{\max}D\}.
\]

An explicit count is

\[
\boxed{|S_D(H)|\le K\left(\left\lceil\frac{c_{\max}D}{v_0}\right\rceil+1\right)
\prod_{i=1}^m\left(\left\lfloor\frac H{v_i}\right\rfloor+1\right).} \tag{2}
\]

For each beta, the admissible time orders lie in an interval of length c_max D/v_0; the displayed deliberately nonsharp integer bound follows. At fixed dimension, weights and D this is O(H^m), whereas the complete row count M(H) is asymptotic to
\(K H^{m+1}/((m+1)!v_0\prod_i v_i)\).
Thus a fixed-degree perturbation acts through a vanishing fraction of the source rows.

## 2. An exact reduced determinant for all coefficients through degree D

Let B be a nonsingular M by M period matrix on the chosen columns. Let J_D be the M by r coordinate-inclusion matrix for S_D(H), where r=|S_D(H)|. Write E_{<=D}=J_D F_D. The matrix determinant lemma gives, in the ring of formal error series modulo terms of total degree greater than D,

\[
\boxed{
\det(B+\eta E(\zeta))
\equiv\det B\;\det\left(I_r+\eta F_D(\zeta)B^{-1}J_D\right)
\pmod{(\zeta)^{D+1}}.}                          \tag{3}
\]

The marker eta counts the number of boundary rows used. Equation (3) follows by first discarding E's error degrees above D, which cannot affect a coefficient of total degree at most D, and then applying det(I+UV)=det(I+VU). It is an exact algebraic reduction, not an approximation in a numerical norm. It permits all such coefficients to be computed on r=O(D H^m) rows. The inverse B^{-1} is still present; no conditioning or smallness has been inferred.

A simple additional divisibility statement follows directly by determinant multilinearity. Under a common scaling zeta_i=t xi_i, let d_(1)<=...<=d_(M) be the sorted row lower orders from (1). For every k,

\[
\boxed{\operatorname{ord}_t[\eta^k]\det(B+\eta E(t\xi))
\ge\sum_{r=1}^k d_{(r)}.}                      \tag{4}
\]

This does not require B to be nonsingular: every term chooses k distinct source rows from E, and each chosen row supplies its own lower power. Extra cancellations can increase the order.

For a fraction k/M(H)->rho in [0,1], weighted-simplex counting converts the guaranteed lower bound in (4) to

\[
\lim_{H\to\infty}\frac1{MH}\sum_{r=1}^k d_{(r)}
=\frac1{c_{\max}}\left[
\rho-\frac{m+1}{m+2}\left(1-(1-\rho)^{(m+2)/(m+1)}\right)\right].
\tag{5}
\]

Indeed, normalized row slack g/H has density (m+1)(1-x)^m on [0,1]. Integrating x over the lowest rho fraction gives the bracketed expression; the ceiling changes the sum by at most M, negligible after division by MH. At small rho the bracket is rho^2/(2(m+1))+O(rho^3); at rho=1 it is 1/(m+2). This is the slack-based guaranteed order, not an exact valuation of the determinant. It does not include outside-row growth or column norm costs.

## 3. The full boundary transport need not have small rank

The fixed-D compression must not be extrapolated to the entire error polynomial. An explicit complete-packet rank formula shows the distinction.

Take one transverse coordinate, v_1>v_0>0, c=v_1-v_0, and w_1=theta v_1 with 0<theta<=1/2. At each center and conserved total order N=s+beta, put

\[
b_N=\min\left(N+1,\left\lceil\frac{H-v_0N}{c}\right\rceil\right)>0
\quad(0\le N<H/v_0).
\]

The source beta indices are 0,...,b_N-1. Impose the genuine polynomial type cap a<=A_H=floor(H/w_1). The available outside indices are b_N,...,min(N,A_H). The one-coordinate Cauchy product formula for R proves

\[
\boxed{\operatorname{rank}R_N
=\min\left(b_N,\max(0,\min(N,A_H)-b_N+1)\right)}
\tag{6}
\]

for every nonzero error zeta. All square minors of maximal possible size are nonzero. This is rank over the complex numbers after nonzero specialization, or rank over the field of rational functions in zeta. It is not automatically the rank of R_N times a particular compatible column matrix.

For all sufficiently large H the cap A_H does not alter the simpler value
\(\min(b_N,N+1-b_N)\). One sufficient threshold is H>=v_1(v_1-v_0)/v_0. To see this, the cap is inactive when N<=A_H. When N>A_H, use A_H>=floor(2H/v_1), hence N>2H/v_1, and

\[
b_N\le\frac H{v_1}-\frac{Hv_0}{v_1c}+1\le\frac H{v_1}\le\frac{A_H+1}{2}.
\]

Thus there are at least b_N outside columns, so both the capped and uncapped rank equal b_N. If the displayed intermediate expression is negative, that N has no source rows and is omitted.

## 4. Its asymptotic rank can approach the whole packet

At fixed v_0,v_1,theta and number K of centers, the complete source count is

\[
M(H)\sim\frac{K H^2}{2v_0v_1}.
\]

Let N/H=x. The limiting rank density per center is zero for x<1/v_1, is (v_1x-1)/c between 1/v_1 and 2/(v_1+v_0), and is (1-v_0x)/c between 2/(v_1+v_0) and 1/v_0. Integrating these three linear pieces gives

\[
\operatorname{rank}R(H)
\sim\frac{K(v_1-v_0)H^2}{2v_0v_1(v_1+v_0)}.
\]

In particular,

\[
\boxed{\frac{\operatorname{rank}R(H)}{M(H)}
\longrightarrow\frac{v_1-v_0}{v_1+v_0}.}        \tag{7}
\]

The fixed-weight limit is taken first. Across weight choices with v_1/v_0 large, this limiting fraction approaches one. Hence the full boundary operator can have rank proportional to the entire source packet even though every fixed-degree truncation has rank O(H) in this one-coordinate case.

## 5. What is still missing

Equations (2)--(4) give a controlled algebraic way to compute and organize low-order aggregate coefficients. Equation (7) rules out treating the full remainder as a uniformly low-rank perturbation merely because its first-order part is low rank. Neither says that high-degree terms dominate after substituting the actual tiny errors: their coefficients and compatible period minors might compensate or cancel. A proof of the endpoint comparison still needs their quantitative total contribution, with the inverse-period or mixed-minor cost explicitly accounted for.

The exact verifier tests boundary coefficient support, formal truncated determinant reduction, and the Cauchy-block rank formula with the column type cap. It also checks finite rank counts against the proved asymptotic formula; no infinite assertion is inferred from finite numerical convergence.
