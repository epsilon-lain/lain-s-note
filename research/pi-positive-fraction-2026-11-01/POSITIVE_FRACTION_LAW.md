# Positive-fraction boundary transport: a uniform quadratic fraction law

2026-11-01. This note treats a positive fraction of the full source rows, rather than the lower-order first-exit examples. It proves a uniform displacement lower bound and a family of nonzero transport coefficients with comparable dimensional order. It does not bound the accompanying period-function minors or decide bad approximability of pi.

All weights and the dimension are fixed before H tends to infinity. We use rational positive weights and admissible H along a common lattice of their denominators. Subsequent large-weight or dimension limits are explicitly separate.

## 1. Sharp cost minimum in one unrestricted weighted bucket

Let c=(c_1,...,c_m) have positive entries, let h>0, and define

\[
B_h=\{\beta\in\mathbb Z_{\ge0}^m:c\cdot\beta<h\},
\qquad b_h=|B_h|,\qquad V=\frac1{m!\prod_i c_i}.
\]

Order all lattice points by their cost c dot beta, resolving ties arbitrarily, and let L(n) be the sum of the n cheapest costs, with L(0)=0. For target sets A of cardinality b_h with exactly k points outside B_h, the smallest possible cost increase is exactly

\[
\boxed{D_h(k)=L(b_h+k)+L(b_h-k)-2L(b_h).}       \tag{1}
\]

Indeed, retain the b_h-k cheapest inside points and take the k cheapest outside points. Every other selection has at least that cost. This assertion concerns sets of lattice points, not their interpolation rank.

Uniformly for 0<=k<=b_h as h tends to infinity,

\[
\boxed{D_h(k)=b_h h F_m(k/b_h)+O((h+1)^m),}    \tag{2}
\]

where

\[
F_m(\rho)=\frac m{m+1}
\left[(1+\rho)^{1+1/m}+(1-\rho)^{1+1/m}-2\right],
\quad0\le\rho\le1.
\]

The implied constant depends on the fixed weights and m, not on k. To justify the error uniformly, comparison with unit lattice cubes gives
\(b_h=Vh^m+O((h+1)^{m-1})\), and weighted integration over those cubes gives
\(L(n)=\frac m{m+1}V^{-1/m}n^{1+1/m}+O(n+1)\).
Substitute n=b_h+k,b_h-k,b_h into (1). Since (b_h/V)^(1/m)=h+O(1), the resulting error is O((h+1)^m), including k=0 and k=b_h.

For clarity, the L(n) estimate also follows by locating the nth cost within O(1) of (n/V)^(1/m): a unit cube changes cost by at most sum c_i, so the counts sandwich each quantile within a fixed-width cost interval. Integrating the cost over the corresponding simplex yields the main term, with O(n+1) boundary and cube errors.

The useful quadratic inequality is

\[
\boxed{F_m(\rho)\ge\rho^2/m.}                 \tag{3}
\]

For 0<=rho<1,
\(F_m''(\rho)=m^{-1}[(1+\rho)^{1/m-1}+(1-\rho)^{1/m-1}]\ge2/m\),
by convexity of x^(1/m-1); for m=1 it is equality. Integrate twice using F_m(0)=F_m'(0)=0, then pass continuously to rho=1. Combining (2)--(3) gives, uniformly in k,

\[
D_h(k)\ge\frac{k^2}{mVh^{m-1}}-O((h+1)^m).   \tag{4}
\]

The admissible h values are bounded away from zero by a fixed lattice unit; finitely many small values are absorbed in the constant. No assertion of a dimension-free k^2 bound is made.

## 2. A lower bound for a positive fraction of the complete packet

Use the source packet

\[
\mathcal R_H=\{(j,s,\beta):0\le j<K,\quad
v_0s+\sum_i v_i\beta_i<H\},
\quad v_i=w_i/\theta>v_0,\quad c_i=v_i-v_0.
\]

Transport preserves the center j and total order N=s+|beta|. In bucket N put h_N=H-v_0N; its source types are
\(B_N=\{\beta:c\beta<h_N,\ |\beta|\le N\}\).
Select exactly |B_N| distinct targets in the same bucket, obeying |a|<=N. Let k_{j,N} be the number outside the original source bucket, and let

\[
k_H=\sum_{j,N}k_{j,N},\quad
M_H=|\mathcal R_H|,\quad
\Gamma=\sum_{\rm targets}a-\sum_{\rm source}\beta,
\quad
\mathcal A=\prod_i\frac{v_i}{c_i}.
\]

If k_H/M_H tends to eta in [0,1], then

\[
\boxed{\liminf_{H\to\infty}\frac{c\cdot\Gamma}{M_HH}
\ge\frac{\eta^2}{(m+1)\mathcal A}.}             \tag{5}
\]

This holds for every such target selection, even before requiring a nonzero Pascal determinant.

**Proof.** In a capped bucket the removed inside points belong to the unrestricted B_{h_N}; the outside targets have cost at least h_N, since they obey the same total-order cap. Thus the displacement is at least D_{h_N}(k_{j,N}). Sum (4), then apply Cauchy--Schwarz:

\[
\sum_{j,N}k_{j,N}^2h_N^{1-m}
\ge\frac{k_H^2}{\sum_{j,N}h_N^{m-1}}.
\]

The source count and denominator satisfy

\[
M_H\sim\frac{KH^{m+1}}{(m+1)!v_0\prod_i v_i},
\qquad
\sum_{j,N}h_N^{m-1}\sim\frac{KH^m}{mv_0}.
\]

The summed O((h_N+1)^m) error is O(H^{m+1})=o(M_HH). Substituting these constants gives exactly (5). The total-order cap was never ignored in the source row count. It only enlarged the comparison problem used for a lower bound.

Suppose now every grouped Pascal coefficient is nonzero, all |zeta_i| lie in (0,1), and every target satisfies wa<=H. The errors zeta_i, as well as the weights, are fixed before H tends to infinity; in particular the quantity tau below is fixed in that limit. Set

\[
\tau=\min_i\frac{-\log|\zeta_i|}{c_i},\qquad
w_*=\min_i w_i.
\]

Nonzero Pascal coefficients imply Gamma_i>=0, so the error-monomial cost is at least tau times c dot Gamma. For a target a, the Euclidean norm of its Pascal column is at most 2^{|a|}, because
\(\sum_{\beta\le a}\binom a\beta^2=\prod_i\binom{2a_i}{a_i}\le4^{|a|}\).
Hadamard, applied separately to the bucket determinants, therefore gives the uniform full transport bound

\[
\boxed{\limsup_{H\to\infty}\frac{\log|\det T_{L,J}|}{M_HH}
\le\frac{\log2}{w_*}-\frac{\tau\eta^2}{(m+1)\mathcal A}.} \tag{6}
\]

There is no factorial counting charge in this particular Pascal norm estimate. Equation (6) bounds the transport coefficient only; the period minor multiplying it must still be estimated.

## 3. Positive-fraction nonzero constructions at the same dimensional order

The lower bound does not restore a dimension-free square law. Fix 0<rho<=1. In an unrestricted bucket B_h, split into first-coordinate fibers y, of lengths
\(n_y=\lceil(h-c' y)/c_1\rceil\).
Set k_y=floor(rho n_y), retain the first n_y-k_y points and take the first k_y outside points. The fiber Pascal product is nonzero, and the weighted displacement is exactly

\[
c\cdot\Gamma=c_1\sum_y k_y^2.
\]

Lattice integration gives

\[
\sum_y n_y^2\sim\frac{2}{(m+1)c_1}b_hh,
\qquad \sum_y k_y\sim\rho b_h,
\]

hence

\[
\boxed{\frac{c\cdot\Gamma}{b_hh}\longrightarrow
\frac{2\rho^2}{m+1}.}                           \tag{7}
\]

The error from the floor in k_y is O(sum n_y)=O(h^m), lower order than h^(m+1). For m=1 this agrees with the sharp bucket value F_1(rho)=rho^2. For m>1 it is a constructive upper bound within the nonzero-Pascal class, rather than an assertion that the unrestricted minimum F_m is attained by independent target sets.

Here is a complete-packet construction respecting the total-order cap. Let c_* = min_i c_i and apply the fiber construction only to buckets with

\[
N\ge\frac{(1+\rho)H+\rho c_1}{c_*+(1+\rho)v_0};             \tag{8}
\]

use identity transport in the remaining buckets. Every chosen target has
\(ca<(1+\rho)h_N+\rho c_1\) and therefore |a|<=N by (8). Assume theta(1+rho)<=1. For all sufficiently large H the column-type cap wa<=H also holds: in the selected buckets

\[
wa\le\theta(1+v_0/c_*)((1+\rho)h_N+\rho c_1),
\]

whose leading maximum divided by H is
\(\theta(1+\rho)(1+v_0/c_*)/(1+(1+\rho)v_0/c_*)<1\).
All omitted small-bucket targets are original source types and already fit.

Define
\(t_0=c_*/(c_*+(1+\rho)v_0)\).
For this complete nonzero transport selection, direct bucket integration yields

\[
\boxed{
\frac{k_H}{M_H}\longrightarrow
\rho\mathcal A t_0^{m+1},\qquad
\frac{c\cdot\Gamma}{M_HH}\longrightarrow
\frac{2\rho^2}{m+2}\mathcal A t_0^{m+2}.}        \tag{9}
\]

The construction can be made at every center. The resulting grouped transport coefficient is nonzero because it is a product of nonzero fiber Pascal determinants and identity blocks. This does not certify a nonzero accompanying period minor, and centers may still occupy the same transverse types.

After the H limit, let c_*/v_0 tend to infinity at fixed m. Then mathcal A and t_0 tend to one. The moved fraction tends to rho, and its cost tends to 2 rho^2/(m+2). Combined with (5), this puts the positive-fraction displacement scale between rho^2/(m+1) and 2 rho^2/(m+2) in this subsequent limit. Both are of order 1/m.

For endpoint-like quality -log|zeta_i|/c_i approaching 2 theta, the corresponding error powers alone remain of order theta/m even though a positive fraction moves. This statement requires those quality assumptions; it is not inferred from non-bad-approximability alone.

## 4. What is left for the main contradiction

A positive moved fraction does force an MH-scale quadratic-fraction penalty at every fixed dimension. It does not force an order-one penalty uniform in dimension: explicit nonzero transport coefficients already demonstrate 1/m dilution. The construction in (9) deliberately says nothing about the period-minor magnitude or collision savings across centers. Those can be decisive and have not been paid for here.

The next estimate must combine this fraction law with the occupancy/independence of the actual period rows, or exploit cancellation between complete target sets. Neither the scalar fraction law nor the local nonzero coefficient construction closes the endpoint comparison.
