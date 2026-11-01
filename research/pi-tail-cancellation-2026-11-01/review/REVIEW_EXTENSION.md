# Independent review of the sharper collision and exact-quality estimates

2026-11-01. This supplements `REVIEW.md`. The reviewed source notes are `SHARP_COLLISION.md`, `EXACT_QUALITY_CRITERION.md`, and the completed `TAIL_CANCELLATION.md`. The statements below concern the displayed estimates and their assumptions, not an independent verification of the manuscript's interpolation theorem.

**Verdict:** the sharper collision bound and anisotropic comparison criterion are valid as stated. The shell-quality conclusion is also valid with fixed dimension and weights before the `H` limit, followed by the specified dimension limit. No scalar estimate is a lower bound for the test determinant.

## 1. Complete summation of distinct Taylor degrees

For `0<r<1`, sorting the distinct nonnegative integers and subtracting `0,1,…,n−1` gives

\[
\sum_{d_1,\ldots,d_n\text{ distinct}}r^{\sum_i d_i}
=n!\,r^{n(n-1)/2}\prod_{k=1}^n(1-r^k)^{-1}.
\]

The factor `n!` is required because the original determinant expansion labels its rows. Its appearance does not cost a main exponential order: after grouping by type, `∏ₐnₐ!≤M!`, and `ln(M!)/(MH)≤ln M/H→0` for the fixed-parameter polynomial row counts used here.

For a Cauchy radius `R₀>7K`, the actual period tests lie in a disk of radius less than `7K`; hence `r=7K/R₀`. The determinant bound uses one factor `M!` for each coefficient determinant and at most another factor `M!` from the distinct-degree sum. With `Σₐnₐ=M`, the resulting leading collision term is

\[
-\frac12\ln(R_0/(7K))\sum_a n_a^2.
\]

The positive correction from replacing `Σₐnₐ(nₐ−1)` by `Σₐnₐ²` is exactly `c_R M`, not a main-order quadratic cost. Together with the denominator bound `∏(1−r^k)≥(1−r)^M`, this yields the stated remainder

\[
\frac{2\ln M+c_R-\ln(1-r)}H.
\]

The maximum in the new entire-function norm is also correct. The single column budget gives

\[
hR_0+|\alpha|\ln(2R_0)
\le (w_0h+w\alpha)
\max\{R_0/w_0,\ln(2R_0)/w_*\}.
\]

It is unnecessary to charge the entire budget separately to both terms. No independence of same-type tests is assumed: their dependence can only make a determinant zero and strengthen its upper bound. Mixed columns incur the explicitly retained norm multiplier and leave the collision coefficient unchanged.

## 2. Exact-quality scalar bounds and the anisotropic split

Write `νᵢ=−ln|π−pᵢ/qᵢ|/wᵢ`. The strict positivity assumptions `νᵢ>0` are essential to the stated anisotropic simplex and to discarding the nonnegative scalar saving in the collision branch. For a surviving constant-shift row at a nonzero center,

\[
|\xi|=\binom a\beta(2j)^{|a-\beta|}
e^{-\sum_i\nu_iw_i(a_i-\beta_i)}.
\]

The prefactor is at least one. It is at most

\[
2^{|a|}[2(K-1)]^{|a-\beta|}
\le[4(K-1)]^{H/w_*},
\]

which verifies both sides of the scalar estimate. At the zero center only `a=β` survives and its scalar is one.

For the arithmetic threshold, `uᵢ=ln(qᵢLᵢ)/wᵢ>0` implies

\[
R_{\rm ar}=u_{\max}-\frac{\theta}{m+2}\sum_i u_i
\ge u_{\max}\left(1-\frac{m\theta}{m+2}\right)>0.
\]

Thus no sign ambiguity occurs in the displayed contradiction criterion. Counting the enlarged low-type simplex gives

\[
D_A(H)\sim\frac{(AH)^m}{m!\prod_i(\nu_iw_i)},
\qquad
L_\nu=\frac{\eta^2K\theta^m\prod_i\nu_i}
{(m+1)v_0A^m}.
\]

Ignoring `wa≤H` only increases the counted capacity, so it produces a valid, possibly weaker, collision lower bound. In the complementary branch, more than `(1−η)M` rows have quality-weighted type above `AH`; subtraction of the original row mean gives the scalar saving `A(1−η)−Dν`. Therefore both strict inequalities in the note's criterion suffice. The polynomial number of row choices contributes `o(1)`, and strict inequalities absorb that remainder.

At fixed background data and fixed arithmetic mean of the `νᵢ`, the **source-row mean quality weight** `Dν` is fixed, whereas the collision expression contains `∏ᵢνᵢ≤ν̄ᵐ`. This is the precise limited sense in which unequal quality is unfavorable for this relaxed criterion. An arbitrary chosen displacement cost does not depend on the mean alone. The criterion does not permit replacing actual approximation qualities by freely chosen equal ones.

## 3. The mean quality in the shell witness

After changing to weighted coordinates `yᵢ=wᵢxᵢ`, the limiting source simplex and radial shell map are symmetric in the `m` coordinates. The unequal original lattice-cell shapes cause only boundary and rounding terms of lower order for fixed weights. Thus the expected weighted displacement per coordinate is the same in the limiting fractional assignment.

An independent check of its radial mean is useful. At a fixed time slice with source radius `L`, the source mean radius is `mL/(m+1)`. The shell between radii `L` and `tL` has mean radius

\[
\frac m{m+1}L\frac{t^{m+1}-1}{t^m-1}.
\]

Subtracting gives

\[
\frac m{m+1}L\frac{t^m(t-1)}{t^m-1}.
\]

The original packet's time average satisfies `E[L]/H→θ(m+1)/(m+2)`. Multiplication by the fraction `(K−1)/K` of non-anchor rows therefore yields exactly

\[
J=\frac{K-1}{K}\frac m{m+2}
\frac{t^m}{t^m-1}(T-\theta).
\]

The quality coefficient becomes `ν̄` after summing the equal coordinate expectations. Decomposing the finite fractional matching into integral injections preserves these edge marginals. Because every such injection already has independent test rows, selecting one with linear cost no greater than its expectation does not require intersecting unrelated favorable events. The deterministic occupancy bound supplies the collision term `c/v₀` for the same injection.

In the explicit `θ=1/2,C=13/10` family, `J→3/20` and `v₀→∞`. If `w_*→∞` and the mean `m^{-1}Σᵢλᵢ/wᵢ→0`, then

\[
\bar\nu
=2-\frac2m\sum_i\frac{w_i-\ln q_i}{w_i}
+\frac1m\sum_i\frac{\lambda_i}{w_i}
\longrightarrow2.
\]

The arithmetic threshold limit follows from the uniform bounds

\[
1-\frac1{w_*}\le u_i
\le1+\frac\Lambda{\theta v_0}.
\]

Together, `w_*,v₀→∞` imply `R_ar→1−θ=1/2`. Thus the claimed quality-score upper bound of `0.3` versus the arithmetic threshold `0.5` is justified under the explicit tuple hypotheses. Increasing a fixed collision coefficient does not change that particular limit, because `c/v₀→0`. Radius choices that change the collision coefficient with dimension require a new check.

The lower bound on the scalar product follows from the prefactor being at least one. It gives no lower bound on the accompanying test determinant, no dominance statement inside the sum, and no guarantee against further cancellation. These distinctions are correctly retained in the reviewed notes.

## 4. Completed tail note

The completed tail-cancellation note agrees with the finite identity and denominator-scaling proof in `REVIEW.md`. Its unchanged-collision combined error

\[
E_{\rm an}^{0}
=100K/w_0+\ln2/v_0+[\ln(800K^2)+\nu]/w_*
\]

has the correct constants: the `w_*^{-1}` contribution is the sum `ln2+ln(2K)+ν+ln(200K)`. The all-tail marker grouping remains valid after the actual-period expansion. The note correctly distinguishes individual tail terms from their zero sums, and it does not substitute an arbitrary formal period into the exponential periodicity identity.

## 5. Optimizing the Cauchy radius within this majorant

The additional radius-optimization calculation is correct. Define positive numbers

\[
d=\frac1{2v_0},\qquad A=\frac{7K}{w_0}.
\]

With `x=R₀/(7K)>1`, the available collision gain after paying only the necessary radial growth is

\[
F(R_0)=\frac{c_R}{v_0}-\frac{R_0}{w_0}
=d\ln x-Ax.
\]

Its derivative is `d/x−A`. Therefore

\[
\sup_{R_0>7K}F(R_0)=
\begin{cases}
d(\ln(d/A)-1),&d>A,\\
-A,&d\le A.
\end{cases}
\]

For `d≤A` the supremum is approached at the excluded endpoint `x=1`; this does not affect its use as an upper bound.

The quantifiers matter. For each fixed set of parameters and each `H`, choose an integral shell matching with quality-linear cost no greater than its expectation. This choice can be made **without reference to the Cauchy radius**. It is independent and has the deterministic occupancy bound for every radius. Consequently, for every `R₀>7K` fixed before taking `H→∞`,

\[
\limsup_{H\to\infty}
\bigl(\mathcal C_H(R_0)-E_{\rm hol}(R_0)\bigr)
\le\bar\nu J+F(R_0)
\le\bar\nu J+\sup F.
\]

Here replacing the full holomorphic maximum by `R₀/w₀` and dropping its positive `ln2/v₀` cost only weakens the upper bound. There is a single matching sequence for all fixed radii, although how large `H` must be can depend on the chosen radius.

On the earlier explicit family,

\[
\frac dA=\frac1{28K^2\theta^m},
\qquad C^2\theta=\frac{169}{200}<1.
\]

Since `K=floor(Cᵐ)`, `ln(d/A)=−m ln(C²θ)+O(1)=O(m)`, whereas `v₀` grows proportionally to `(13/12)ᵐ`. Hence `sup F=O(m/v₀)→0`. For mean quality tending to two, the score after the analytic cost still has upper limit at most `0.3`, below the arithmetic threshold tending to `0.5`.

This excludes rescue by choosing dimension-dependent radii `R₀=R₀(m)` **fixed before each H limit**, within the displayed termwise score and holomorphic majorant. It is not a statement about the size of the actual determinant. The argument does not interchange `sup_R` with `limsup_H`, does not cover an unrestricted `R₀(H)` by the same remainder argument, and does not exclude a genuinely stronger function estimate or another regrouping of the constant-shift sum.
