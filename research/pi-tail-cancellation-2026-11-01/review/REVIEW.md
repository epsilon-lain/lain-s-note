# Independent review: exact tail cancellation and denominator scaling

2026-11-01. This review checks the finite algebraic identities and revised bounds against the uploaded manuscript, §3. It does not independently verify Theorem 2.1 or prove a statement about the bad approximability of π.

**Verdict:** the exact tail-cancellation argument is valid. The weak cutoff condition, including equality, is sufficient. A further column-and-row denominator scaling removes the factor `m` from the proposed arithmetic error bound. Both conclusions require their stated row and column scope.

## 1. A finite unitriangular identity

Fix positive weights `v₀,w₁,…,wₘ`, a positive number `θ`, and `H>0`. The row packet at each center is

\[
I_H=\{(s,\beta)\in\mathbb Z_{\ge0}^{1+m}:
v_0s+w\beta/\theta<H\}.
\]

Let `L(t)=ln(1+t)` as a formal power series, and let `τᵢ(t)` have order at least the positive integer `Tᵢ`, with

\[
v_0T_i\ge w_i/\theta.
\]

For fixed polynomials `P` and center coordinates `cⱼᵢ`, form coefficient rows of

\[
P(1+t,c_{j1}+L(t)+u_1+z\tau_1(t),\ldots,
c_{jm}+L(t)+u_m+z\tau_m(t)).
\]

Here `z` is an independent formal variable; every extracted coefficient is polynomial in `z`. Write the resulting full-row matrix as `A_z`, and write `A_0` for the same matrix without tails. Exact expansion gives

\[
(A_z)_{j,s,\beta}
=\sum_{a\ge\beta}\binom a\beta
\sum_{k=0}^s z^{|a-\beta|}
[t^k]\tau(t)^{a-\beta}\,
(A_0)_{j,s-k,a}.
\]

An indicated term can be nonzero only when

\[
k\ge\sum_iT_i(a_i-\beta_i).
\]

Consequently its new row index satisfies

\[
v_0(s-k)+wa/\theta
\le v_0s+w\beta/\theta<H.
\]

Thus every row on the right remains in the same packet. The term `a=β,k=0` has coefficient one. Every other contributing term has `k≥1`, hence strictly smaller time index `s-k<s`. Ordering rows first by increasing `s` makes the change of rows lower unitriangular, block by block in the center label `j`:

\[
A_z=U_zA_0,\qquad \det U_z=1.
\]

This identity holds for every fixed column family, not only for monomials. For any square minor using **all packet rows** and a fixed set of columns,

\[
\det A_z=\det A_0.
\]

Neither full rank nor the interpolation theorem is needed for this identity; if both determinants are zero, it remains true.

## 2. Why every positive total tail degree cancels

Apply the constant translation from the rational centers to the actual logarithmic periods after the preceding expansion, or before it. These are translations of the same transverse variables, so their finite coefficient expansions agree. A summand of source formula (3.13), after replacing every tail by `zτᵢ`, contains the power `z^{|d|}`. A determinant term contains `z^{Σρ|dρ|}`.

The determinant is independent of `z`. Therefore, for each integer `D≥1`, the sum of all terms with total tail degree `Σρ|dρ|=D` vanishes exactly. In particular the whole determinant equals the sum using only `dρ=0` for every row. In that case the tail coefficient is `[t^k]1`, forcing `k=0`.

This is cancellation **between determinant terms**. It does not assert that each tail-containing term vanishes individually, nor that the terms without tails have no further cancellation. It applies at the true period using the source identity `e^{hjω}=1`; it does not impose that identity on an arbitrary formal replacement of `ω`.

## 3. Improved scalar estimate

Only the scalar

\[
\xi_{j,s,\beta,a}
=\binom a\beta\prod_i[j(r_i-\omega)]^{a_i-\beta_i}
\]

remains. Under the source assumption `|π−pᵢ/qᵢ|≤qᵢ^{−ν}`, and with `wᵢ=ceil(ln qᵢ)`,

\[
|j(r_i-\omega)|\le 2K e^\nu e^{-\nu w_i}.
\]

Since `binom(a,β)≤2^{|a|}`, `|a−β|≤|a|`, and `|a|≤H/w_*`, this gives

\[
\log|\xi|
\le-\nu w(a-\beta)
+H\frac{\ln2+\ln(2K)+\nu}{w_*},
\quad w_*=\min_iw_i.
\]

Zero scalars satisfy the corresponding modulus inequality. The `ν/F₀` and scalar `ln2/v₀` costs disappear, and `ln4` becomes `ln2`. The source holomorphic determinant estimate is still applicable with `ℓ=s`; its own `ln2/v₀` term is not removed by this argument.

## 4. Stronger arithmetic estimate by two-sided denominator scaling

The minimal admissible truncation is

\[
T_i=\left\lceil\frac{w_i}{\theta v_0}\right\rceil,
\quad G_i(t)=\sum_{1\le k<T_i}(-1)^{k+1}t^k/k,
\quad L_i=\operatorname{lcm}(1,\ldots,T_i-1),
\]

where `Lᵢ=1` for `Tᵢ=1`. The equality case in the cutoff is valid because the original row inequality is strict and the transformed weight is no greater.

For a monomial column `YʰXᵅ`, scale the column by `∏ᵢ(qᵢLᵢ)^{αᵢ}` and a row `(j,s,β)` by `∏ᵢ(qᵢLᵢ)^{−βᵢ}`. An entry is zero unless `α≥β`; otherwise it becomes

\[
\binom\alpha\beta[t^s](1+t)^h
\prod_i[q_iL_i(jr_i+G_i(t))]^{\alpha_i-\beta_i}.
\]

This is a Gaussian integer: `qᵢLᵢjrᵢ=2ijpᵢLᵢ`, and `LᵢGᵢ` has integer coefficients. Therefore a nonzero fixed monomial minor obeys the exact inequality

\[
\log|\Delta_H|
\ge-\sum_{\text{columns }\alpha}\sum_i\alpha_i\ln(q_iL_i)
+\sum_{\text{rows }\beta}\sum_i\beta_i\ln(q_iL_i).
\]

There is no entrywise common denominator factor `D_H`. With `Λ=4ln2`, the elementary lcm bound gives

\[
\ln L_i\le\Lambda(T_i-1)
<\frac{\Lambda w_i}{\theta v_0}.
\]

The column sum is consequently at most `MH(1+Λ/(θv₀))`. The row sum is at least `MH(b̄−θ/w_*)`, after simply discarding its nonnegative lcm part. Thus

\[
\boxed{\frac{\ln|\Delta_H|}{MH}
\ge-(1-\bar b)-\frac{\Lambda}{\theta v_0}
-\frac{\theta}{w_*}.}
\]

This verifies the proposed derivation and improves the initially proposed `Λm/(θv₀)+θ/w_*` bound. The factor `m` was a consequence of the source's entrywise denominator clearing, not a necessary cost of these minors.

For a sharper asymptotic form define

\[
\kappa_i=\frac{\ln(q_iL_i)}{w_i},\qquad
\kappa_{\max}=\max_i\kappa_i.
\]

The exact normalized column cost is at most `κmax`. On the full packet the row simplex moments give, with parameters fixed before `H→∞`,

\[
\frac1{MH}\sum_\rho\sum_i\beta_i\ln(q_iL_i)
\longrightarrow\frac{\theta}{m+2}\sum_i\kappa_i.
\]

Hence the lower bound can use the asymptotic arithmetic threshold

\[
\kappa_{\max}-\frac{\theta}{m+2}\sum_i\kappa_i.
\]

This sharper form keeps the favorable row-denominator contribution.

### The same threshold for appropriately pre-scaled common mixed columns

The common-mixed-column extension is also valid with the following explicit scaling. Let `A,B` be the original and test matrices on all `N` admissible monomials, and put

\[
Q=\operatorname{diag}\left(\prod_i(q_iL_i)^{\alpha_i}\right).
\]

Assuming both matrices have full row rank, the previously proved bounded mixing lemma supplies an integer matrix `C∈{0,1,2}^{N×M}` for which both `AQC` and `BQC` are nonsingular. After the row scaling `∏ᵢ(qᵢLᵢ)^{−βᵢ}`, the entries of `AQC` are Gaussian integers. Therefore

\[
\frac{\ln|\det(AQC)|}{MH}\ge b_{\kappa,H},
\qquad
b_{\kappa,H}=\frac1{MH}\sum_\rho\sum_i\beta_i\ln(q_iL_i).
\]

On the analytic side each new column is a sum of at most `N` columns with coefficient modulus at most `2exp(κmax H)`. Consequently the entire-function norm bound acquires the factor `2Nexp(κmax H)`. The collision argument is unchanged, and its normalized exponent acquires `κmax+ln(2N)/H`. With fixed parameters `ln N/H→0`. The comparison threshold is exactly

\[
\kappa_{\max}-b_{\kappa,H}.
\]

Thus the improved estimate also applies to this explicitly pre-scaled common-column extension. The positive arithmetic row contribution is retained, while the column scaling is paid once on the analytic side. It must not be charged a second time on the arithmetic side. The determinant-invariance identity in §1 is valid for these mixed columns as well.

## 5. Two boundary witnesses

**Below the cutoff.** Take `m=K=1`, `v₀=w₀=w₁=1`, `θ=1/2`, `H=2`. The complete packet is `(s,β)=(0,0),(1,0)`. Use columns `1,X` and any constant center `c`. With the full logarithm the matrix is

\[
\begin{pmatrix}1&c\\0&1\end{pmatrix},
\]

whose determinant is one. Choosing `T=1` gives `G=0` and determinant zero. Here `v₀T=1<2=w₁/θ`. Thus a blanket invariance assertion without the cutoff condition is false. The condition is sufficient uniformly, not necessary for every individual column choice.

**An incomplete row packet.** With the same weights and `θ`, choose the valid boundary cutoff `T=2`, so `G=t`, and take `H=3`. Retain only the row `(s,β)=(2,0)` and column `X`. The logarithmic entry is `−1/2`, whereas the truncated entry is zero. The omitted row `(0,1)` belongs to the full packet and is precisely the row needed for the triangular operation. Thus the theorem does not apply to an arbitrary sparse subset of source rows.

## 6. Resulting scope

The finite identity is unconditional algebra. The use of a nonzero source minor still depends on a valid interpolation/rank input. The arithmetic and scalar improvements are genuine but do not establish a positive endpoint comparison gap. In particular they neither remove the independently constructed low-cost sparse selections nor prove that the remaining constant-shift terms fail to cancel. Those are separate questions.
