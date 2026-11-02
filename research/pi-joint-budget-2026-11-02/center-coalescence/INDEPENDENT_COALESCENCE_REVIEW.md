# Independent review of the center-coalescence bound

2026-11-02. Reviewed `CENTER_COALESCENCE.md`, concentrating on the weighted row change, time differentiation, arbitrary mixed columns, and the limiting normalization. The coalescence theorem and its stated necessary weight-ratio consequence are valid under the stated assumptions. No proof of a period-determinant magnitude estimate or a Diophantine conclusion follows.

This review uses \(V_H\) for the forced valuation and \(\mathcal V_H\) for the weighted polynomial column space, as in the final main note.

## 1. The maximal-weight anchor preserves exactly the complete row packet

Fix a time s and write R=H−v₀s. The spatial derivative indices form the complete lower packet

\[
\mathcal D_R=\{\beta\ge0:\sum_i v_i\beta_i<R\}.
\]

With v_m=maxᵢv_i and ξ_m≠0, set

\[
T=X_m/\xi_m,\qquad U_i=X_i-\xi_iX_m/\xi_m\quad(i<m).
\]

The derivative changes, in both directions, are explicitly

\[
\partial_{X_i}=\partial_{U_i}\ (i<m),\qquad
\partial_{X_m}=\xi_m^{-1}\partial_T
-\sum_{i<m}(\xi_i/\xi_m)\partial_{U_i},
\]

\[
\partial_{U_i}=\partial_{X_i},\qquad
\partial_T=\sum_i\xi_i\partial_{X_i}.
\]

Expanding an old derivative monomial replaces some derivatives of cost v_m by ones of cost at most v_m. Thus its new expansion stays inside the packet with weights (v₁,…,v_{m−1},v_m) on (U₁,…,U_{m−1},T). Expanding the inverse transformation has the same property. The two finite-dimensional derivative spaces therefore coincide, and the change between their indexed bases is invertible.

Factorial normalization of spatial jets only conjugates this change by invertible diagonal matrices. The change is independent of Z, is performed separately at each time s and each center, and does not mix in a missing low-time row. Consequently the determinant's order in Z is unchanged.

The maximal-weight condition is essential to this argument. For example, with weights (5,2), direction (1,1), and cutoff R=3, anchoring the second coordinate would send the permitted derivative ∂X₂ to ∂T−∂U₁. The latter summand has weight 5 and is absent. One cannot simply choose the cheapest anchor to obtain a stronger denominator in the final bound.

This is a statement about the row packet. The spatial coordinate change need not preserve the weighted column budget, and the proof does not require it to do so.

## 2. The time operator does commute with tangent translation

For the entries used in the note,

\[
[z^su^\beta]P(e^z,jZ\xi+z\mathbf1+u),
\]

time differentiation is the autonomous operator

\[
D=Y\partial_Y+\sum_i\partial_{X_i}.
\]

In the new coordinates this is

\[
\boxed{D=Y\partial_Y+\xi_m^{-1}\partial_T
+\sum_{i<m}(1-\xi_i/\xi_m)\partial_{U_i}.}
\]

All coefficients of the spatial derivatives are constant; Y∂Y is independent of T and U. Hence D commutes with ∂T, every ∂U_i, and translation in T. The normalized row is D^s∂X^β/(s!β!) evaluated at the center. After the row change, fixing (s,ν) leaves a polynomial in T whose consecutive normalized derivatives are evaluated at T=jZ,U=0,Y=1.

The time direction need not agree with ξ. Its transverse normal components are present in the displayed operator and cause no commutation problem. No alignment of those directions was silently assumed.

## 3. The wedge calculation does not assume independence of other blocks

For a fixed (s,ν), the tangent orders are ℓ=0,…,n−1, with

\[
n=\max\left\{0,\left\lceil\frac{H-v_0s-\sum_{i<m}v_i\nu_i}{v_m}\right\rceil\right\}.
\]

Write the associated row-vector function as f(T)=∑_q c_q T^q. Its normalized ℓ-jet at jZ is

\[
\sum_{q\ge\ell}\binom q\ell(jZ)^{q-\ell}c_q.
\]

In the wedge of the Kn rows in this block, any repeated Taylor index q repeats a coefficient vector and gives zero. Distinct q have total at least Kn(Kn−1)/2; subtracting the total derivative order K n(n−1)/2 gives the divisibility

\[
Z^{\binom K2n^2}.
\]

This is a divisibility assertion for the whole block wedge in the exterior algebra of the fixed column space. It remains true when its coefficient vectors satisfy relations. Wedge multiplication across different (s,ν) blocks therefore multiplies the powers of Z. Cross-block dependence may increase the order or annihilate the determinant, but cannot lower it. Thus

\[
\operatorname{ord}_0P_C(Z\xi)\ge
V_H:=\binom K2\sum_{s,\nu}n_{s,\nu}^2.
\]

For finite polynomial columns the T expansions are finite. There is no convergence or interchange-of-infinite-sums issue.

## 4. The degree argument is valid for arbitrary fixed mixed columns

Columns must be fixed independently of the formal center variable Z. A column may be any polynomial in the original allowed space, not necessarily a raw monomial or a polynomial homogeneous in transverse degree.

Let b_c be its largest weighted transverse degree. Extracting u^β from a monomial X^α gives zero unless α≥β coordinatewise. Every resulting center monomial r^γ has γ≤α−β, even after the time coefficient is extracted. Consequently each entry has weighted center degree at most b_c−w·β. A determinant term uses each column and each row once, proving

\[
\deg_wP_C(r)\le E_H:=\sum_c b_c-\sum_{L_H}w\cdot\beta
\le M_HH-\sum_{L_H}w\cdot\beta.
\]

If the radial polynomial P_C(Zξ) is nonzero and has order at least \(V_H\), some monomial of the multivariate polynomial contributing to a surviving radial coefficient has ordinary total center degree at least \(V_H\). Its weighted degree is at least \(w_{\min}V_H\). Hence

\[
w_{\min}V_H\le E_H.
\]

This reasoning does not confuse weighted degree with the degree after radial specialization. Division by w_min is precisely what connects them. It also does not require all lower-degree multivariate monomials to survive the chosen specialization.

For columns with individual transverse degrees α, the additional exact bound deg P_C(Zξ)≤∑ᵢD_i follows in the same way, with D_i=∑_Cα_i−∑_Lβ_i.

## 5. The lattice normalization and necessary inequality check out

The sums indexed by (s,ν) have m variables. At fixed m, K and positive weights, the elementary simplex integral

\[
\int_{u_i\ge0,\ \sum_{i=1}^m u_i<H}
(H-\sum_i u_i)^p\,du
=\frac{p!H^{m+p}}{(m+p)!}
\]

with p=1,2 gives

\[
M_H\sim\frac{KH^{m+1}}{(m+1)!v_0\prod_i v_i},
\qquad
\sum_{s,\nu}n_{s,\nu}^2
\sim\frac{2H^{m+2}}{(m+2)!v_0(\prod_{i<m}v_i)v_m^2}.
\]

Ceilings and strict cutoffs affect lower-order terms. Therefore

\[
\frac{V_H}{M_HH}\longrightarrow
\frac{K-1}{(m+2)v_m}.
\]

Also, because w_i=θv_i,

\[
\frac{\sum_{L_H}w\cdot\beta}{M_HH}
\longrightarrow\frac{m\theta}{m+2}.
\]

Combining the preceding bounds and v_m=w_max/θ yields exactly

\[
\frac{w_{\max}}{w_{\min}}
\ge\frac{\theta(K-1)}{m(1-\theta)+2}
\]

whenever the stated nonzero source minors exist for unbounded H. The limits are fixed-parameter statements. A later high-dimensional or large-weight limit requires the prescribed order of limits; no uniform simultaneous-limit error estimate has been proved here.

## 6. Sharpness and scope

For m=1, the complete Hermite staircase with columns (Y−1)^sX^α, α<K n_s, has determinant

\[
\prod_s\prod_{i<j}((j-i)Z)^{n_s^2}.
\]

It attains the claimed order exactly, with columns inside the original weighted space for sufficiently large H when w₀<v₀ and Kw₁<v₁. This confirms sharpness of the order theorem in that region; it does not supply a legal staircase in the critical high-dimensional regime.

The zero-time count Kθ^m≤1 is a separate necessary condition already used earlier. It is not a new consequence to be advertised as resolving the original interpolation assumptions. The coalescence inequality gives a different restriction, which can be compatible with strongly separated original weights.

The reviewed theorem controls a forced factor at center coalescence Z=0. The remaining quotient's magnitude at the actual period is still uncontrolled. Nonvanishing, coalescence order, and a quantitative actual-period lower bound remain distinct statements.

The accompanying `verify_center_coalescence.py` was also read and run independently. It passed all eight complete-packet polynomial checks, including two-variable directions (1,2) and (2,3), which do not coincide with the original time direction. Its two staircase equality cases and nine lattice-count checks passed. The finite polynomial degree bounds and row-factorial normalization used by that script agree with the proof reviewed above.
