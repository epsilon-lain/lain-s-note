# Exact aggregate period polynomials and controlled positive time surplus

2026-11-02. This note replaces a single selected transport term by the complete, unmarked transport sum on the **same** columns. It supplies exact finite cancellation and noncancellation certificates, a counterexample to universal real-rootedness, and an infinite-height family whose full aggregate can be calculated in closed form. None establishes a determinant magnitude bound, the critical high-dimensional endpoint, or any claim about π's bad approximability.

## 1. The complete sum has an exact reciprocal-polynomial form

Fix a complete source packet L and M polynomial columns C, either raw Y^hX^α or shifted (Y−1)^hX^α. All columns obey the original weighted budget. Put

\[
B^Z_{j,s,\beta;h,\alpha}
=[z^su^\beta]F_h(e^z)(jZ+z+u)^\alpha,
\quad P_C(Z)=\det B^Z_{L,C},
\quad D=\sum_C\alpha-\sum_L\beta,
\]

where F_h(Y)=Y^h or (Y−1)^h consistently. The full transport matrix replaces u by u−ζz. Its entry is

\[
[z^su^\beta]F_h(e^z)(j\omega+(1-\zeta)z+u)^\alpha.
\]

Factoring (1−ζ) from the second coordinate gives a column factor (1−ζ)^α, a row factor (1−ζ)^{−β}, and center ω/(1−ζ). Therefore the **complete sum**, including every target selection, is exactly

\[
\boxed{\mathcal F_C(\zeta)
=(1-\zeta)^D P_C\!\left(\frac{\omega}{1-\zeta}\right),
\qquad\omega=2\pi i.}
\tag{1}
\]

This is a polynomial identity after cancellation of the displayed denominators; no numerical substitution for ω is used. If P_C(Z)=∑_d p_d Z^d, then

\[
[\zeta^g]\mathcal F_C
=(-1)^g\sum_{d\le D-g}\binom{D-d}{g}p_d\omega^d.
\tag{2}
\]

In particular the aggregate degree is D−ord_{Z=0}P_C. A nonzero coefficient polynomial in ω does not vanish at the actual period, by the known transcendence of 2πi. It does not identify which particular target selection supplied a surviving contribution. There is also no outside-row marker in (1): the coefficient at error degree g sums selections with potentially different outside-row counts.

## 2. A complete six-row sum with exact cancellation

Take K=2, H=6, w₀=1,w₁=2,v₀=2,v₁=8. The complete source packet is

\[
L=\{(j,s,0):j=0,1,\ s=0,1,2\}.
\]

Use raw columns, in the indicated order,

\[
(h,\alpha)=(0,0),(1,0),(2,0),(0,1),(1,1),(2,2).
\]

Every column satisfies h+2α≤6. The exact ordinary-jet determinant is

\[
P_C(Z)=Z^3(Z+3),\qquad D=4.
\]

It is nonzero both at the genuine rational center r=44i/7 and at ω. Equation (1) gives

\[
\boxed{\mathcal F_C(\zeta)=\omega^3(\omega+3)-3\omega^3\zeta.}
\tag{3}
\]

For a row (j,s,0), select a∈{0,…,s} and target (j,s−a,a). There are exactly 36 complete target selections. Direct enumeration finds 11 nonzero target minors: one at Γ=0, four at Γ=1, and six at Γ=2. The six Γ=2 minors are individually nonzero and cancel exactly. Writing each center's chosen types as (a₀,a₁,a₂), their signed contributions are:

| Center 0 types | Center 1 types | Contribution to [ζ²] |
|---|---|---:|
| (0,0,0) | (0,0,2) | Z² |
| (0,0,0) | (0,1,1) | −Z² |
| (0,0,1) | (0,1,0) | Z² |
| (0,0,2) | (0,0,0) | −Z² |
| (0,1,0) | (0,0,1) | −Z² |
| (0,1,1) | (0,0,0) | Z² |

All target minors at Γ≥3 vanish individually in this example. Thus the Γ=2 disappearance is genuine cancellation, not a rank-zero explanation. The verifier checks the complete sum coefficient by coefficient against (3).

## 3. The genuine 32-row certificate has a nonzero complete degree-30 coefficient

Use K=2,H=24,w₀=1,w₁=2,v₀=2,v₁=16 and the complete source L={2s+16β<24}. The standalone fixture records the 32 selected columns from the previous genuine 22/7 certificate. The same index set works in raw and shifted bases.

The selected adjacent-interval transport packet has 22 outside rows and Γ=30. In both column bases its exact normalized period minor is

\[
\boxed{\det B^Z_{J,C}=-\frac{1697809}{147}Z^{32}.}
\tag{4}
\]

The **source** period polynomials differ between raw and shifted bases, but both satisfy

\[
D=62,\qquad \operatorname{ord}_0P_C=27,
\qquad\deg P_C=43.
\]

Therefore the complete aggregate has degree 35 in ζ. Its degree-30 coefficient is nonzero in both bases. For the raw basis the exact coefficient is

\[
\begin{aligned}
[\zeta^{30}]\mathcal F_C
={\omega^{27}\over13716864000}\big(&70079133260195\omega^5
+3287549503281792\omega^4\\
&+64369355710576640\omega^3
+646180039983434880\omega^2\\
&+3204437332426176000\omega
+6039420937570944000\big).
\end{aligned}
\tag{5}
\]

This is nonzero at ω, not merely at a formal-period sample. The full exact source polynomials and the shifted degree-30 coefficient are stored as rational coefficient lists in `aggregate_period_verification.json`. The verifier additionally evaluates each source polynomial exactly in Q(i) at r=44i/7 and confirms nonzero real/imaginary components.

The conclusion is about the complete coefficient at error degree 30. It does **not** say the individual selection in (4) can be identified as surviving cancellation, nor does it prove a nonzero coefficient with outside-row marker exactly 22. Those are different statements.

### A root-location shortcut fails on these genuine examples

Write P_C(Z)=Z²⁷R_C(Z). In both bases R_C has degree 16 and is squarefree. An exact rational Sturm calculation gives

\[
V(-\infty)=9,\qquad V(+\infty)=7.
\]

Hence R_C has exactly two real roots and fourteen nonreal roots. The reciprocal polynomial has the corresponding nonreal roots as well. Universal real-rootedness of these source or reciprocal period polynomials is therefore false even on the complete, genuine weighted packet and columns used here. The result does not rule out weaker sector or size bounds.

## 4. An infinite-height family with prescribed nonconstant aggregate

The preceding examples are finite. The following family works for every sufficiently large height in an explicit one-dimensional parameter regime.

Take positive fixed weights satisfying

\[
w_0<v_0,\qquad Kw_1<v_1,
\]

and the complete source L_H={v₀s+v₁β<H}. Write

\[
C_t=\left\lceil\frac{H-v_0t}{v_1}\right\rceil,
\qquad Q_{t,\alpha}=(Y-1)^tX^\alpha,
\qquad 0\le\alpha<KC_t.
\]

The original Hermite staircase has determinant

\[
P_0(Z)=\prod_{t}\prod_{0\le i<j<K}((j-i)Z)^{C_t^2},
\qquad D_0=\deg P_0,
\]

and zero time surplus. Here rows are ordered by time, then center and type, and columns by time and α; another consistent ordering only changes the common overall sign. At a flat adjacent pair C_{t−1}=C_t=B, replace the single column

\[
\boxed{Q_{t,KB-1}\quad\longmapsto\quad Q_{t-1,KB}.}
\tag{6}
\]

Keep the new column in the deleted column's position. Then

\[
\boxed{P_1(Z)=KB\,P_0(Z),\qquad D_1=D_0+1,
\qquad\Delta_{\mathrm{time}}=1.}
\tag{7}
\]

### Proof of the exact replacement coefficient

Set

\[
A_Z(X)=\prod_{j=0}^{K-1}(X-jZ)^B,
\quad X^{KB}=A_Z(X)+R_Z(X),\quad\deg R_Z<KB.
\]

The component (Y−1)^{t−1}R_Z is a combination of existing columns at level t−1. After subtracting it, the new column is (Y−1)^{t−1}A_Z(X). Its source rows of time below t vanish: for time t−1 every type β<B vanishes at all centers. At time t its data are the ordinary Hermite data of A'_Z, since the possible next coefficient of (e^z−1)^{t−1} multiplies A_Z and vanishes in those rows. The polynomial A'_Z has degree KB−1 and leading coefficient KB.

The time-t diagonal block of the original matrix is invertible Hermite interpolation. Therefore the new column's coordinate along the deleted column is exactly KB. Cramer's rule proves (7), for every nonzero Z. This proof works at the actual period directly.

### Multiple separated replacements

Choose times t₁<⋯<t_r with t_{i+1}−t_i≥2 and C_{t_i−1}=C_{t_i}=B_i. Perform (6) at each time. In coordinates of the original staircase, the new column belonging to t_j has zero coordinates at all times below t_j−1. Thus at every previously selected deleted column, whose time is at most t_j−2, its coordinate is zero. The submatrix of replacement coordinates at the deleted positions is triangular, with diagonal K B_i. Consequently

\[
P_r(Z)=\left(\prod_{i=1}^rKB_i\right)P_0(Z),
\quad D_r=D_0+r,\quad\Delta_{\mathrm{time}}=r.
\]

The entire aggregate is now explicit:

\[
\boxed{\mathcal F_r(\zeta)
=(1-\zeta)^r\left(\prod_{i=1}^rKB_i\right)P_0(\omega).}
\tag{8}
\]

Every coefficient at error degrees 0,…,r is nonzero. This is a complete determinant identity, so no unproved noncancellation assertion is needed.

### Budget, availability, and scale

Let η=max(w₀/v₀,Kw₁/v₁)<1. Each original column has weighted degree less than ηH+(K−1)w₁. Replacement (6) increases the degree by w₁−w₀, so the sufficient condition

\[
H>\frac{Kw_1}{1-\eta}
\]

puts all changed columns inside the original weighted space. They have integer coefficients. They are shifted polynomial columns; they are not generally a determinant-one change of a same-size raw-column subset. No raw arithmetic prefactor is claimed for them.

If v₁≥2v₀ and H>v₀, the last two time levels always have C=1. Hence at least one replacement is available for every sufficiently large H. Several separated flat levels give the stated multiple-replacement version.

For fixed parameters, r=O(H), whereas M_H=Θ(H²). With ζ fixed and ζ≠1,

\[
\frac{\log|(1-\zeta)^r|}{M_HH}\longrightarrow0.
\]

Thus this family certifies a genuine nonconstant aggregate, but does **not** produce the main-scale error gain required for near-total escape. Its exact value also shows why nonconstancy alone is insufficient for the proposed endpoint argument.

## 5. Standalone verification

Run `python verify_aggregate_period.py`. It uses only the Python standard library and the accompanying `genuine_32_fixture.json`; no other scratch folder is imported or read.

The verifier uses integer Bareiss determinants, rational forward-difference interpolation with a proved transverse degree bound, and two additional exact evaluation points. It exhausts all 36 small target selections and compares the full sum with (1). It also checks the 32-row source and selected-target polynomials, the genuine rational-center evaluation, and the rational Sturm chain. All coefficients are retained in the output JSON.

For the multiple-replacement theorem it additionally checks a complete 56-row packet with H=25,K=2,v₀=2,v₁=8,w₀=w₁=1 and replacement times 6 and 10. The exact results are

\[
P_0(Z)=Z^{72},\qquad P_2(Z)=8Z^{72},\qquad
D_2=74,\qquad\mathcal F_2(\zeta)=8\omega^{72}(1-\zeta)^2.
\]

These tests establish finite identities and check interfaces of the general proof. They do not replace the argument for arbitrary height, yield an actual-period magnitude bound, or certify a prescribed outside-row marker coefficient.
