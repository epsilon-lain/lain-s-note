# Independent review of reciprocal resummation

The reciprocal-polynomial identity, its homogeneous-slice zero criterion, and the radial coefficient theorem are correct under the stated fixed algebraic-column assumptions. The distinction between exact nonvanishing and a quantitative estimate is essential and is retained in the note.

## 1. Polynomial support and valuations

For \(g_c(Y)X^{\alpha_c}\), the order of \(g_c(e^z)\) at \(z=0\) equals the order of \(g_c(Y)\) at \(Y=1\), because \(e^z-1=z+O(z^2)\). In each spatial factor, choosing \(d_i\) additional powers for \(z\) leaves center degree \(\alpha_{c,i}-\beta_i-d_i\). Thus every nonzero permutation term has exponent \(D_i-\delta_i\), with \(\delta_i\ge0\), and total additional spatial order at most \(\sum_Ls-\sum_c\nu_c\).

This proves the displayed support bounds. In particular \(D_i\ge0\) whenever the determinant is nonzero. The reversed expression is an ordinary polynomial, rather than a Laurent polynomial. Algebraicity of its coefficients follows from finite coefficient extraction from \(g_c(e^z)\), whose coefficients are algebraic multiples of rational factorial coefficients.

For raw monomials \(Y^h\), the individual valuation is zero. Improving that bound by combining columns requires a basis change inside the selected span. The positive-time-surplus note supplies the intrinsic filtration for this purpose. Its support consequences are valid without confusing the original raw subset with a same-index shifted subset; any determinant multiplier from a basis change must still be tracked if one uses its numerical value or arithmetic normalization.

## 2. Scaling direction and normalization

Put \(b_i=(1-\zeta_i)^{-1}\), \(r_i=b_i\omega\). Directly,

\[
jr_i+z+u_i=b_i\bigl(j\omega+b_i^{-1}z+b_i^{-1}u_i\bigr).
\]

Factoring \(b_i^{\alpha_{c,i}}\) from column \(c\) and taking the coefficient of \(u_i^{\beta_i}\) contributes \(b_i^{-\beta_i}\) from row \((j,s,\beta)\). Therefore

\[
P_C(b\omega)=b^D F_C(\zeta).
\]

The direction of the prefactor is correct. Reversing the center polynomial gives exactly

\[
F_C(\zeta)=\sum_\delta c_\delta
\omega^{|D|-|\delta|}(1-\zeta)^\delta.
\]

The initially excluded values \(\zeta_i=1\) are restored by polynomial identity. Factoring the unit-triangular inside shear block gives the normalized boundary determinant without any need to invert the source period matrix. Returning to the original center determinant retains \(b^D\).

## 3. Exact zero criterion and radial theorem

For fixed error multi-index \(\gamma\), group its coefficient by total degree \(d=|\delta|\):

\[
[\zeta^\gamma]F_C
=(-1)^{|\gamma|}\sum_d\omega^{|D|-d}
\frac{\partial^\gamma Q_d(1)}{\gamma!}.
\]

Different \(d\) give different nonnegative integer powers of \(\omega\). Since the other coefficients are algebraic, transcendence of \(\omega\) implies that this value is zero exactly when every homogeneous-slice coefficient is zero. A transcendental number cannot be a root of a nonzero polynomial over the algebraic numbers: otherwise it would be algebraic over an algebraic extension of the rationals, hence algebraic over the rationals.

Along the common radial line, the coefficient of \(t^k\) is

\[
(-1)^k\sum_{d\ge k}\binom dk Q_d(1)\omega^{|D|-d}.
\]

If \(d_*\) is the largest degree with \(Q_{d_*}(1)\ne0\), then the term with exponent \(|D|-d_*\) is nonzero for every \(0\le k\le d_*\). No other summand has that exponent. Every coefficient in this range therefore survives at the actual period. If all \(Q_d(1)\) vanish, the entire radial restriction is zero; this exceptional case is correctly separated.

The implication “zero first radial derivative implies a constant radial aggregate” is consequently valid, including the zero constant. It applies to the full signed aggregate, not to every exterior term separately. In multiple transverse dimensions it does not establish a nonzero derivative in every direction: homogeneous-slice directional derivatives can vanish even when another directional derivative survives.

The algebraic-coefficient hypothesis is material. If column coefficients were allowed to be arbitrary functions of the actual transcendental period, the distinct-power argument could fail. The note's fixed algebraic polynomial columns satisfy the hypothesis.

## 4. Finite checks and their scope

The supplied exact rational polynomial verifier was rerun; all fourteen cases passed, including genuine transverse derivative rows and cases whose common-radial restriction is identically zero.

An additional independent SymPy calculation used center-major source rows \((j,s,0,0)\), \(j=0,1\), \(s=0,1\), and columns

\[
\sqrt2,\quad X_1,\quad YX_1,\quad X_2.
\]

Direct coefficient extraction, rather than the supplied dictionary-polynomial implementation, gave

\[
P_C(r_1,r_2)=\sqrt2\,r_1(r_1-r_2),\qquad
F_C(\zeta)=\sqrt2\,\omega^2(\zeta_1-\zeta_2).
\]

This independently checks the scaling, sign, algebraic-coefficient extension, and radial/normal distinction. The source packet can be realized by \(K=2\), \(H=2\), \(v_0=1\), \(v_1=v_2=4\), with \(w_0=1/4\), \(w_1=w_2=1\), \(\theta=1/4\); the listed columns all fit the budget. Formal center variables were not treated as actual exponential periods.

This review checks the general reciprocal theorem and the small exact verifier. It does not independently reproduce the separate 32-row interpolation certificate cited in the aggregate companion work.

**Audit outcome:** the theorem passes. It removes ambiguity about exact complete-coefficient vanishing once the finite center polynomial is known. It does not provide a lower bound on a nonzero coefficient, control its growth with height or dimension, or establish that the aggregate at a particular tuple of approximation errors is nonzero.
