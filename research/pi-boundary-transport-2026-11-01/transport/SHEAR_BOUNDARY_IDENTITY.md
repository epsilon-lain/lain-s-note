# Exact removal of the internal part of a shear

2026-11-01. This is an algebraic continuation of the grouped Hardy Gram calculation. It gives an exact boundary remainder and a quantitative analytic bound for its individual rows. It does not establish a useful relative determinant bound in the high-dimensional regime, and it does not prove a new approximation property of π.

## 1. Aligning the centers

Write

\[
F(z,u)=P(e^z,z+u_1,\ldots,z+u_m),\qquad
\omega=2\pi i,\quad r_i=\omega+\delta_i,\quad
\lambda_i=\delta_i/\omega.
\]

After the previously established exact removal of logarithmic tails, the source rows are coefficient tests of
\(F(j\omega+\ln(1+t),j\delta+v)\).
For each fixed \((j,\beta)\), replacing the complete consecutive logarithmic time jets by ordinary time jets is a triangular row operation with diagonal one. We therefore work with

\[
A_{j,s,\beta}=[t^s v^\beta]F(j\omega+t,j\delta+v).
\]

The following identities apply to each column separately. Define

\[
\widetilde F(z,v)=F(z,v+\lambda z)
=P(e^z,(1+\lambda_1)z+v_1,\ldots,(1+\lambda_m)z+v_m),
\]

\[
B_{j,\ell,a}=[t^\ell v^a]\widetilde F(j\omega+t,v).
\]

Since \(\lambda j\omega=j\delta\), substitution gives the exact finite identity

\[
\boxed{
A_{j,s,\beta}
=\sum_{\substack{\gamma\ge0\\|\gamma|\le s}}
\binom{\beta+\gamma}{\beta}(-\lambda)^\gamma
B_{j,s-|\gamma|,\beta+\gamma}.}
\tag{1}
\]

One way to verify (1) is to use
\(F(j\omega+t,j\delta+v)=\widetilde F(j\omega+t,v-\lambda t)\)
and expand each transverse monomial. In particular, there is no extra time factorial or multinomial coefficient in (1): those coefficients are already accounted for by the normalized Taylor tests.

The centers of every \(B\) row are now \((j\omega,0)\). The column functions have changed to the displayed \(P(e^z,(1+\lambda_i)z+v_i)\); only the longitudinal part of each transverse argument is scaled. These functions are analytic objects used to estimate the same source matrix; they are not new rational or Gaussian-integral columns. Arithmetic lower bounds must still refer to the original matrix.

## 2. The weighted packet and its boundary

Put \(b_i=w_i/\theta\), and let the original full row packet be

\[
L=\{(s,\beta):q(s,\beta):=v_0s+b\cdot\beta<H\}.
\]

Assume in this note that

\[
c_i=b_i-v_0\ge0\qquad(1\le i\le m).
\tag{2}
\]

Every transport step in (1) changes the weighted row cost by

\[
q(s-|\gamma|,\beta+\gamma)-q(s,\beta)=c\cdot\gamma\ge0.
\tag{3}
\]

Let \(N\) be at least the largest \(s+|\beta|\) in \(L\), and introduce the finite ambient index set
\(\Gamma_N=\{(s,\beta):s+|\beta|\le N\}\).
Transport preserves \(s+|\beta|\), so it never leaves this set. Write \(O=\Gamma_N\setminus L\).

For a parameter vector \(x\), let \(T(x)\) have entries

\[
T(x)_{(s,\beta),(s-|\gamma|,\beta+\gamma)}
=\binom{\beta+\gamma}{\beta}x^\gamma.
\]

All unspecified entries are zero. Direct binomial convolution, or composition of the substitutions \(v\mapsto v+xt\), gives
\(T(x)T(y)=T(x+y)\).
Ordering indices by the time order makes these matrices triangular with diagonal one. By (3), transport starting outside \(L\) never reenters it. Therefore

\[
T(-\lambda)=
\begin{pmatrix}U&V\\0&W\end{pmatrix},
\quad \det U=1,
\quad U^{-1}=T(\lambda)_{LL}.
\tag{4}
\]

This holds independently at every center \(j\). For a square matrix containing all source rows, left multiplication by these blocks preserves its determinant. We obtain

\[
\boxed{\widehat A:=U^{-1}A=B_L+RB_O,\qquad R=U^{-1}V.}
\tag{5}
\]

Thus every internal transport term has been removed exactly. The remaining terms use only rows outside the weighted packet. Equation (2) is essential to the proof; with both positive and negative cost increments, a transport path could leave the packet and return, and the displayed inverse block formula would need modification.

## 3. Explicit boundary coefficients

For a source row \(\rho=(s,\beta)\), define its positive distance in row weight to the boundary by
\(d_\rho=H-v_0s-b\cdot\beta\).
For a target \(\tau=(s-|\gamma|,\beta+\gamma)\) with \(|\gamma|\le s\) and \(c\cdot\gamma\ge d_\rho\), the coefficient in (5) is

\[
\boxed{
R_{\rho,\tau}
=\binom{\beta+\gamma}{\beta}(-\lambda)^\gamma
\sum_{\substack{0\le\kappa\le\gamma\\c\cdot\kappa<d_\rho}}
(-1)^{|\kappa|}\binom\gamma\kappa.}
\tag{6}
\]

Indeed, the intermediate row in the product \(T(\lambda)_{LL}T(-\lambda)_{LO}\) is indexed by \(\kappa\). Its membership in \(L\) is precisely the strict inequality in (6), and

\[
\binom{\beta+\kappa}{\beta}
\binom{\beta+\gamma}{\beta+\kappa}
=\binom{\beta+\gamma}{\beta}\binom\gamma\kappa.
\]

The sum in (6) is a finite, explicitly computable alternating binomial sum over a weighted simplex. It is not generally equal to zero or to a single monomial. If \(c_i=0\) and \(\gamma_i>0\), it is zero: summing over \(\kappa_i\) gives \((1-1)^{\gamma_i}\). Consequently directions of zero cost increase disappear entirely from the boundary remainder.

For one transverse coordinate with \(c>0\), set

\[
r_\rho=\left\lceil d_\rho/c\right\rceil-1.
\]

If \(r_\rho\ge s\), the row has no remainder. Otherwise all its remainder indices satisfy \(r_\rho<\gamma\le s\), and the partial alternating-binomial identity gives

\[
\boxed{
R_{\rho,(s-\gamma,\beta+\gamma)}
=\binom{\beta+\gamma}{\beta}(-\lambda)^\gamma
(-1)^{r_\rho}\binom{\gamma-1}{r_\rho}.}
\tag{7}
\]

This identifies the exact first possible error order of each row. A row far from the boundary requires more powers of the approximation error before an unresolved term can appear. Rows with \(s=0\) have no remainder at all, even if they are close to the boundary.

### A direct formula using the original period column functions

For raw monomial columns \(P=Y^hX^\alpha\), the changed functions can be eliminated by an explicit diagonal scaling. Assume \(r_i\ne0\), and put

\[
t_i=r_i/\omega=1+\lambda_i,
\qquad \zeta_i=1-t_i^{-1}=\frac{r_i-\omega}{r_i}.
\]

Let \(A^r\) denote the ordinary local jets
\([z^s u^\beta]e^{hz}\prod_i(jr_i+z+u_i)^{\alpha_i}\),
and let \(B^0\) denote the same expressions with \(r_i\) replaced by \(\omega\). The latter are precisely the original period coefficient rows, without any changed transverse column functions. Define diagonal matrices on source rows and raw columns by

\[
(D_{\rm row})_{\rho,\rho}=t^\beta,
\qquad(D_{\rm col})_{P,P}=t^\alpha.
\]

Then, as a rectangular matrix identity before selecting columns,

\[
\boxed{D_{\rm row}A^rD_{\rm col}^{-1}
=T(-\zeta)_{L,\Gamma_N}B^0.}
\tag{7a}
\]

Indeed, factor \(t_i^{\alpha_i}\) from each transverse monomial and use
\(jr_i+z+u_i=t_i(j\omega+z+u_i/t_i-\zeta_i z)\).
Extracting \(u^\beta\) contributes \(t^{-\beta}\), and expanding the remaining shift gives \(T(-\zeta)\). This is a formal local identity. It is valid even at surrogate nodes where no exponential periodicity is asserted; the interpretation as the original global period functions uses the actual \(\omega=2\pi i\).

Removing the internal block as in (4)–(6) gives, for any square selection of raw monomial columns,

\[
\boxed{
\det A^r
=t^{\sum_P\alpha_P-\sum_\rho\beta_\rho}
\det\bigl(B^0_L+R(\zeta)B^0_O\bigr).}
\tag{7b}
\]

Here \(R(\zeta)\) is exactly (6) with \(\lambda\) replaced by \(\zeta\), independently at each center. This expresses the entire unresolved effect through the weighted boundary and leaves the original compatibility between period coefficient rows intact.

For a rectangular family followed by a mixing matrix \(C\), (7a) remains valid with \(D_{\rm col}C\) on the right after undoing the displayed scaling. The scalar prefactor in (7b) is a raw-square-minor statement; it must not be extracted unchanged from a mixture of different raw columns.

## 4. A uniform analytic bound for each boundary row

Fix one column and one center. Suppose \(\widetilde F\) is holomorphic in a neighborhood of the closed polydisc

\[
|z-j\omega|\le h,\qquad |v_i|\le\rho_i,
\]

where \(h,\rho_i>0\), and has modulus at most \(D\) there. Put \(I=\{i:c_i>0\}\) and

\[
z_i=\frac{2h|\lambda_i|}{\rho_i}\qquad(i\in I).
\]

For any \(a\ge0\) satisfying \(z_i e^{ac_i}<1\) for all \(i\in I\), the remainder \(E_\rho=(RB_O)_\rho\) obeys

\[
\boxed{
h^s\rho^\beta |E_\rho|
\le D\,e^{-a d_\rho}
\prod_{i\in I}(1-z_i e^{ac_i})^{-\beta_i-1}.}
\tag{8}
\]

Here and below an empty product is one; if \(I\) is empty, the actual remainder is identically zero. To prove (8), Cauchy's coefficient bound gives
\(|B_{s-|\gamma|,\beta+\gamma}|\le D h^{-(s-|\gamma|)}\rho^{-(\beta+\gamma)}\).
The absolute value of the alternating sum in (6) is at most \(2^{|\gamma|}\), and it vanishes if a zero-cost direction occurs. Therefore the normalized remainder is bounded by

\[
D\sum_{\substack{\gamma\ge0\\c\cdot\gamma\ge d_\rho}}
\binom{\beta+\gamma}{\beta}z^\gamma.
\]

We have dropped the finite restriction \(|\gamma|\le s\), making an upper bound. Insert
\(1\le e^{a(c\cdot\gamma-d_\rho)}\), then sum the product of negative-binomial series. This proves (8). The factor \(e^{-a d_\rho}\) must be kept together with the displayed \(\beta\)-dependent product and the polydisc norm \(D\); it is not by itself a determinant saving.

For genuine approximations \(|\pi-p_i/q_i|\le q_i^{-2}\),
\(|\lambda_i|\le 1/(\pi q_i^2)\). If \(w_i\) is approximately \(\ln q_i\), \(b_i=w_i/\theta\), and \(v_0\) and the radii are fixed, this upper error bound suggests that a fixed \(a<2\theta\) becomes admissible for sufficiently large weights. Much smaller approximation errors may allow larger \(a\); no maximal admissible value is inferred from the upper bound alone. This is only a scale calculation. It does not choose the radii, bound \(D\), or control the accumulated \(\beta\)-dependent factor in the actual high-dimensional family.

The same estimate applies directly to the remainder in (7b), replacing \(\lambda\) by \(\zeta\) and taking \(D\) to bound the original period column function. In that version \(|\zeta_i|\le q_i^{-2}/|p_i/q_i|\). This formulation avoids an additional change in the period column functions.

## 5. What this does and does not transfer

The rows of \(B_L\) are complete ordinary jets at the true period, so their evaluation Gram matrix has the exact Cauchy product already established. Its column coefficient functions are now those of \(\widetilde F\), whose compatibility across transverse types is still known explicitly. Formula (5) is a way to retain that compatibility while isolating the failure of shear closure in the row packet.

A determinant expansion of \(B_L+RB_O\) contains terms with different numbers and positions of boundary rows. Bounding each entry by (8), or asserting that \(\lambda\) is tiny, does not show that these terms are negligible relative to \(\det B_L\). In particular, this note does not contradict the previously proved conditional relative instability in the critical high-dimensional regime. Any further determinant argument must control the mixed row minors, the column norm, and arithmetic height together.

The next concrete target is therefore the mixed exterior-product norm of the compatible rows in (5), with the exact coefficients (6), rather than a uniform condition-number comparison between two full Gram matrices.

## 6. Verification

`verify_shear_boundary.py` uses exact rational arithmetic. It checks transport composition, the invariant block decomposition, the boundary formula, its one-dimensional simplification and zero-cost cancellation, and the original analytic substitution independently on polynomial columns. It separately checks (7a) using finite Taylor coefficients of the local exponential columns, and verifies (7b) on square minors at several nonzero rational surrogate centers. It also checks finite positive majorants against an exact rational instance of (8), using integer costs and \(e^a=2\). These checks accompany the general proofs; they do not verify the manuscript's interpolation theorem or infer a statement about π from surrogate rational nodes.
