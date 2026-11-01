# Independent review of the shear boundary and compatible-column route

2026-11-01. Reviewed `structured-gram-transfer/SHEAR_BOUNDARY_IDENTITY.md`,
its exact verifier, the proposed scalar boundary budget, and the compatible
column example supplied by the parent researcher. This review separates
proved identities from the still missing determinant estimate.

## 1. Shear identities: accepted with the displayed cost assumption

The substitution

\[
F(j\omega+t,j\delta+v)
=\widetilde F(j\omega+t,v-\lambda t)
\]

gives the transport coefficient
$\binom{\beta+\gamma}{\beta}(-\lambda)^\gamma$ without additional time
factorials. Transport composition is the binomial identity
$T(x)T(y)=T(x+y)$.

When $c_i=w_i/\theta-v_0\ge0$, transport cannot decrease the source row cost.
Thus the outside-to-inside block is zero, and the internal block really has
inverse $T(+\lambda)_{LL}$. Its determinant is one. The finite ambient set
$s+|\beta|\le N$ is valid because this quantity is preserved by transport.

The outside coefficient

\[
R_{\rho,\tau}
=\binom{\beta+\gamma}{\beta}(-\lambda)^\gamma
\sum_{\substack{\kappa\le\gamma\\c\cdot\kappa<d_\rho}}
(-1)^{|\kappa|}\binom\gamma\kappa
\]

has the correct sign, strict inequality and combinatorial factors. The
one-dimensional partial-binomial formula is correct, including the case of
an integral ratio $d_\rho/c$. If $c_i=0$ and $\gamma_i>0$, summation over
$\kappa_i$ cancels that coefficient exactly.

These statements require the cost hypothesis. They are not a claim that
arbitrary weighted packets or selected sparse row sets are shear invariant.

## 2. The analytic row bound: accepted, with endpoint wording clarified

The polydisc Cauchy estimate and the negative-binomial generating function
prove the displayed Chernoff bound. Dropping $|\gamma|\le s$ is an upper
bound, and discarding zero-cost directions is justified by exact cancellation.

Two wording corrections were sent to the author and are incorporated in the
reviewed final version:

1. The sheared function is precisely
   $P(e^z,(1+\lambda_i)z+v_i)$. A substitution
   $X_i\mapsto(1+\lambda_i)X_i$ alone would also multiply $v_i$ by
   $1+\lambda_i$ and is therefore not identical without an accompanying
   change of transverse coordinates.
2. The assumption $|\pi-p_i/q_i|\le q_i^{-2}$ alone does not say that the
   **largest** admissible Chernoff parameter is approximately $2\theta$.
   It guarantees an endpoint scale near $2\theta$ in the fixed-radius,
   large-weight regime; much smaller errors can permit larger parameters.
   An asymptotic equality requires the normalized error exponent to approach
   $2$ (and the other displayed cost/radius assumptions).

Neither issue changes the exact shear or boundary formulas. The final
polydisc bound must keep its column norm and beta-dependent product; its
isolated exponential factor is not already a determinant gain.

The added direct raw-column identities (7a)--(7b) also pass review. With
$t_i=r_i/\omega\ne0$ and $\zeta_i=1-t_i^{-1}$, factoring $t_i$ out of
$jr_i+z+u_i$ leaves $j\omega+z+u_i/t_i-\zeta_i z$. Extracting $u^\beta$
therefore gives exactly the stated row and column diagonal factors and the
transport $T(-\zeta)$. The determinant exponent is
$\sum_P\alpha_P-\sum_\rho\beta_\rho$. The restriction of that scalar
prefactor to raw square minors is essential and is correctly stated; mixing
different raw columns retains the diagonal matrix inside the mixing.

## 3. Scalar coefficient budget: valid and slightly sharpened

Assume first $c_i>0$ and $0<|\lambda_i|<1$. Put

\[
\tau=\min_i\frac{-\log|\lambda_i|}{c_i}>0.
\]

For an outside index with $w(\beta+\gamma)\le H$, one has
$c\cdot\gamma\ge d_\rho$ and
$|\lambda|^\gamma\le e^{-\tau d_\rho}$. Moreover

\[
\left|\sum_{\kappa\in S}(-1)^{|\kappa|}
\binom\gamma\kappa\right|\le2^{|\gamma|},
\qquad
\binom{\beta+\gamma}{\beta}2^{|\gamma|}
\le3^{|\beta+\gamma|}.
\]

The second inequality follows termwise from the multinomial expansion of
$(1+2)^{\beta_i+\gamma_i}$. Therefore the proposed factor $4^{H/w_*}$ can be
replaced by $3^{H/w_*}$:

\[
\boxed{|R_{\rho,\tau'}|
\le3^{H/w_*}e^{-\tau d_\rho},\qquad
\sum_{\tau'}|R_{\rho,\tau'}|
\le N_A3^{H/w_*}e^{-\tau d_\rho}.}
\]

The reuse of the letter tau for a target row should be avoided in exposition;
the exponent parameter above is scalar. Zero-cost directions can be omitted
because they cancel, and zero shift coordinates can be omitted because their
positive powers vanish. If no active directions remain, the remainder is
zero. These are coefficient bounds, not bounds for the outside analytic rows.

For the full fixed-parameter simplex, the mean row slack satisfies

\[
\frac{1}{MH}\sum_\rho d_\rho\longrightarrow\frac1{m+2}.
\]

This is the ordinary uniform-simplex slack mean in $m+1$ dimensions; the
center multiplicity cancels. If one forms the determinant term using the
boundary remainder in **every** row, multiplying these scalar estimates
gives the exponent contribution
$-\tau/(m+2)+\log3/w_*+o(1)$. A mixed determinant term receives only the gaps
of the rows actually replaced by boundary rows. No saving for all rows may
be assigned to such a mixed term. Norms of the outside rows and their mixed
minors remain to be bounded.

## 4. Compatible columns: the determinant really is $r^n$

Use columns $Y^h$ and $XY^h$, $0\le h<n$, and rows

\[
[t^s]P(1+t,jr+\log(1+t)),\qquad
j=0,1,\quad0\le s<n.
\]

Let $B_{s,h}=\binom hs$ and
$L_{s,h}=[t^s](1+t)^h\log(1+t)$. With the indicated block ordering,

\[
A_n(r)=\begin{pmatrix}B&L\\B&L+rB\end{pmatrix}.
\]

The matrix $B$ is upper triangular with diagonal one. Subtracting the first
block row from the second proves

\[
\boxed{\det A_n(r)=r^n.}
\]

The two transverse coefficient functions of the column $XY^h$ are
$f_0(z)=ze^{hz}$ and $f_1(z)=e^{hz}$; hence $f_0=z f_1$. The corresponding
$Y^h$ coefficient is also $e^{hz}$. They cannot be chosen independently in
the direct sum of Hardy spaces. At the genuine period, for fixed
$r=\omega+\delta\ne0$,

\[
\frac{|\det A_n(r)|}{|\det A_n(\omega)|}
=|1+\delta/\omega|^n,
\]

whose logarithm is only linear in $n$. This is consistent with a quadratic
lower growth rate for the generic evaluation Gram ratio: the latter bounds
all transverse coefficient choices and discards their compatibility.

The example uses only beta zero and is a kernel/column subproblem; it is not
the complete fixed-weight source simplex as $n$ grows. Rational checks of
the polynomial identity do not assert exponential periodicity at rational
surrogate nodes.

## 5. Verification and remaining proof obligation

`verify_review.py` uses a separate direct coefficient calculation for the
compatible columns and independently evaluates the weighted alternating
boundary sums. It passes 315 outside coefficient checks, 56 row-sum checks,
and 18 compatible-column determinant checks ($1\le n\le6$, including zero,
complex rational and imaginary rational values of $r$).

The final shear verifier was also rerun successfully: 1,478 boundary
coefficient checks, 5,535 polynomial substitution entries, 5,535 local
exponential substitution entries, 6 raw square determinant prefactors, and
100 positive-series majorant checks. These counts are finite verification
coverage, not a proof of the subsequent high-dimensional determinant bound.

The missing step is an asymptotic bound for the sum of the mixed determinants
of $B_L+RB_O$ which uses the shared coefficient functions, the boundary
coefficients and the arithmetic height simultaneously. The present results
do not establish such a bound, a contradiction at the endpoint, or bad
approximability of pi.
