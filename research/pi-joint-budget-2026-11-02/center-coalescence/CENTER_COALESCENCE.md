# Center coalescence forces a positive degree cost

2026-11-02. This note strengthens the elementary finite-difference bound on the degree needed to distinguish K centers. It applies to the complete source packet and any chosen polynomial columns. It gives a necessary weight-separation condition, not a sufficient interpolation theorem or a quantitative lower bound at the actual period.

## 1. Fixed source and a line of centers

Fix K>=2, m>=1, positive rational weights with v_i=w_i/theta and 0<theta<=1. The source rows are

\[
L_H=\{(j,s,\beta):0\le j<K,\ v_0s+v\cdot\beta<H\},
\qquad M_H=|L_H|.
\]

Columns belong to

\[
\mathcal V_H=\operatorname{span}\{Y^hX^\alpha:w_0h+w\cdot\alpha\le H\}.
\]

Fix a direction xi=(xi_1,...,xi_m) with xi_m!=0, and order the weights so v_m=max_i v_i. Fix M_H selected columns independently of Z. Let P_C(Z xi) be their square source determinant, with entries

\[
[z^su^\beta]P_c(e^z,jZ\xi+z\mathbf1+u).
\tag{1}
\]

The entries in (1) are formal-center polynomials in Z. At common actual-period centers take xi=1 and Z=2*pi*i. At rational displaced centers the same formal polynomial is the already defined normalized source matrix. No exponential periodicity is assumed at arbitrary Z.

For a nonzero polynomial P, ord_0 P denotes its least power with nonzero coefficient; the zero polynomial satisfies every lower bound on this order.

## 2. Separate tangent derivatives from normal derivatives

Use fixed coordinates

\[
T=X_m/\xi_m,\qquad U_i=X_i-(\xi_i/\xi_m)X_m\quad(i<m).
\]

The centers become T=jZ,U=0. This coordinate change depends on the fixed direction xi but not on Z. At each time s the complete normalized derivative packet is invariant: replacing an expensive X_m derivative by a U_i derivative lowers its weight, and the inverse tangent derivative uses only coordinates of weight at most v_m. Both derivative changes are invertible and have determinant independent of Z. We make no assertion that this greatest-weight anchor preserves the original column budget; only the source row space is being changed.

At a fixed time s and normal derivative multi-index nu in Z_(>=0)^(m-1), the permitted tangent orders are precisely l=0,...,n_(s,nu)-1, where

\[
n_{s,\nu}=\max\left\{0,
\left\lceil\frac{H-v_0s-\sum_{i<m}v_i\nu_i}{v_m}\right\rceil\right\}.
\tag{2}
\]

For each column, taking its fixed time and normal derivative first leaves a polynomial function f_(s,nu;c)(T). Time differentiation corresponds to the operator

\[
\mathcal D=Y\partial_Y+\xi_m^{-1}\partial_T+
\sum_{i<m}(1-\xi_i/\xi_m)\partial_{U_i}.
\]

It commutes with translation in T and with its derivatives. Thus the transformed block consists exactly of n consecutive normalized tangent jets of f at each of the K nodes jZ. Its other coefficients and the normal time components do not affect the following Taylor argument.

## 3. A confluent block has a forced order at zero

For a block of n tangent derivatives at K centers, expand

\[
\frac{f^{(l)}(jZ)}{l!}
=\sum_{q\ge l}\binom ql(jZ)^{q-l}[T^q]f,
\qquad 0\le l<n.
\]

In the wedge of its Kn rows, repeated Taylor indices q give zero. Distinct nonnegative q have sum at least Kn(Kn-1)/2, while the sum of the derivative orders is K n(n-1)/2. Hence every nonzero wedge term has Z-degree at least

\[
\binom K2 n^2.
\]

This lower bound does not require the Taylor coefficient vectors from different blocks to be independent. Further dependencies can only increase the order or make the entire determinant zero. Multiplying the block wedges gives

\[
\boxed{\operatorname{ord}_{Z=0}P_C(Z\xi)
\ge V_H:=\binom K2\sum_{s,\nu}n_{s,\nu}^{\,2}.}
\tag{3}
\]

The row coordinate change has nonzero determinant independent of Z, so it does not alter the order. For a fixed finite column family all Taylor expansions here are finite in T.

For m=1, the previously constructed shifted Hermite staircase has determinant

\[
\prod_s\prod_{0\le i<j<K}((j-i)Z)^{n_s^2},
\]

and therefore attains (3). That staircase obeys the original budget in its stated region w_0<v_0 and K w_1<v_1 at sufficiently large H. Sharpness in this region is not a claim that it fits the critical high-dimensional parameters.

## 4. Leading order and a necessary separation condition

At fixed dimension and weights, weighted-lattice integration yields

\[
M_H\sim\frac{K H^{m+1}}{(m+1)!v_0\prod_i v_i},
\]

\[
\sum_{s,\nu}n_{s,\nu}^2
\sim\frac{2H^{m+2}}
{(m+2)!v_0(\prod_{i<m}v_i)v_m^2}.
\]

The strict cutoff and ceiling in (2) change only lower-order terms. Consequently

\[
\boxed{\frac{V_H}{M_HH}\longrightarrow
\frac{K-1}{(m+2)v_m}.}
\tag{4}
\]

For columns of single transverse degree alpha, the ordinary center degree is at most sum_i D_i, with D_i=sum_C alpha_i-sum_L beta_i. Thus every nonzero source determinant satisfies sum_i D_i>=V_H. This strengthens the elementary order M_H(K-1)/2 from finite differences across each center label; the latter misses overlap of consecutive tangent derivatives.

A useful conclusion also holds for arbitrary mixed polynomial columns. If b_c is the maximum weighted transverse degree in column c, coefficient extraction gives

\[
\deg_w P_C(r_1,\ldots,r_m)
\le E_H:=\sum_c b_c-\sum_Lw\cdot\beta
\le M_HH-\sum_Lw\cdot\beta.
\]

Every monomial of total center degree at least V_H costs at least w_min V_H. If P_C(Z xi) is nonzero, (3) therefore forces

\[
w_{\min}V_H\le E_H.
\tag{5}
\]

The complete source has mean weighted transverse degree m theta/(m+2). If nonzero source determinants exist for unbounded H in the fixed direction xi, combine (4)--(5) with v_m=w_max/theta to obtain

\[
\boxed{\frac{w_{\max}}{w_{\min}}
\ge\frac{\theta(K-1)}{m(1-\theta)+2}.}
\tag{6}
\]

This condition is necessary only. It can be weak when the original weights are already strongly separated; it does not contradict an interpolation theorem with separation hypotheses. The more refined selected-column value E_H may be used in (5) when available.

## 5. Parameter eligibility still includes the zero-time rows

There is also a simpler, previously used necessary count condition. At s=0 all rows factor through P(1,X), whose spatial column dimension has leading term H^m/(m! product w_i). The source has K theta^m times that leading row count. Hence asymptotic full source rank requires

\[
K\theta^m\le1.
\tag{7}
\]

The original construction uses a strict version. This is separate from the total row-to-column ratio. For example, if theta_m=1-c/m+o(1/m) with fixed c, then (7) bounds K_m by e^(c+o(1)); one cannot let K_m grow arbitrarily while assuming the old nonzero-source hypothesis still applies.

Equations (3)--(6) retain derivative-overlap information absent from that simple count. They do not resolve quantitative period estimates: a nonzero determinant can have the forced factor Z^V times a polynomial whose size at 2*pi*i remains uncontrolled. Positive coalescence order at Z=0 is distinct from vanishing at the actual period.

## Verification

`verify_center_coalescence.py` constructs complete one- and two-transverse-dimensional packets, chooses genuine allowed raw and shifted columns, and recovers their exact center polynomials by rational interpolation with a proved degree bound. It verifies the forced zero coefficients, a staircase equality case, and the fixed-parameter lattice normalizations. Formal integer centers are only polynomial evaluation points; no arbitrary evaluation point is called an exponential period.
