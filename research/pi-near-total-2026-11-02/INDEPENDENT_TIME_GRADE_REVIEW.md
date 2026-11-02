# Independent review of the same-column time-grade cancellation

This review checks the proposed finite-dimensional theorem independently. The result is valid with the midpoint diagonal weights stated below. It proves an exact cancellation for a specific family of columns; it does not establish a Diophantine conclusion for pi or extend the cancellation to the critical high-dimensional parameter regime.

## 1. Conventions and the support argument

Use ordinary coefficient jets (derivatives divided by their factorials). For a shifted column

\[
P_{h,\alpha}(Y,X)=(Y-1)^hX^\alpha,
\]

the actual-period entry at the row \((j,t,a)\) is

\[
B_{j,t,a;h,\alpha}
=[z^tu^a](e^z-1)^h\prod_i(j\omega+z+u_i)^{\alpha_i}.
\]

Consequently an entry is zero unless both \(h\le t\) and \(a_i\le\alpha_i\) for every coordinate. The same time support holds with rational centers. For a square period minor on a row set \(J\) and a fixed column set \(\mathcal C\), a nonzero determinant therefore requires

\[
\sum_{J}t\ge\sum_{\mathcal C}h,
\qquad
\sum_Ja_i\le\sum_{\mathcal C}\alpha_i.
\]

These necessary conditions follow from one nonzero permutation term in the determinant; cancellation among other terms cannot invalidate them.

For the complete shear transport from source rows \(L\), a nonzero entry preserves \(j\) and \(N=s+|\beta|=t+|a|\), and satisfies \(a_i\ge\beta_i\). If a square transport minor \(T_{L,J}\) is nonzero, select one nonzero permutation term and match its source and target rows. Then

\[
\Gamma_i=\sum_Ja_i-\sum_L\beta_i\ge0,
\qquad
\sum_Jt=\sum_Ls-|\Gamma|.
\]

In particular, the joint transport-period contribution can be nonzero only if

\[
|\Gamma|\le\Delta_{\rm time}
:=\sum_Ls-\sum_{\mathcal C}h,
\qquad
\Gamma_i\le D_i:=\sum_{\mathcal C}\alpha_i-\sum_L\beta_i.
\]

There is a useful simplification in the proof of strictness. If \(|\Gamma|=0\), all nonnegative increments in the selected permutation term are zero. Thus every matched pair has \(a=\beta\), then \(t=s\), and also the same \(j\). Hence the row sets coincide: \(J=L\). Therefore

\[
J\ne L,\quad \det T_{L,J}\ne0
\quad\Longrightarrow\quad |\Gamma|>0.
\]

This step needs no weighted-cost comparison and no separate positivity assumption on the cost increments. It follows directly from entrywise transport support.

If \(\Delta_{\rm time}=0\), every nonidentity **transport-supported** target has zero period minor on these same columns. This qualifier matters: the assertion is not that every unrelated target row set has zero period minor. Cauchy–Binet now leaves just the source term, with its unit triangular transport determinant.

## 2. Hermite staircase and the direction of the inverse

For one spatial variable, take

\[
d_s=\left\lceil\frac{H-v_0s}{v_1}\right\rceil,
\qquad 0\le s<H/v_0,
\]

with \(K\) centers and source rows \(0\le\beta<d_s\). At time \(h=s\), choose the \(Kd_s\) columns

\[
(Y-1)^sX^\alpha,\qquad 0\le\alpha<Kd_s.
\]

There are exactly as many columns as source rows at each time, so \(\Delta_{\rm time}=0\). Ordering source and column blocks by increasing time makes \(B\) block lower triangular. Its diagonal time block is

\[
\left(\binom{\alpha}{\beta}(jr)^{\alpha-\beta}\right)_{(j,\beta),\alpha},
\]

because \([z^s](e^z-1)^s=1\). This is an ordinary-coefficient confluent Vandermonde block. Its determinant, with center-major row order, is

\[
\prod_{0\le i<j<K}((j-i)r)^{d_s^2}.
\]

Thus \(B\) is invertible for \(r\ne0\). No additional factorials occur with the specified normalization. With raw derivative jets there would be nonzero row factors, which must be tracked consistently.

For completeness the total spatial degree surplus is

\[
D=\frac{K(K-1)}2\sum_s d_s^2,
\]

which is exactly the exponent of \(r\) in this determinant.

The inverse orientation is as follows. The entry of \(B^{-1}\) with column-label time \(h\) as its row index and source-row time \(s'\) as its column index vanishes when \(h<s'\). If \(C=B_O\) is the outside-period matrix, then \(C_{t,h}=0\) for \(h>t\). Consequently

\[
X:=CB^{-1},\qquad X_{t,s'}=0\text{ whenever }s'>t.
\]

Write the inside transport block as \(U=I+U_{<}\). Its off-diagonal entries move to a strictly smaller time. Its finite inverse has the same triangular support. The outside transport block also strictly lowers time, so

\[
R=U^{-1}T_{LO},\qquad R_{s,t}=0\text{ unless }t<s.
\]

It follows that

\[
N=RX\quad\text{satisfies}\quad N_{s,s'}=0\text{ unless }s'<s.
\]

Hence \(N\) is strictly block lower triangular and nilpotent. For every complex \(\eta\), and every shear parameter for which the boundary representation is defined,

\[
\det(B+\eta RC)=\det(I+\eta N)\det B=\det B.
\]

This exact identity concerns the boundary-normalized determinant. Any prefactor such as \(b^D\) in the original rational-jet identity remains present; it must not be discarded when returning to that determinant.

## 3. Column budget and arithmetic

Set \(v_1=w_1/\theta\). The chosen columns satisfy

\[
w_0s+w_1\alpha
<K\theta H+(w_0-K\theta v_0)s+w_1(K-1)
\le H\max(K\theta,w_0/v_0)+w_1(K-1).
\]

Thus \(K\theta<1\) and \(w_0<v_0\) give a sufficient budget condition for all

\[
H>\frac{w_1(K-1)}{1-\max(K\theta,w_0/v_0)}.
\]

The strict ceiling convention is correct: \(\beta<(H-v_0s)/v_1\) gives exactly \(d_s\) nonnegative integers, including when the quotient is an integer.

The shifted columns have an integral triangular expansion in the raw columns, with the same spatial exponent \(\alpha\) in every summand. Thus the established arithmetic denominator clearing still applies. In this staircase the expansion also stays inside the selected raw staircase: \(d_h\) is nonincreasing in \(h\), so lowering \(h\) does not remove an allowed \(\alpha\). The basis transformation has determinant one.

For a general selected raw-column family, one must check this closure before identifying its determinant with that of a shifted-column family. The full weighted column space is closed under lowering \(h\), but an arbitrary subset need not be.

## 4. Correct support-adapted Gram weights

For the established product-Gram bound, put

\[
\widetilde R=E^{1/2}RD^{-1/2},\qquad
\widetilde X=D^{1/2}XE^{-1/2}.
\]

The correct midpoint grading is

\[
E_{s,s}=\varepsilon^{4s},\qquad
D_{t,t}=\varepsilon^{4t+2},\qquad 0<\varepsilon<1.
\]

On the respective nonzero supports,

\[
\widetilde R_{s,t}
=\varepsilon^{2(s-t)-1}R_{s,t},\qquad
\widetilde X_{t,s'}
=\varepsilon^{2(t-s')+1}X_{t,s'}.
\]

Both exponents are at least one. Therefore both scaled matrices tend to zero, and

\[
\det(I+\widetilde R\widetilde R^*)
\det(I+\widetilde X^*\widetilde X)\longrightarrow1.
\]

The left-hand determinant squared is exactly \(|\det(I+RX)|^2=1\), so the infimum of the valid product-Gram bound over positive diagonal weights is exactly one.

The earlier candidate \(E_{s,s}=\varepsilon^{2s}\), \(D_{t,t}=\varepsilon^{2t}\) is insufficient: equality \(t=s'\) may leave a nonzero limit in \(\widetilde X\). The midpoint shift repairs that defect.

This limit is a fixed finite-matrix statement. It does not supply a bound uniform in \(H\), dimension, or shear size for a prescribed \(\varepsilon\). The exact determinant cancellation itself needs no limiting argument.

## 5. Why this does not close the critical regime

For a fixed-dimensional weighted source simplex, as \(H\to\infty\),

\[
\frac{\sum_Ls}{MH}\longrightarrow\frac1{(m+2)v_0}.
\]

Every admissible shifted column has \(h\le H/w_0\), so

\[
\frac{\Delta_{\rm time}}{MH}
\ge \frac1{(m+2)v_0}-\frac1{w_0}+o(1).
\]

If \(w_0>(m+2)v_0\), no selected column family can be time-saturated once \(H\) is sufficiently large. The asymptotic is taken with the dimension and weights fixed; applying it to a changing-dimensional family requires choosing sufficiently large \(H\) for each member or proving a uniform remainder separately.

This excludes this particular exact cancellation mechanism in that parameter range. It does not exclude cancellations with positive time surplus, other choices of determinant identity, or other Diophantine approaches.

## 6. Independent exact sanity check

An independent SymPy calculation used \(K=2\), source times \(s=0,1,2\), \(d_s=1\), and columns \((Y-1)^hX^\alpha\) with \(h=0,1,2\), \(\alpha=0,1\). The center difference \(r\), shear \(\xi\), and multiplier \(\eta\) were symbolic. The calculation checked the ordinary jets directly and returned

\[
\det B=r^3,\qquad N\ne0,\qquad N^3=0,\qquad
\det(I+\eta N)=1.
\]

It also checked every entry of the claimed supports of \(X\) and \(N\). Here \(r\) is a formal center difference in the period-normalized jet formula; no claim that an arbitrary real center has exponential period was used. The proof above, rather than this finite example, establishes the general statement.

**Audit outcome:** the theorem, staircase construction, inverse orientation, arithmetic transfer, and critical-regime limitation are correct with the qualifications above. The Gram statement requires the midpoint weights, which have been incorporated here.
