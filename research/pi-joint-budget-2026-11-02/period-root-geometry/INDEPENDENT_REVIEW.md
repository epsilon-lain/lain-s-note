# Independent review of period-root geometry

The reflection identity, packet and budget checks, explicit root obstructions, and conditional relative magnitude estimate are correct. All eight explicit certificates and all fifty members of the six-row family were rerun successfully. An additional independent symbolic test covers transverse derivative rows and an odd reflection sign, which the fifty-case family does not exercise.

## 1. Reflection sign and the use of completeness

For a fixed \((j,\beta)\), multiplication of every column function by \(e^{Az}\) transforms the time rows by

\[
\text{new row}_{s,\beta}
=\sum_{k=0}^s\frac{A^{s-k}}{(s-k)!}
\text{old row}_{k,\beta}.
\]

Every row on the right is available in the complete weighted packet. This row transformation is triangular with diagonal one, so it does not change the determinant. It reduces the reflection \(h\mapsto A-h\) to the formal reflection \(h\mapsto-h\).

Applying \((z,u)\mapsto(-z,-u)\) to the original entry at center \(-Z\) contributes \((-1)^{s+\beta}\) from coefficient extraction and \((-1)^\alpha\) from the spatial polynomial. Therefore

\[
P_{C^{\rm ref}}(Z)
=(-1)^{\sum_L(s+\beta)+\sum_C\alpha}P_C(-Z),
\]

with the reflected column kept in its original position. The sign is correct. The intermediate negative exponents are only formal exponential jets; legality of the final nonnegative reflected columns is a separate budget check and is made in the note.

An independent test used \(K=2\), \(H=3\), \(v_0=1\), \(v_1=2\), \(w_0=1/4\), \(w_1=1/2\), retaining every row with \(s+2\beta<3\). There are eight rows, including two with \(\beta=1\). With increasing-time row order, take columns

\[
(h,\alpha)=(0,0),(0,1),(0,2),(0,4),(1,0),(1,1),(2,0),(2,1).
\]

Their reflection at \(A=2\) also obeys the budget. Direct SymPy coefficient extraction gives

\[
P_C(Z)=2Z^7,\qquad P_{C^{\rm ref}}(Z)=2Z^7,
\]

while the predicted sign is \(-1\). Thus the test verifies the genuinely signed relation \(P_{C^{\rm ref}}(Z)=-P_C(-Z)\), including nonzero transverse derivative orders.

## 2. Complete packets and exact root examples

For each displayed packet, \(v_1\ge H\) excludes all \(\beta\ge1\), and the strict inequality \(v_0s<H\) gives exactly the stated time levels. No permitted source rows have been dropped. The script reconstructs these packets from the inequalities and checks each distinct raw column against \(h+2\alpha\le H\).

The six-row family formula follows from the stated row subtraction. The upper \(\alpha=0\) block has ordinary-coefficient determinant one. The remaining determinant is \(Z^3\) times

\[
(Z+2\partial_a)\frac{(l-k)(a-k)(a-l)}2,
\]

which expands to the formula in the note. The fifty enumerated choices fit the genuine eight-height packet and their reflected budgets are valid.

The explicit factors correctly establish the geometric claims. In the first twelve-row case, the squared roots have opposite signs, producing nonzero roots on all four coordinate half-rays. A closed half-plane through the origin containing a pair of opposite real roots would have to have the real axis as its boundary; it cannot also contain both nonzero imaginary roots. Thus no such half-plane, or sector of width at most \(\pi\), contains that root set.

For the off-axis example, the negative discriminant gives two nonreal conjugate squared roots. Taking both square roots of each produces roots in all four open quadrants. For the remaining symmetric split, the stated positive discriminant is smaller than \(335^2\), so both squared roots are strictly negative. All these conclusions use exact signs and integer comparisons, not numerical root locations.

The eight-row reciprocal polynomial is correctly \(1+3x^2\). Its positive coefficients do not prevent zeros on the imaginary axis. This refutes that global shortcut without showing smallness at the actual reciprocal period.

## 3. Product ratio and the conditional logarithmic bounds

Factor the nonzero source polynomial as

\[
P_C(Z)=aZ^\ell\prod_{j=1}^{n-\ell}(Z-\rho_j),\qquad n\le D.
\]

Substitution into the exact reciprocal identity gives

\[
\mathcal F_C(\zeta)
=a\omega^\ell(1-\zeta)^{D-n}
\prod_j(\omega-\rho_j+\rho_j\zeta).
\]

Dividing by \(\mathcal F_C(0)=P_C(\omega)\ne0\) proves the displayed ratio formula. The power \(D-n\) is nonnegative. The zero roots of \(P_C\) have already canceled from this relative expression and are correctly omitted from the separation sum.

For any complex \(u\) with \(|u|\le1/2\), the triangle inequality gives

\[
-2|u|\le\log|1+u|\le |u|.
\]

Apply this to \(u=-\zeta\), repeated \(D-n\) times, and to \(u=\rho_j\zeta/(\omega-\rho_j)\) for every nonzero root. The two smallness hypotheses in the note ensure every factor is nonzero and the sum of these inequalities is exactly its conditional bound with \(\mathcal S_C\).

The bound is relative to \(|P_C(\omega)|\). It does not produce a baseline lower bound, and the hypotheses do not supply a uniform estimate of \(\mathcal S_C\). Both missing ingredients remain quantitative, even though transcendence ensures the denominator itself is nonzero.

## 4. Scope and verification

The standard-library verifier passed all eight explicit polynomial certificates, the fifty six-row identities and reflections, and exact nonvanishing at \(44i/7\). Its degree-bounded rational interpolation is exact, and its supplementary sample points are consistency checks rather than a numerical inference about roots.

The displayed nonzero roots all have modulus below five. For example the largest squared real root in the four-ray case is less than \((57+88)/8<25\); the off-axis roots have fourth power of their modulus equal to \(1604/9<625\). The remaining displayed factors satisfy the same bound immediately. Since \(|2\pi i|>6\), these finite examples are uniformly separated from the actual period by more than one.

**Audit outcome:** the examples rule out the stated universal root-location shortcuts. They do not construct roots approaching \(2\pi i\), small actual-period determinants, or an obstruction to a more restricted column-family theorem. The conditional ratio bound is valid with its two smallness hypotheses and supplies no unconditional magnitude conclusion.
