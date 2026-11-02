# Independent review of the multiweight endpoint identity

2026-11-02. Independently reviewed `MULTIWEIGHT_ENDPOINT_MOMENTS.md`, including the final columnwise Gram bound. The complete packet, raw-column budgets, exact frequency and power factors, independent endpoint symbols, Gram estimate, and finite comparison scope all check. No mathematical blocker was found. This is an independent agent review, not external peer review.

## 1. Complete source and legal columns

With H=max(2n+max h_l,5), v_0=H/n, v_1=H, and theta=2/H, the strict source inequality is precisely s/n+beta<1. It contains exactly beta=0 and 0<=s<n at each of the two centers. The zero-frequency polynomial block has n columns, and the positive-frequency blocks have total size sum r_l=n, so the raw matrix is square.

For a positive-frequency column, a<r_l<=n gives h_l+2a<=max h_l+2n-2<H. The other budget and strict theta inequalities also hold. In the three balanced n=4,6,8 examples, the stronger same-height claim is valid: the zero-frequency block has degree at most n-1, and the h=1,2 blocks have degree at most n/2-1, so all their columns also lie in the earlier H=2n budget. This permits the stated same-source, same-allowed-height comparison.

## 2. The determinant prefactor and exponent are correct

Eliminating the zero-frequency columns leaves the normalized Taylor-remainder entries

\[
\frac1{s!(n-s-1)!}\int_0^Z
(Z-x)^{n-s-1}e^{hx}(\partial_x+h)^n x^a\,dx.
\]

Within a frequency block, the polynomial operator (partial+h)^n is triangular on degrees 0,...,r_l-1 and has h^n on every diagonal entry. Its block determinant is therefore h_l^(n r_l). These transformations do not mix different frequency blocks and do not introduce an additional sign.

The row polynomials have degrees n-1,...,0. Reversing their order contributes (-1)^(n(n-1)/2), while their leading coefficients contribute the same sign. The signs cancel. The product of the row factorial denominators is A_n^2. After this change, the entries are integrals of x^(i+a)e^(h_l x).

Changing x to Zu gives Z^(i+a+1) per entry. Factoring row and column powers yields

\[
\sum_{i=0}^{n-1}i+n+
\sum_l\sum_{a=0}^{r_l-1}a
=\frac{n(n+1)}2+\sum_l\frac{r_l(r_l-1)}2.
\]

This proves the displayed prefactor C/A_n^2 and exponent E exactly. For a single positive frequency h=1 with r=n, they reduce to the previous Hankel formula.

The independent-endpoint-symbol statement is also justified correctly. Each remainder-entry identity is affine in the endpoint symbol for its own frequency. The identity can first be established for that single frequency, and then the row transform and separate block transforms can be applied with all endpoint symbols independent. One must not infer independent symbols merely from a global identity among e^(h_l Z), since integer frequencies have multiplicative relations. The note uses the valid entrywise argument instead.

## 3. The sharper Gram bound needs no frequency-count factor

Let psi_i be the real orthonormal row polynomials of degree below n, and let phi_(l,a) be the real orthonormal basis of the l-th column polynomial block. After the positive Gram normalizations, one column has entries

\[
\int_0^1\psi_i(u)e^{h_lZu}\phi_{l,a}(u)\,du.
\]

Because psi_i is real, these are exactly the complex L^2 projection coefficients of e^(h_lZu)phi_(l,a). Bessel's inequality bounds their squared sum by the squared norm of that function. Since phi_(l,a) has norm one, the column norm is at most exp(h_l max(0,Re Z)). This does not require orthogonality between different frequency blocks.

Hadamard's determinant inequality multiplies the individual column bounds. Undoing the row Gram normalization and the block diagonal column normalization therefore gives

\[
|\det M(Z)|\le
\left(\det J_n\prod_l\det J_{r_l}\right)^{1/2}
\exp\left(\max(0,\operatorname{Re}Z)\sum_lh_lr_l\right).
\]

In particular the normalized determinant is at most one on the imaginary axis. The previously considered common operator-norm estimate for the concatenated blocks would have introduced an unnecessary factor s^(n/2); the final columnwise argument correctly avoids it. This is an absolute upper bound and supplies no positive lower bound or nonvanishing theorem.

## 4. Finite algebra and the retained limitations

An independent symbolic computation of the original 8-by-8 raw matrix at n=4, frequencies 0,1,2 and positive profile (2,2), gives exactly

\[
-Z^4(Z^2+15)/9.
\]

The denominator is 9 and its nonzero roots are plus or minus i*sqrt(15). The stated exclusion of products of the old unscaled P_j is valid: order four permits only the index lists (4) and (2,2), whose degrees are eight and four rather than six. This does not exclude rescaling or every other algebraic representation, and the note does not claim that it does.

For n=2 and profile (1,1) at frequencies 1,2, the two time-one raw rows are literally identical when the endpoint symbols are set to one. The period polynomial is therefore identically zero. This verifies the note's warning that adding a frequency does not guarantee a nonzero period determinant.

The n=4,6,8 root comparisons concern selected positive imaginary branches only. The exact brackets put each new selected branch below the corresponding previously certified branch, and both lie below 2pi. They do not locate every complex root, and do not show that the new polynomial's nearest root must be farther away. The primitive-value comparisons are a different statement and are kept separate. The values are smaller in the three checked cases but still larger than one.

The n=6 primitive leading coefficient is 2. Therefore the earlier exact rational-value denominator theorem for unit-leading even polynomials cannot automatically be applied when the rational denominator is even. The final note explicitly retains this restriction and makes no general arithmetic-content or rational-approximation claim.

## Independent checks

Independent exact symbolic arithmetic checked the n=4 raw polynomial above and the n=2 repeated-row identity. Nine additional raw/moment determinant comparisons were checked with independent rational endpoint symbols, across profiles (n;r_h)=(2;1,1), (3;2,1), and (4;1,1,2), at three nonzero rational centers each. All passed. These review checks used exact SymPy arithmetic; the accompanying main verifier remains independently runnable with the Python standard library only.

The main verifier additionally records the complete-packet checks, exact root brackets, primitive-value enclosures, and squared rational versions of the final Gram upper bounds. The all-size identities and norm inequality follow from the proofs, not from those finite computations.
