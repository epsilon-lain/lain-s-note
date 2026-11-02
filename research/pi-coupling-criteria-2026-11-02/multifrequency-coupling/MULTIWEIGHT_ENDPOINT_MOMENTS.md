# More than two frequencies: an exact multiweight moment determinant

2026-11-02. This note changes the frequency structure while keeping a complete two-center source packet and legal raw monomial columns. It gives an all-size multiweight moment identity and a positive-Gram upper bound, together with small exact examples whose roots are not obtained merely by multiplying the earlier two-frequency polynomials. It does not prove a general primitive-content formula or a new rational approximation bound.

## 1. Complete packet and raw columns

Let n>=2. Choose distinct positive integer frequencies h_1,...,h_s and positive integers r_1,...,r_s with sum r_l=n. The 2n raw columns are

\[
X^a\quad(0\le a<n),\qquad
Y^{h_\ell}X^a\quad(1\le\ell\le s,\ 0\le a<r_\ell).
\]

Rows are (j,t,0), j=0,1 and 0<=t<n, ordered by center then time. Columns are ordered first by the zero-frequency block, then the displayed positive-frequency blocks, with increasing a within each block. A legal complete realization is

\[
H=\max(2n+\max h_\ell,5),\quad w_0=1,\quad w_1=2,
\quad v_0=H/n,\quad v_1=H,\quad\theta=2/H.
\]

The strict source inequality v_0*t+v_1*beta<H gives exactly beta=0 and t<n. All columns obey h+2a<=H. Moreover v_0>w_0, 2theta<1 and 2(w_0/v_0)theta<1. These parameters are not asserted to be height-optimal or to reach the original critical density regime.

Write

\[
\mathcal P_{n,\boldsymbol r}(Z)
=\det\bigl([z^t]e^{hz}(jZ+z)^a\bigr),
\]

over these rows and columns. Define the actual exponential endpoint determinant by replacing e^(hz) with e^(h(jZ+z)), and call it mathcal F. They agree at every integer exponential period 2*pi*i*k. They are different functions at general Z.

## 2. An exact multiweight moment identity

Let A_j=product_(a=0)^(j-1) a! and put

\[
E=n(n+1)/2+\sum_{\ell=1}^{s}r_\ell(r_\ell-1)/2,
\qquad C=\prod_{\ell=1}^{s}h_\ell^{nr_\ell}.
\]

Form an n by n matrix with rows i=0,...,n-1 and columns (ell,a), 0<=a<r_ell:

\[
M(Z)_{i,(\ell,a)}
=\int_0^1u^{i+a}e^{h_\ell Zu}\,du.
\]

Then

\[
\boxed{\mathcal F_{n,\boldsymbol r}(Z)
=\frac{C Z^E}{A_n^2}\det M(Z).}
\tag{1}
\]

To prove this, eliminate the zero-frequency n-column block, which is [I;T(Z)]. Each remaining entry becomes the order-n Taylor remainder of g(x)=e^(hx)x^a, differentiated t times and divided by t!:

\[
\frac1{t!(n-t-1)!}\int_0^Z(Z-x)^{n-t-1}g^{(n)}(x)\,dx.
\]

Within a frequency block, g^(n)(x)=e^(hx)(partial_x+h)^n x^a. The latter polynomial transformation preserves degrees below r_l and has determinant h_l^(n*r_l). Reversing the row polynomials (Z-x)^(n-t-1) into the monomial basis gives two cancelling reversal/sign factors and the denominator A_n^2. Substituting x=Zu contributes one Z per column, the sum of row degrees, and the sum of remaining column degrees. These are precisely the exponent E in (1). Thus there is no omitted sign or frequency multiplier.

There is also a formal version with independent endpoint symbols y_l, not only y_l=e^(h_l Z). Define

\[
\mu_{\ell,0}(Z,y_\ell)=\frac{y_\ell-1}{h_\ell Z},
\qquad
\mu_{\ell,k}(Z,y_\ell)=\frac{y_\ell}{h_\ell Z}
-\frac{k}{h_\ell Z}\mu_{\ell,k-1}(Z,y_\ell).
\]

Replace each moment in (1) by mu_(ell,i+a), and each raw endpoint column multiplier by y_l^j. The identity still holds. One may prove the single-frequency remainder-entry identity first as a polynomial identity in its endpoint symbol; it is affine in that symbol. The common row transformation and the separate frequency-block column transformations then preserve it for independent y_l. Setting every y_l=1 gives exactly mathcal P. Apparent negative powers cancel because the original raw determinant is a polynomial.

For s=1, h_1=1 and r_1=n, equation (1) recovers the previous two-frequency Hankel formula. For multiple frequencies, its columns carry genuinely different exponential weights.

## 3. Positive Gram normalization survives with multiple weights

Let J_m=[1/(i+j+1)]_(0<=i,j<m), the positive Gram matrix of monomials on [0,1]. Its determinant is det J_m=A_m^4/A_(2m). Normalize M(Z) on the left by J_n^(-1/2), and on the right by the block diagonal matrix with blocks J_(r_l)^(-1/2). For the resulting matrix G(Z), each column in frequency block h_l is the L^2[0,1] projection of e^(h_l Z u) times one real orthonormal polynomial of norm one. Its Euclidean norm is therefore at most

\[
\|G(Z)_{:,(\ell,a)}\|_2
\le e^{h_\ell\max(0,\operatorname{Re}Z)}.
\tag{2}
\]

To justify the projection statement, expand the normalized row basis in L^2[0,1]. Its coefficients against the weighted column polynomial are exactly the entries of that column. Bessel's inequality bounds the coefficient norm by the full function norm, and multiplication by e^(h_l Z u) has norm at most the displayed exponential. The basis polynomials are real, so the coefficient formula remains valid for complex Z.

Hadamard's determinant inequality now gives the sharper columnwise bound

\[
\boxed{|\det M(Z)|\le
\left(\det J_n\prod_{\ell=1}^{s}\det J_{r_\ell}\right)^{1/2}
\exp\left(\max(0,\operatorname{Re}Z)\sum_{\ell=1}^s h_\ell r_\ell\right).}
\tag{3}
\]

In particular the normalized determinant has magnitude at most one on the imaginary axis. There is no extra factor s^(n/2): estimating each column before taking the determinant is better than applying one common operator-norm bound to the concatenation.

At omega=2*pi*i this combines with (1) to give

\[
|\mathcal P_{n,\boldsymbol r}(\omega)|
\le\frac{C(2\pi)^E}{A_n^2}
\left(\det J_n\prod_\ell\det J_{r_\ell}\right)^{1/2}.
\tag{4}
\]

This is an upper bound only. The multiweight determinant need not have the earlier Vandermonde-square positive integral, and no uniform positive lower bound is claimed. In particular n=2 with profile r_1=r_2=1 and frequencies 1,2 has mathcal P identically zero. A new frequency by itself does not guarantee nonvanishing.

## 4. A finite root change that is not a product of old blocks

Take n=4, frequencies 1,2 and r_1=r_2=2. Exact computation gives

\[
\boxed{\mathcal P_{4,(2,2)}(Z)
=-\frac{Z^4(Z^2+15)}9.}
\tag{5}
\]

Its nonzero roots are plus or minus i*sqrt(15). Its exact coefficient denominator is 9, and the cleared polynomial is primitive with unit leading coefficient up to sign.

This polynomial cannot be a nonzero scalar times a product of the previous unscaled two-frequency polynomials P_j, j>=2. Its order at zero is four. Since ord_0 P_j=j, the only possible index partitions are 4 and 2+2; but P_4 has degree eight and P_2^2 has degree four, whereas (5) has degree six. This is a precise exclusion of that product mechanism, not a claim of impossibility of every algebraic transformation or rescaling.

Thus changing the frequency profile can change root locations even though merely adding normal-coordinate copies in the preceding block extension could not.

## 5. Small balanced comparisons do not show better roots near pi

For n=4,6,8 take r_1=r_2=n/2. These balanced examples also fit the earlier H=2n, v_0=2, v_1=2n, theta=1/n packet, so the comparison can use the same complete source and allowed weighted height. The companion exact certificates isolate a positive imaginary root in each of the following intervals; decimals here are exact rational endpoints:

| n | Imaginary part of the certified root | Least coefficient denominator | Primitive absolute value at 2pi i |
|---|---|---|---|
| 4 | (3.87298, 3.87299) | 9 | (38150,38151) |
| 6 | (6.26801, 6.26803) | 291600 | (67362199,67362200) |
| 8 | (6.28317604, 6.28317605) | 1666598976000000 | (55905687002,55905687003) |

Each selected root is farther from 2pi i than the previously certified near-period root of P_n with the same 2n source rows. That comparison is proved by disjoint exact rational root brackets and a Machin interval for pi. It does not assert that every complex root of the new polynomial is farther away.

On the other hand, all three primitive actual-period values are smaller than the corresponding earlier primitive values (approximately 279078, 5.2566*10^9, and 4.6368*10^13). The exact verification uses outward integer enclosures, not these rounded diagnostics. The values remain greater than one. Smaller primitive values and a nearer selected root are therefore distinct criteria; neither is a rational approximation theorem.

The complete n=6 and n=8 polynomials and their exact coefficients are in the verification output. The n=6 primitive leading coefficient is 2, so the unit-leading-coefficient rational-denominator theorem must not silently be applied to it for even denominators. No all-size minimal coefficient denominator, nonvanishing, root convergence, or optimal profile has been established here.

## Verification and remaining gap

The standard-library verifier checks seven complete legal packets and 63 independent raw/moment determinant identities, allowing independent endpoint symbols. It reconstructs the three balanced polynomials exactly, checks extra raw evaluations, their least denominators and primitive contents, isolates root branches with rational Sturm arithmetic, and compares actual-period values using rigorous pi intervals. The previous branch fixture retains its immutable source commit and source hash.

The structural identities (1)--(4) are proved for every displayed size and profile. What remains is to combine the new frequency dependence with a sufficiently strong general arithmetic-content theorem and a useful actual-period lower or local comparison estimate. Those steps are not supplied by the finite examples.


## Prior comparison sources

The original single-weight identity is in the [two-frequency Hankel analysis](../../pi-period-mechanism-2026-11-02/raw-root-mechanism/HANKEL_PERIOD_MAGNITUDE.md). The old root-branch data come from the [immutable certified coefficient and root records](../../pi-period-mechanism-2026-11-02/root-separation-next/raw_near_period_verification.json), with the fixture recording their exact source commit and SHA-256. The companion [coupled raw certificate](../coupled-packet-search/COUPLED_RAW_CERTIFICATE.md) concerns spatial-profile coupling instead; the two constructions should not be conflated.
