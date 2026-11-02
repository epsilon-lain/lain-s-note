# Exact root-geometry obstructions for genuine period minors

2026-11-02. The complete source packet and original raw weighted columns do not force Hurwitz stability, a zero-free imaginary axis, or containment of all nonzero roots in a half-plane through the origin. The counterexamples below are exact small determinants, not numerical root plots. A conditional product estimate identifies the additional root separation information a magnitude argument would actually need.

These results do not show that roots can approach 2πi arbitrarily closely, that actual-period determinants become arbitrarily small, or that a critical high-dimensional estimate fails. π's bad approximability remains unresolved by this investigation.

## 1. Setup and genuine packets

Use the coefficient-normalized period rows

\[
B^Z_{j,s,\beta;h,\alpha}
=[z^su^\beta]e^{hz}(jZ+z+u)^\alpha,
\qquad P_C(Z)=\det B^Z_{L,C}.
\]

In every example K=2, w₀=1,w₁=2, all source rows satisfying v₀s+v₁β<H are retained, and each selected raw column satisfies h+2α≤H. All source packets below happen to have β=0, because v₁≥H; this is an exact consequence of the finite weights, not a discarded portion of the packet.

| Packet | H | v₀ | v₁ | Complete source rows |
|---|---:|---:|---:|---|
| Six rows A | 6 | 2 | 8 | j=0,1; s=0,1,2; β=0 |
| Six rows B | 8 | 3 | 10 | j=0,1; s=0,1,2; β=0 |
| Eight rows | 8 | 2 | 10 | j=0,1; s=0,…,3; β=0 |
| Twelve rows | 12 | 2 | 14 | j=0,1; s=0,…,5; β=0 |

Rows are ordered by center then time. Columns are ordered by α and then by the listed h values. All determinants below use this convention. The verifier evaluates each polynomial exactly at the genuine rational center r=44i/7 and confirms nonzero. At the actual period ω=2πi the polynomials are also nonzero. None of their algebraic roots equals ω; for the displayed cases one may alternatively check the explicit roots have modulus below 5, while |ω|>6.

## 2. Positive and negative real roots already occur in six rows

In packet A, take h values

| α | First column family | Second column family |
|---:|---|---|
| 0 | 0,1,2 | 0,1,2 |
| 1 | 0,1 | 1,2 |
| 2 | 2 | 0 |

The exact polynomials are respectively

\[
\boxed{P_-(Z)=Z^3(Z+3),\qquad P_+(Z)=Z^3(Z-3).}
\tag{1}
\]

Thus neither universal left-half-plane nor universal right-half-plane root containment is possible. These are common genuine raw-column minors, not arbitrary polynomial combinations.

More generally, use α=0 frequencies 0,1,2; α=1 frequencies k,l; and α=2 frequency a. Direct calculation gives

\[
\boxed{P_{k,l,a}(Z)=\frac{l-k}{2}Z^3
\left((a-k)(a-l)Z+4a-2k-2l\right).}
\tag{2}
\]

To derive it, subtract each center-0 row from the corresponding center-1 row. The α=0 columns disappear in the lower block, while their upper three-by-three determinant is one. Factoring Z from each lower row leaves the Taylor-coefficient columns

\[
e^{kz},\quad e^{lz},\quad (Z+2z)e^{az}.
\]

Their determinant is (Z+2∂ₐ) applied to the ordinary Vandermonde determinant (l−k)(a−k)(a−l)/2, proving (2).

The verifier exhausts all 50 choices 0≤k<l≤4, 0≤a≤4 on packet B. Every chosen column obeys h+2α≤8. It verifies (2) and the reflection identity below for every case.

## 3. Frequency reflection structurally reverses the root geometry

Let the complete source packet be arbitrary in one transverse dimension, and let A be an integer at least every selected h. Reflect each raw frequency by h↦A−h, keeping its column position and α unchanged. Whenever both original and reflected column families obey the budget,

\[
\boxed{P_{C^{\mathrm{ref}}}(Z)
=(-1)^{\sum_L(s+\beta)+\sum_C\alpha}P_C(-Z).}
\tag{3}
\]

Multiplying every column function by e^{Az} induces a determinant-one triangular row transformation: at each fixed (j,β), it only uses smaller time orders, all of which belong to the complete source packet. It therefore suffices to reflect h to −h. Replacing z,u by −z,−u contributes the row factors (−1)^{s+β}, and replacing jZ−z−u by −(−jZ+z+u) contributes the column factors (−1)^α. This proves (3).

The intermediate formal identities do not assert that negative-frequency columns belong to the allowed polynomial space. Only the final nonnegative reflected columns are used, and their weighted budgets must be checked. The examples and the exhaustive finite family do check them. Thus a one-sided root principle would require an additional hypothesis not preserved by this legal reflection.

## 4. Imaginary-axis roots occur in genuine six- and eight-row minors

Packet B with α=0 frequencies 1,3; α=1 frequencies 0,4; and α=2 frequencies 1,3 gives

\[
P_C(Z)=-8Z^3(3Z^2+4),
\]

so its nonzero roots are ±2i/√3.

An especially simple monic example uses the eight-row packet:

| α | Frequencies h |
|---:|---|
| 0 | 0,1,2,3 |
| 1 | 0,3 |
| 2 | 1,2 |

Its polynomial and reciprocal polynomial are

\[
\boxed{P_C(Z)=Z^4(Z^2+3),\qquad
Q_C(x)=x^6P_C(1/x)=1+3x^2.}
\tag{4}
\]

Thus P_C has imaginary-axis roots ±i√3. Even the reciprocal polynomial with **positive coefficients** in (4) has imaginary-axis roots ±i/√3. Positivity of real coefficients cannot be used as a positive lower bound on the imaginary axis where the aggregate is evaluated.

This does not mean this particular polynomial is small at the actual argument 1/ω. Its zeros are a definite distance away. It rules out a global zero-free-axis shortcut.

Exchanging the frequency pairs to α=1:1,2 and α=2:0,3 instead gives P_C(Z)=Z⁴(Z²−13), with opposite real roots ±√13.

## 5. One complete twelve-row minor has roots on all four coordinate half-rays

On the twelve-row packet choose

| α | Frequencies h |
|---:|---|
| 0 | 0,1,2,3,4,5 |
| 1 | 1,4 |
| 2 | 0,2,3,5 |

The exact result is

\[
\boxed{P_C(Z)=\frac{Z^6}{4}(4Z^4-57Z^2-276).}
\tag{5}
\]

The squared nonzero roots are

\[
Z^2=\frac{57\pm\sqrt{7665}}8.
\]

Since 57²<7665, one is positive and the other negative. Hence this **single** allowed minor has nonzero roots on the positive and negative real axes and on the positive and negative imaginary axes. No closed half-plane whose boundary passes through the origin contains all its nonzero roots. In particular a universal sector of angular width at most π is impossible.

For an off-axis version, choose α=1 frequencies 2,3 and α=2 frequencies 0,1,4,5, retaining the same six α=0 frequencies. Then

\[
\boxed{P_C(Z)=\frac{Z^6}{9}(9Z^4-148Z^2+1604).}
\tag{6}
\]

The discriminant of the quadratic in Z² is 148²−4·9·1604=−35840. Its two squared roots are nonreal conjugates; their square roots occupy all four open quadrants. No numerical root approximation is needed for this conclusion.

For completeness the remaining symmetric frequency split, α=1:0,5 and α=2:1,2,3,4, gives

\[
P_C(Z)=\frac{Z^6}{36}(36Z^4+335Z^2+380),
\]

whose squared roots (−335±√57505)/72 are both negative. These three symmetric twelve-row splits are all verified exactly.

## 6. The magnitude question reduces to quantitative separation, not coefficient support

Let a nonzero genuine source polynomial factor over the complex numbers as

\[
P_C(Z)=aZ^\ell\prod_{j=1}^{n-\ell}(Z-\rho_j),
\qquad n=\deg P_C\le D,
\]

where nonzero roots are repeated with multiplicity. The exact aggregate formula implies

\[
\boxed{\frac{\mathcal F_C(\zeta)}{\mathcal F_C(0)}
=(1-\zeta)^{D-n}
\prod_j\left(1+\frac{\rho_j\zeta}{\omega-\rho_j}\right).}
\tag{7}
\]

The denominator is nonzero because P_C has rational coefficients and ω is transcendental. Define the root-separation statistic

\[
\mathcal S_C=(D-n)+\sum_j\frac{|\rho_j|}{|\omega-\rho_j|}.
\]

If |ζ|≤1/2 and |ζ|·|ρ_j|/|ω−ρ_j|≤1/2 for every root, then the elementary inequalities log(1−x)≥−2x and log(1+x)≤x for 0≤x≤1/2 give

\[
\boxed{-2|\zeta|\mathcal S_C
\le\log\left|\frac{\mathcal F_C(\zeta)}{\mathcal F_C(0)}\right|
\le|\zeta|\mathcal S_C.}
\tag{8}
\]

This is a conditional relative estimate. It neither proves the required separation nor supplies an absolute lower bound for |P_C(ω)|. It identifies two genuinely quantitative missing ingredients: a bound on the roots' weighted proximity to ω and a baseline determinant estimate. Knowing that each allowed Taylor coefficient is nonzero supplies neither ingredient.

The examples rule out several simple global geometric assumptions that might have supplied such control. They do not exclude a narrower root theorem for a specially chosen column family, nor do they construct roots converging to ω.

## 7. Reproducible exact verification

Run `python period-root-geometry/verify_period_root_geometry.py`. The script is self-contained and uses only the Python standard library. It verifies eight explicit genuine raw-column certificates and all 50 members of (2), including their frequency reflections.

Each determinant is computed by integer Bareiss elimination after factorial scaling of time rows. Exact rational interpolation uses the proved degree bound D, and two additional evaluation points are checked. The script then verifies the displayed polynomial coefficients, exact rational-center nonvanishing, and the integer sign/discriminant inequalities proving the root claims. Complete rows, columns, coefficients and rational evaluations are recorded in `period_root_geometry_verification.json`.

The fixed small heights and weights are not a successive-large-weight or high-dimensional asymptotic construction. No general determinant-size or Diophantine conclusion is claimed.
