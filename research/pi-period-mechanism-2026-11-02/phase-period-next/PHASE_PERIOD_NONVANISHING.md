# Same-column actual-period nonvanishing for finite phase-dilated packets

2026-11-02. This note closes the rank check for sixteen explicit finite members of the phase-dilation construction. The complete source, the phase target, and one common set of genuine raw columns are retained. Both source and target minors are nonzero at the actual period. This does **not** give a theorem for all sufficiently large heights, a high-dimensional separated-weight construction, an actual determinant magnitude estimate, or a new rational approximation of pi.

## 1. Exact algebraic certificate criterion

For a finite row set E of labels (j,t,a), and raw columns Y^h X^alpha, define the formal polynomial matrix

\[
B_E(Z)_{j,t,a;h,\alpha}
=[z^tu^a]e^{hz}\prod_{i=1}^m(jZ+z+u_i)^{\alpha_i}.
\]

Each entry belongs to Q[Z]. Fix the row and column sets independently of Z. For a square choice C put P_(E,C)(Z)=det B_(E,C)(Z).

**Certificate criterion.** Suppose one common column set C has nonzero source and target determinants after evaluation at a Gaussian rational number r. Then P_(L,C) and P_(J,C) are nonzero rational polynomials. Since omega=2*pi*i is transcendental, both are nonzero at omega as well. No numerical approximation to pi or to polynomial roots is needed.

It suffices to certify the rational evaluation by a nonzero residue under any valid map from its denominator-localized Gaussian-integer ring into a finite field. This is a one-sided certificate: a zero residue would be inconclusive, but a nonzero residue proves the original determinant is nonzero.

The verifier uses r=44i/7 and p=1,000,033. It proves primality by trial division through floor(sqrt(p)) and explicitly finds a field element I with I^2=-1. The map i->I, 1/7->7^(-1) is valid. Rows are multiplied by t! product_i a_i! before evaluation, so the entries have no factorial denominators. The displayed parameter ranges also have all factorial indices below p.

Writing d=|alpha-a|, the scaled entry is zero when some alpha_i<a_i, and otherwise is

\[
\left(\prod_i\alpha_i^{\underline{a_i}}\right)
\sum_{k=0}^{\min(t,d)}\binom tk h^{t-k}
 d^{\underline k}(jZ)^{d-k}.
\tag{1}
\]

Formula (1) follows by extracting u first and differentiating e^(hz)(jZ+z)^d. The rising complexity of a numerical pi evaluation plays no role in the certificate.

## 2. The fixed model and common raw columns

The complete source is

\[
L_H=\{(j,s,\beta):0\le j<2,\ s+v\cdot\beta<H\},
\qquad v_i=w_i/\theta.
\]

All ambient columns satisfy w0*h+w dot alpha<=H. At N=s+|beta| use integer dilation 2 and the binary phase

\[
r_i(N)=\lfloor N/2^{i-1}\rfloor\bmod2.
\]

Put N0=ceil(2*theta*H/min_i w_i)+m. At N<N0 preserve the source row. Otherwise set a_i=2*beta_i+r_i(N) and t=N-|a|. This is exactly the preceding phase-cycled construction, including its low-order correction.

The verifier lists all source rows and allowed raw columns by deterministic integer loops. Gaussian elimination of the source evaluation at 44i/7 modulo p selects a raw column basis C. It then computes **both** square determinant residues on that very same C. It neither replaces raw columns by another family nor certifies the target on unrelated columns. Every selected column is recorded in the JSON output.

Both parameter families satisfy exactly

\[
v_0=1=2K\theta^m w_0,\qquad K=2,
\]

so their asymptotic source-to-column density is 1/2. They also satisfy the phase construction's strict time-budget hypothesis

\[
w_0>[m+2](m+1)v_0
\]

for L=2. Their finite source-to-column fractions need not equal 1/2, and that asymptotic relation is not substituted for an exact finite count.

| Family | m | theta | w0 | Spatial weights w | Source spatial weights v |
|---|---:|---:|---:|---|---|
| A | 1 | 1/28 | 7 | (1) | (28) |
| B | 2 | 1/16 | 64 | (1/2,4) | (8,64) |

These are finite rational model weights. In particular they are not claimed to equal ceil(ln q_i) for the actual approximation denominators, nor to satisfy the original high-dimensional interpolation theorem's separation hypotheses.

## 3. Sixteen complete finite certificates

| Family | Heights H | Complete source / target sizes M |
|---|---|---|
| A | 20, 28, 30, 32, 36, 40, 48, 56, 64 | 40, 56, 64, 72, 88, 104, 136, 168, 216 |
| B | 12, 16, 20, 24, 28, 32, 40 | 32, 48, 72, 96, 128, 160, 240 |

Every case has two nonzero determinant residues on one common raw basis. The criterion in Section 1 therefore proves actual-period nonvanishing for both minors in every listed case.

The checks also verify:

- the target has M distinct rows, every target time is nonnegative, and N is conserved;
- every coordinate moves upward and every target satisfies its transverse column-universe budget;
- all intrinsic raw-column time cuts for C, using d_r=sum_alpha(max(k_alpha-r,0));
- nonnegative coordinate degree surplus for C relative to both the target and source;
- every changed bucket's exact Pascal determinant.

For dilation 2, the last check has the particularly simple formula

\[
\det[\binom{2\kappa+r}\beta]_{\beta,\kappa\in B_N}
=2^{\sum_{\beta\in B_N}|\beta|}>0.
\]

It is checked by independent integer Bareiss elimination, not inferred from entrywise support. Thus at every nonzero shear error parameter the entire transport minor is nonzero. Combined with actual-period target nonvanishing, each listed case supplies a genuine nonzero term of the unreduced Cauchy-Binet expansion. No sign or size bound for the complete sum follows.

## 4. Why the earlier zero-center rigidity does not contradict this

The earlier anchor result concerned the original logarithmic-tail expansion. Its filtration inequality v0*T_i>v_i forced a nontrivial target to move to strictly smaller weighted order inside the same finite zero-center packet. Injectivity then forced identity transport on the whole zero-center block.

The current shear expansion has t=s-|a-beta| and a>=beta. With v_i>v0, its weighted target order is

\[
v_0t+v\cdot a
=v_0s+v\cdot\beta+\sum_i(v_i-v_0)(a_i-\beta_i),
\]

which increases. The old closure premise fails, so that particular anchor-rigidity proof cannot be applied. The concrete common-column certificates show directly that nontrivial zero-center shear movement need not kill the actual-period target minor. Different expansions of the same determinant can have different individual terms and different cancellations.

## 5. What remains open

The finite checks establish that the phase construction's support and transport mechanism is not automatically destroyed by actual-period rank, even on a common raw basis. They do not establish persistence as H grows with fixed parameters. They also do not cover m->infinity with strongly separated log-denominator weights and a target mean above the endpoint collision threshold.

A common nonzero finite minor can still be too expensive arithmetically or too large analytically. Even a successful infinite rank theorem would leave its quantitative size and complete-sum cancellation to be estimated. The scalar endpoint obstruction from the preceding round is unchanged.

Run `python verify_phase_period.py` in this directory. It uses only the Python standard library and writes `phase_period_verification.json` containing the exact selected raw columns, modular determinant residues, cut counts, degree surpluses, and escape counts. The modular certificate and transcendence argument are exact; no floating-point rank decision is used.
