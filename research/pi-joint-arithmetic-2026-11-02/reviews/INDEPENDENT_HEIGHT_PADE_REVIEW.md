# Independent audit: reflection height, finite factors, and the Padé arithmetic benchmark

2026-11-02. Reviewed the supplied `arithmetic_route/REFLECTION_HEIGHT.md`, its factor fixtures and verification programs, and `literature/PADE_LINEAR_DEGREE_BENCHMARK.md`. The included reviews were not used as evidence. The input directory was not modified.

**Verdict:** the constant-term calculation, sharp leading height scale, finite factor/irreducibility certificates, and Padé arithmetic-cost calculation pass. One minor domain qualification should be added to Padé equation (13): its displayed deduction uses `b>=2`, equivalently `q>=3` for reduced `p/q`. The two remaining denominator cases can be treated separately. This qualification does not affect any asymptotic conclusion. No substantive arithmetic blocker was found.

This review accepts the earlier all-size normalization of the reflection factors, and the earlier nearest-root scale, as explicitly stated inputs. It independently checks the new constant/height proof and certifies the supplied finite factors from the defining determinants. It does not establish general irreducibility or verify the bibliographic/current-literature assertions in the literature section.

## 1. Reflection constant terms, including their signs

Set `n=2r`, `B_n=A_n/A_r^2`, and retain the note's definitions of `L_+` and `L_-`. At `Z=0`, the original row indexed by `s` vanishes precisely when `s` is odd for `L_+`, or even for `L_-`. There are exactly `r` such rows. Therefore the coefficient of `Z^r` is obtained by taking the linear term from each vanishing row and the constant term from the others. This is a determinant coefficient calculation, so there is no extra factor `r!`.

The entries used in that determinant are correctly stated:

\[
2\,2^{a-s}/(s-a)!,\qquad
a\,2^{a-s-1}/(s-a+1)!.
\]

Extracting `2^{-s}/s!` from row `s` and `2^a a!` from column `a` gives total extracted product one. The remaining row factors `2` and `1/2` also cancel. The resulting Pascal/shift determinant therefore has no missing normalization.

For the unit lower Pascal matrix `P` and upper shift `J`, the first `n-1` rows of `PJP^{-1}` are `-e_s+e_(s+1)` and its final diagonal entry is `n-1`. Deleting the coordinate pivots of the complementary parity creates no sign: the deleted row and column index sets coincide. The resulting determinants are exactly

\[
[Z^r]L_+=(-1)^{r-1}(n-1),\qquad
[Z^r]L_-=(-1)^r.
\]

Combining with the specified monic normalizations gives the signed, rather than merely absolute, formulas

\[
F_+(0)=(-1)^{r(r-1)/2+r-1}(n-1)B_n,
\]
\[
F_-(0)=(-1)^{r(r+1)/2+r}B_n.
\]

For `n=4,6,8,10,12`, the pairs of signs are respectively `(+,-),(-,-),(-,+),(+,+),(+,-)`. The stored coefficients agree with these signs and magnitudes.

## 2. Coefficient and Mahler heights

The entry bound `||entry||_1<=2(3/2)^a` follows by summing the binomial expansion in column `a` and bounding the remaining exponential coefficient by one. Determinant expansion gives `n!` products, each bounded by the product of its column bounds. Removing the exact monomial `Z^r` changes no coefficient norm. Consequently

\[
\|F_\pm\|_1\le B_n\,n!\,2^n(3/2)^{n(n-1)/2}
\]

is valid. Its extra logarithmic cost is `O(n^2)`, uniformly in the sign.

The lower bounds do not require assumptions about root distribution. Since `F_±` are monic,

\[
H(F_\pm)\ge|F_\pm(0)|\ge B_n,
\qquad M(F_\pm)\ge|F_\pm(0)|\ge B_n.
\]

The standard coefficient norm estimate `M(f)<=||f||_2<=||f||_1` supplies the required upper Mahler bound. Thus both claimed logarithms equal `log B_n+O(n^2)`.

Independently substituting `r=n/2` into the factorial-product asymptotic gives

\[
\log B_n=\frac{n^2}{4}\log n
+\left(\frac{\log2}{4}-\frac38\right)n^2+O(n\log n).
\]

Hence the leading term `(n^2/4)log n` is correct for both coefficient height and Mahler measure. The lower and upper bounds establish this all-size asymptotic; it is not inferred from the finite fixtures.

Conditional on the cited root distance `log delta_n=-2n log n+O(n)`, the ratio in section 3 is indeed `(8+o(1))/n`. It describes the height of the full reflection factor. The note correctly refrains from transferring that height to an unknown irreducible subfactor containing the chosen root.

## 3. Independent finite factor and irreducibility certification

A separate audit implementation, `finite-height/independent_finite_certificates.py`, uses standard-library integer and finite-field arithmetic. It does not call SymPy factorization or accept the stored `Q_irreducible` flags.

For each sign and `n=4,6,8,10,12`, it evaluates the defining reflection determinant after integral row scaling at `D+1` distinct arguments, where `D=n(n-1)/2` is an elementary upper degree bound for the raw determinant. It compares these evaluations with the signed normalized stored polynomial. The total of **330 exact Bareiss determinant evaluations** certifies complete polynomial identities by interpolation, not only sampled agreement: both compared polynomials have degree at most `D`.

For each stored modular factorization, the implementation checks the prime, the squarefree reduction, and the exact product identity. It verifies each modular factor using Rabin's irreducibility criterion with independently implemented finite-field operations. This certifies **22 modular factorizations and 71 modular irreducible factors**.

The degree-subset intersection argument is valid: a proper monic rational factor is an integer polynomial by Gauss's lemma, and reduction preserves its degree at every prime. Its degree must be a subset sum of the modular irreducible factor degrees. The recorded intersections exclude every proper degree for all ten supplied factors. Accordingly all ten are irreducible over `Q`.

The result is restricted to those ten factors. Nothing here proves irreducibility for later dimensions. The independently reconstructed coefficients and certificates are recorded in `finite-height/independent_finite_certificates.json`, with an additional detailed review in `finite-height/FINITE_CERTIFICATE_REVIEW.md`.

The supplied height and modular-irreducibility programs were also rerun successfully with their output redirected to a separate audit directory. In particular, the supplied Pascal checks through `n=40` pass. This rerun supplements the independent derivations above.

## 4. Padé normalization, value, root scale, and arithmetic cost

The exponential Padé remainder has the stated sign and factorial normalization. At `omega=2*pi*i`, symmetry of the beta density contributes the factor `-chi_m`. For odd `m`, the additional power of `i` gives

\[
C_m(\omega)=\frac{(2\pi)^{2m}m!}{2(2m+1)!}\chi_m>0.
\]

Stirling expansion therefore gives the stated linear term:

\[
\log|C_m(\omega)|=-m\log m+(1+2\log\pi)m+O(\log m).
\]

The coefficient ratio in section 4 is correct, so the maximum coefficient is exactly the constant coefficient, `(2m-1)!/(m-1)!`. Its height logarithm is

\[
m\log m+(2\log2-1)m+O(\log m).
\]

The integer coefficient proof is valid. At a reduced imaginary rational argument with `b>=2`, multiplying by `b^(m-1)` leaves the leading term as a unit modulo every prime dividing `b`. Because the degree is even, the value is rational rather than merely Gaussian rational. Its reduced denominator is exactly `b^(m-1)` and its absolute value is at least the reciprocal. Rewriting the polynomial in `Z^2` does not reduce this denominator exponent.

For the baseline criterion, multiplication by `b^(m-1)` and division of its logarithm by `m`, along `m/b->c`, gives

\[
\log(e\pi^2/c).
\]

The threshold `c=e*pi^2` is therefore correct. The explicit conclusion is for `c` strictly above or below that threshold; behavior at equality is not inferred from the displayed first-order asymptotic.

The derivative majorant is conservative but valid. For `b>=2`, subtracting the period value from the nonzero rational value gives a gap at least `1/(2b^(m-1))`. Since `|ia/b-omega|=2|p/q-pi|`, the denominator `4` in equation (13) is correct. Choosing odd `m~cb`, with `c>e*pi^2`, yields

\[
\log\left(\text{displayed lower bound}\right)
=-2cb\log b+O_c(b).
\]

Thus this particular baseline-plus-height argument is much weaker than a fixed multiple of `q^{-2}`. The note appropriately says that this is not an impossibility theorem for other uses of the Padé family.

**Minor qualification:** equation (13) should explicitly assume `b>=2` (or `q>=3`). Lemma (3) alone does not justify its nonzero-value step at denominator one. This does not alter the asymptotic argument; the cases `q=1,2` are separate finite-denominator cases.

## 5. Padé local geometry and verification limits

Although not needed for the height-cost calculation, the local coefficient error estimates and stated root scale were checked for consistency. The product deficit `k(k-1)/(2m)` controls both existing and omitted coefficients. Differentiation of its majorant gives the displayed derivative bound. Near `omega`, the limiting derivative is `-1`, so the displacement has the stated positive imaginary sign for odd `m`.

The conservative threshold `m>=10001` can be verified without numerical pi or exponential values. Use `e^4<81` and

\[
\cosh(1/2)-1<\frac{1/8}{1-1/48}=\frac6{47}.
\]

Then `6/47+972/10001<1/4` proves the derivative estimate, while `1296/10001<41/47` suffices for the unit-disk Rouché comparison. The remainder is bounded by `(8/(m+1))^(2m+1)<1/2`. These bounds justify the small-disk argument and the factors `4/5` and `4/3`. The resulting logarithmic root scale has the correct linear term `2+2 log(pi/2)`.

The optional Hermite–Biehler reduction is mathematically consistent conditional on the cited Bessel stability theorem: the transformation `B_m(Z)=Z^m y_m(2/Z)` preserves the left half-plane, and the odd part on the imaginary axis is a nonzero scalar multiple of `x C_m(ix)`. Source verification is outside this audit.

The independent script `independent_pade_checks.py` passed **1,333 exact remainder-series coefficient checks, 1,085 coefficient-majorant checks, and 1,890 full-denominator checks**, together with the rational threshold estimates above. The supplied numerical Padé verifier also passed. Numerical root-finding was not used as proof of uniqueness, simplicity, irreducibility, or the asymptotic statements.

The supplied literature's claims about recent papers and current irrationality-measure bounds were not adopted as mathematical inputs here and require the separate source review. Neither audited family currently proves bad approximability of pi.
