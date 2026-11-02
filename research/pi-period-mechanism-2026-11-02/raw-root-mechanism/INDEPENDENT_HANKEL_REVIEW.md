# Independent review of the actual-period Hankel estimate

2026-11-02. Independently reviewed `HANKEL_PERIOD_MAGNITUDE.md`, read and ran `verify_hankel_magnitude.py`, and separately derived the determinant transformation and the uniform lower-bound polynomial.

**Result: passed.** The displayed two-sided estimate holds for every integer n>=2 in the specified complete one-dimensional two-frequency family. It gives actual-period nonvanishing and magnitude, not just finite numerical evidence. Its row-factorial normalization result has the stated scope; it is not a claim about the least possible cleared denominator or about the critical high-dimensional construction.

## 1. Exact determinant identity and sign

The distinction between the formal-period polynomial P_n(Z) and the actual endpoint determinant F_n(Z) is essential and is correctly maintained. Their equality at omega=2*pi*i follows from exp(omega)=1; it is not assumed away from that point.

With frequency-first columns and center-first rows, the actual endpoint matrix is

\[
\begin{pmatrix}I&E\\T&e^ZET\end{pmatrix}.
\]

Its determinant is det(e^ZET-TE), with no ordering sign. The entries after elimination are normalized derivatives of order-n Taylor remainders. The operator (1+partial)^n on polynomials of degree below n is triangular with diagonal entries one. Its inverse is therefore a legitimate determinant-one column change.

In the resulting moment matrix, the row polynomials (Z-t)^(n-1-s) have decreasing degrees. Reversing their degrees and multiplying their leading signs each contribute (-1)^(n(n-1)/2), so their product is +1. The row and column factorial normalizations give exactly A_n^(-2). Finally, the substitution t=Zu contributes Z^(n^2), proving the sign and constant in

\[
F_n(Z)=Z^{n^2}A_n^{-2}
\det\left[\int_0^1u^{i+j}e^{Zu}\,du\right]_{i,j<n}.
\]

I also obtained the same identity independently by writing the remainder of e^x x^a as

\[
\frac{x^n}{(n-1-a)!}\int_0^1(1-u)^{n-1-a}e^{ux}\,du
\]

and dividing the local jet function by x^n. The corresponding jet transformation is triangular with determinant Z^(-n^2); the reversed (1-u)-polynomial basis again has determinant +1. The derivation holds first for Z!=0 and extends to Z=0 by analyticity.

Retaining a separate symbol for exp(Z) before setting it to one is also justified. After clearing powers of Z, any failure of the formal identity would be a finite nonzero polynomial combination of distinct exp(kZ). Looking at the greatest exponent on the positive real axis rules that out. Thus the Laurent-moment formula for P_n is a valid algebraic identity.

## 2. Positive measure and exact moments

Expanding both Vandermonde determinants under the integral gives the factor 1/n! exactly. The Hilbert determinant equals A_n^4/B_n by the Cauchy determinant formula. Reflection u_i -> 1-u_i preserves the normalized Vandermonde-squared density and makes X=sum(u_i)-n/2 symmetric. Therefore the phase factor e^(int/2) and the real cosine expectation in the note are correct.

The Schur-moment formula follows by multiplying one alternant by a Schur polynomial, integrating with the other alternant, and evaluating the resulting Cauchy determinant ratio. Iterated one-cell addition gives the coefficient k!/product(hook length). For n>=6 every partition needed for moments through order six has length at most n. The verifier cross-multiplies the generic rational formulas as polynomial identities in n; this establishes the algebraic identities on the whole generic range rather than fitting finitely many values. For n=1,...,5 the length-restricted partition calculation avoids invalid zero-over-zero substitutions.

As an additional independent check, I formed the truncated moment-generating matrix

\[
\mu_{i+j}(t)=\sum_{k=0}^6\frac{t^k}{(i+j+k+1)k!}
\quad\text{in }\mathbb Q[t]/(t^7)
\]

and computed its determinant by subset dynamic programming for n=1,...,8. Normalizing by the constant term and centering gives all 24 displayed second, fourth and sixth moments exactly, without using Schur polynomials or hook lengths. This additional finite check supports the implementation; the universal conclusion uses the polynomial identities and derivation above.

## 3. Uniform positive lower bound

The sixth-order cosine polynomial is a global lower bound. The sixth derivative of the difference is 1-cos(x)>=0; its derivatives through order five vanish at zero. Six integrations and evenness prove the inequality for all real x, without a restriction on the size of X.

The resulting averaged polynomial L_n(t) is decreasing for t>0. The stated Cauchy--Schwarz inequality m4^2<=m2*m6 makes the derivative's quadratic factor positive. Equivalently, before averaging, the factor 1-u/6+u^2/120 is positive for every real u; hence substituting the upper bound 44/7 for 2*pi is legitimate even though the fourth-order term has positive sign.

I independently obtained

\[
L_n(44/7)-\frac16=
\frac{2841056n^8-35378248n^6+131969540n^4-186040404n^2+44118375}
{235298(4n^2-25)(4n^2-9)(4n^2-1)^2}.
\]

For n=2 this is 41809379/1111783050>0. For n>=3 the denominator is positive, and substitution n=x+3 gives exactly the nine positive coefficients displayed in the note. Thus the lower bound is proved for all n>=2, not extrapolated from the finite certificates. Together with the positive-density upper bound this proves U_n/6<|P_n(omega)|<=U_n and the stated fixed complex phase.

## 4. Factorial cost and asymptotic scope

The elementary asymptotic log A_n=(n^2/2)ln n-(3n^2/4)+O(n ln n), used also at 2n, gives

\[
\ln U_n=-n^2\ln n+
\left(\frac32+\ln(\pi/2)\right)n^2+O(n\ln n).
\]

The ratio |P_n(omega)|/U_n stays between 1/6 and 1, so taking its logarithm changes this expression by only O(1). Multiplication by the exact row-factorial factor A_n^2 then leaves ln(pi/2)n^2+O(n ln n). The sign of this leading coefficient is positive.

The note correctly distinguishes this sufficient row clearing from a least denominator after coefficient gcd cancellation. At the actual transcendental period the cleared value is not an integer. Additional denominator powers are required at a rational formal center. The result consequently supplies no contradiction with integrality and no rational approximation theorem.

The complete source has 2n rows while the ambient raw-column count is (n+1)^2. Its density tends to zero and the transverse dimension remains one. Nothing here establishes a corresponding estimate for the original high-dimensional half-density family.

## 5. Verification run

The standard-library verifier exited successfully with:

- 24 exact endpoint/Hankel identities and 8 Hilbert product identities;
- 3 universal polynomial moment identities;
- 15 centered-moment checks for the small-n cases;
- 35 independent Hilbert-series raw-moment checks;
- the positive-coefficient certificate for every n>=3, the separate n=2 case, and the Machin enclosure used for pi<22/7.

No correction to the final theorem or proof was required. This is an independent agent review, not external peer review.

## 6. Review of the finite least-denominator supplement

The subsequently added section 6 and its `least_denominator_diagnostic()` routine were reviewed separately. The routine reads the companion's certified rational polynomial coefficients for n=4,...,10 and records the fixture SHA-256. I ran that routine and matched every displayed denominator, row-multiplier quotient, and outward integer interval in the new table.

For a rational polynomial, the least common multiple of its reduced coefficient denominators is exactly the least positive integer multiplier making all coefficients integral. The additional gcd-one check confirms that each resulting integer polynomial is primitive in these seven cases. These are coefficient statements, not assertions about the least denominator of its value at an arbitrary rational argument.

The new actual-period evaluation independently encloses 2*pi using Machin intervals with 80 and 30 arctangent terms. The parity of the certified coefficients makes i^(-n) x^(-n) P_n(ix) a real polynomial. Exact rational interval Horner evaluation, followed by multiplication by the positive x^n interval, encloses its absolute value without a floating-point cancellation calculation. The interval excludes zero in each case. Scaling by the least coefficient denominator yields a lower endpoint greater than one for all seven cases.

For n=10 the least coefficient denominator is 489407583931698134384640000000000, and the certified cleared absolute value lies in [60061286897514303,60061286897514304]. The table therefore checks a genuinely smaller clearing than A_n^2 in this case. It does not establish a least-denominator asymptotic, monotonicity, or a Diophantine conclusion. The argument remains at the transcendental actual period; its cleared value is not asserted to be an integer.

The reviewed coefficient fixture had SHA-256 `ff5c3d51a87a9ea7e5527f0b87552205c4732a2d4a63a1bc56543de3a7a927ec`. Its polynomial certification belongs to the companion verifier; the new diagnostic correctly documents this dependency rather than treating the file's status string as an independent mathematical proof.

The final moment-verifier addition also passed: recursive corner removal counts all 30 Young-diagram cell-addition paths for partitions of sizes zero through six and agrees with k!/product(hook length). This checks the finite hook coefficients used here without relying on an unverified general hook formula.

No correction to either addition was required.
