# Independent review of rational-specialization cost

2026-11-02. Independently audited `RATIONAL_SPECIALIZATION_COST.md`, including the exact denominator, disk derivative estimate, comparison constants, height threshold, and resulting asymptotic lower bound. No substantive error was found. The only suggested addition is an explicit justification of the parity statement; that justification is given below.

## 1. Parity and the exact value denominator

The commutator representation is P_n(Z)=det(E T(Z)-T(Z)E), where T(-Z)=T(Z)^(-1) and det T(Z)=1. Therefore

\[
E T^{-1}-T^{-1}E
=-T^{-1}(ET-TE)T^{-1},
\qquad P_n(-Z)=(-1)^nP_n(Z).
\]

For even n=2r, R_n=d_nP_n is consequently an even integer polynomial. Its degree N_n=2r² is even, and its leading coefficient is (-1)^r. These facts justify treating R_n(ia/b) as a real rational number.

After multiplication by b^N_n, every term except the leading one is divisible by every prime dividing b. The leading term is

\[
(-1)^r i^{N_n}a^{N_n}=a^{N_n},
\]

since i^(2r²)=(-1)^r. Because gcd(a,b)=1, the numerator is coprime to b, hence to b^N_n. Its reduced denominator is exactly b^N_n. For b>=2 the numerator cannot be zero. The resulting lower bound |R_n(ia/b)|>=b^(-N_n) is valid without a nonvanishing hypothesis. Evenness does not halve this denominator: rewriting in Z² replaces the rational denominator by b² at the same time as it halves the degree.

The reduction a=2p/g, b=q/g, g=gcd(2,q), is correct for a reduced p/q, including negative p. The b=1 exception is properly separated.

## 2. Derivative estimate

For |Z|<=8 the entrywise bounds 9^a and a*9^(a-1) follow directly from the finite coefficient expansion. Rows at j=0 have zero derivative and cause no additional term. In a determinant permutation, the sum of the spatial exponents of all 2n columns is

\[
2\sum_{a=0}^{n-1}a=n(n-1)=D_n.
\]

Thus the product rule bounds each differentiated permutation product by D_n 9^(D_n-1). Summing the (2n)! permutations and using the nonnegative binomial term D_n9^(D_n-1)<=10^D_n yields exactly the stated M_n=d_n(2n)!10^D_n. No extra factor 2n is missing; it is already included through the sum of the column degrees.

The line segment lies inside the disk when its length is at most one, because 2pi+1<8. The use of a complex line integral introduces precisely the segment length, with no additional constant.

## 3. Comparison inequality and its range

Under d_nU_n<=1/(2b^N_n), the value at the actual period is at most one half of the rational nonzero-value lower bound. Hence

\[
M_n|2ip/q-2\pi i|\ge\frac1{2b^{N_n}}.
\]

Dividing by 2M_n proves the claimed

\[
\left|\pi-\frac pq\right|
\ge\frac1{4b^{N_n}d_n(2n)!10^{D_n}}.
\]

The factor 4 is correct. If the segment length exceeds one, the error exceeds 1/2, while the displayed right side is smaller than 1/2, so the same conclusion holds. Consequently the note correctly removes the proximity assumption after deriving the local inequality.

For each fixed b, the negative term -(n²/2)log n dominates +(n²/2)log b in log(b^N_n d_nU_n), proving eventual satisfaction of the condition. Replacing 2pi by 44/7 gives a rational sufficient certificate and also terminates for each fixed b.

## 4. Height threshold and the resulting scale

Writing C=3/4+log(pi/sqrt(2)), the established estimates give

\[
\frac{\log(d_nU_n)}{N_n}
=-\log n+2C+o(1),\qquad
e^{2C}=\frac{e^{3/2}\pi^2}{2}=\kappa.
\]

The same limit holds for the actual value, because it differs from d_nU_n by a factor between 1/6 and 1. For n the least even integer at least cb, one has n/b -> c, and therefore

\[
\frac{\log(b^{N_n}d_nU_n)}{N_n}\longrightarrow\log(\kappa/c).
\]

This proves eventual success for c>kappa, and proves eventual failure of the analogous actual-value baseline condition for c<kappa. No conclusion at c=kappa is asserted. The result concerns only this baseline comparison, as the note correctly emphasizes.

The rational upper bound 44/7 has its own slightly larger threshold e^(3/2)(22/7)²/2. This does not invalidate the note: the analytic statement with c>kappa concerns the true condition, whereas the exact rational search is separately asserted to terminate for each fixed b. One should not identify the two thresholds when reporting numerical certificates.

Finally, n=cb+O(1) and

\[
\log d_n=N_n\log n+O(n^2),\qquad
\log((2n)!)=O(n\log n),\qquad D_n\log10=O(n^2).
\]

The logarithm of the explicit right side is consequently

\[
-N_n\log b-N_n\log n+O_c(b^2)
=-c^2b^2\log b+O_c(b^2).
\]

The coefficient c², rather than c²/2, is correct: the rational evaluation and the coefficient clearing each contribute a leading logarithm. Since q/2<=b<=q, the exp(-O(q² log q)) description is valid as an asymptotic lower-bound scale. It neither estimates the actual error from above nor supplies the fixed positive multiple of q^(-2) required for bad approximability.

## 5. Independent exact checks

Using the existing certified coefficients for n=4,6,8,10, exact rational arithmetic independently confirmed the degree, parity, unit leading coefficient, and integer coefficients. It then checked the reduced denominator identity and nonzero-value lower bound at 4,420 rational evaluations: for each n, every reduced p/q with 3<=q<=30 and -2q<=p<=2q, plus 355/113. In particular the n=10 example has denominator exactly 113^50.

For each of those four polynomials, the stronger directly computable coefficient majorant sum_k k|[Z^k]R_n|8^(k-1) was also verified to be no greater than the proposed M_n. These finite checks supplement the all-dimension algebraic arguments above; the threshold and asymptotic conclusions were checked from those arguments, not fitted from numerical data.
