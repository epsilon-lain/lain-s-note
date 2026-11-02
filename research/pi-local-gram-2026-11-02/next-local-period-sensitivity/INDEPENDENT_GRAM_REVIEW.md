# Independent review of the normalized Gram disk argument

2026-11-02. Independently reviewed `LOCAL_PERIOD_DISC.md`, equations (1)--(10). The moment identity, complex operator-norm estimate, disk constants, exact rational-value exclusion, and resulting exp(-O(q log q)) scale check out. No mathematical blocker was found. This is an independent agent review, not external peer review or a claim about the optimal irrationality measure of pi.

## 1. The formal endpoint correction has the right sign

With the endpoint exponential replaced by an independent variable Y, integration by parts gives

\[
\mu_0(Z,Y)=(Y-1)/Z,\qquad
\mu_k(Z,Y)=Y/Z-k\mu_{k-1}(Z,Y)/Z.
\]

These moments are affine in Y. Their Y derivatives satisfy nu_0=1/Z and nu_k=1/Z-k*nu_(k-1)/Z, whose finite expansion is exactly the nu_k in the note. Thus

\[
\mu_k(Z,1)=\mu_k(Z,e^Z)+(1-e^Z)\nu_k(Z).
\]

The correction vanishes at the exponential period, but its derivative need not vanish. Keeping this correction is essential; replacing the derivative of the formal polynomial by the derivative of the actual endpoint determinant would be invalid. The proof retains the complete correction and controls its determinant, so it does not discard higher perturbation orders.

## 2. The Gram norm estimate works for complex vectors

Let v(u)=(1,u,...,u^(n-1))^T, let J be its positive real Gram matrix, and put psi(u)=J^(-1/2)v(u). Then psi has real entries and

\[
\int_0^1\psi(u)\psi(u)^T\,du=I.
\]

For complex unit vectors x,y, both x*psi and psi^T y have L^2 norm one. The second assertion follows from the reality of psi:

\[
\int_0^1|\psi(u)^T y|^2du=y^*Iy=1.
\]

Consequently Cauchy--Schwarz bounds

\[
\left|x^*\left(\int_0^1e^{Zu}\psi(u)\psi(u)^Tdu\right)y\right|
\le e^{\max(0,\operatorname{Re}Z)}.
\]

Taking the supremum over x,y proves the stated Euclidean operator-norm bound. No Hermitian symmetry of H_n(Z) is being assumed. The moment matrix may be singular at other Z; neither its inverse nor a lower bound for its determinant is needed.

The inverse-Hilbert diagonal identity is correct. The cofactor quotient from the Cauchy determinant gives

\[
(J^{-1})_{ii}=(2i+1)
\prod_{k\ne i}\left(\frac{i+k+1}{i-k}\right)^2,
\]

which equals the two squared binomial factors displayed in the note. Bounding the binomials by powers of two yields the stated generous B_n=2n*2^(8n). Positivity of J^(-1) justifies bounding its norm by its trace, at most n B_n.

For 5<=|Z|<=8, the finite reversed sum for nu_k is bounded by k!|Z|^(-k-1)e^|Z|. The inequalities e<3, |Z|>=5, and k<=2n-2 give the stated uniform entry bound 3^8(2n-2)!. Multiplying the norm bounds yields precisely

\[
L_n=2\cdot3^8n^3 2^{8n}(2n-2)!.
\]

## 3. The disk and constants are sufficient

For n>=2 and r_n=1/(8nL_n), the bound L_n>=n gives r_n<=1/32. This explicitly justifies e^(r_n)<2 in the estimate |1-e^Z|<=2r_n. The weaker statement r_n<=1 alone would not justify that estimate; the actual radius is well within the required range.

The disk remains in 5<|Z|<8 because 6<2pi<44/7. Factoring the positive Gram square root out of the determinant contributes exactly det J; no condition number or additional determinant factor is lost. The singular-value inequality |det M|<=||M||_2^n then gives the complete disk bound from the note.

Substituting the exact radius gives

\[
\frac{n^2r_n}{2\pi}<\frac1{48},\qquad
nr_n\le\frac18,\qquad
2nr_nL_n=\frac14.
\]

Thus 19/48 is a valid upper bound for the exponent. It is smaller than 1/2, and e^(1/2)<2. The bound |P_n(Z)|<2U_n therefore holds on the entire closed disk. Cauchy's estimate on radius r_n/2 disks inside it gives the optional derivative bound 4U_n/r_n, only on the stated smaller disk.

## 4. Rational exclusion and the asymptotic improvement

For even n>=4, the previously proved primitive polynomial R_n=d_nP_n is even and has unit leading coefficient. At ia/b with gcd(a,b)=1 and b>=2, multiplying by b^(n²/2) produces an integer coprime to b. Its value is therefore nonzero, with absolute value at least b^(-n²/2). This rational-value floor is independent of how small the polynomial is at the transcendental center.

Under 2d_nU_n<b^(-n²/2), the entire closed disk is excluded for such a rational center. The distance in the imaginary coordinate is twice the original approximation error, so

\[
\left|\pi-\frac pq\right|>r_n/2
=\frac1{32\cdot3^8n^4 2^{8n}(2n-2)!}.
\]

The factor 32 is correct. No prior proximity assumption on p/q is required. Cases b=1, equivalently q=1 or 2 for reduced p/q, are explicitly outside this denominator argument and can be treated separately.

The previous coefficient-denominator asymptotic gives

\[
\frac{\log(d_nU_n)}{n^2/2}
=-\log n+\log\kappa+o(1),
\qquad \kappa=e^{3/2}\pi^2/2.
\]

For fixed c>kappa and n the least even integer at least cb, the strict exclusion condition holds for all sufficiently large b. Stirling's estimate gives

\[
\log(2/r_n)=2n\log n+O(n),
\qquad \log(r_n/2)=-2cb\log b+O_c(b).
\]

Since q/2<=b<=q, this is an exp(-O(q log q)) lower-bound scale. It genuinely improves the preceding global-derivative comparison for this construction. It remains weaker than any fixed inverse-power lower bound, and supplies no bad-approximability conclusion. The note correctly separates the analytic threshold kappa from the slightly larger threshold obtained by using 44/7 in a fully rational sufficient certificate.

## Independent finite checks

Separate exact Fraction arithmetic checked 104 inverse-Hilbert diagonal identities for n=2,...,14, including their B_n bounds; 180 closed-form nu_k values against the endpoint-derivative recurrence for k=0,...,29 and Z in {-8,-6,-5,5,6,8}; and 13 radius/reciprocal/exponent substitutions. All passed. These checks supplement the all-size argument; they are not its basis. The main companion verifier supplies the reproducible formal-determinant and rational-height certificates.
