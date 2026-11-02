# Independent review of integer frequency dilation

2026-11-02. Reviewed `../frequency-dilation-next/FREQUENCY_DILATION_CONTENT.md`, including its completed central-limit section. No mathematical blocker found in the scaling, primitive arithmetic, or fixed-frequency comparison. This review does not address the unproved optimization with jointly varying frequency and dimension.

## The cited central limit theorem applies

Checked the primary source: Jonathan Breuer and Maurice Duits, *Central Limit Theorems for Biorthogonal Ensembles and Asymptotics of Recurrence Coefficients*, [arXiv:1309.6224v2](https://arxiv.org/pdf/1309.6224v2), Theorem 2.6. Its hypotheses hold for the probability measure du on [0,1]: compact support, essential support equal to that interval, and positive density almost everywhere. The test function f(u)=u is real and continuously differentiable on the real line. The theorem centers by the exact expectation and supplies a Gaussian variance from the positive-index Fourier coefficients.

Here the transformed test function is

\[
f\bigl((1+\cos\theta)/2\bigr)
=\frac12+\frac14e^{i\theta}+\frac14e^{-i\theta}.
\]

Therefore the only nonzero positive-index coefficient is 1/4, and the variance is 1/16. Reflection of every coordinate about 1/2 gives exact mean n/2 and a symmetric centered law. Distributional convergence then gives the characteristic-function limit

\[
\chi_n(t)\longrightarrow\exp(-t^2/32)
\]

at each fixed real t. No additional scaling by a power of n is appropriate. This is an application of the cited theorem, not a newly proved central limit theorem.

The subsequently added direct variance check also has the correct normalization. The shifted orthonormal Legendre recurrence has off-diagonal coefficient \(a_j=j/[2\sqrt{4j^2-1}]\). For the projection onto indices 0,...,n-1, the trace difference \(\operatorname{Tr}(P_nJ^2P_n-(P_nJP_n)^2)\) contains exactly the one crossing edge, with value \(a_n^2=n^2/[4(4n^2-1)]\). Its limit is 1/16, agreeing with the Fourier calculation rather than twice that value.

## Arithmetic and the primitive ratio

Under x=bz, every raw entry acquires b^(s-a). The total row and column exponents agree, giving exactly P_(n,b)(Z)=P_n(bZ). There is no leftover determinant factor.

For a primitive integer polynomial R_n(Z)=sum a_k Z^k of degree delta_n and lowest nonzero degree n, the dilated coefficient gcd is

\[
g_n(b)=\gcd_k(a_kb^k),\qquad b^n\mid g_n(b)\mid b^{\delta_n}.
\]

The second divisibility follows by multiplying a Bezout combination of the a_k by b^delta_n, with each term expressed as an integer multiple of a_k b^k. It is stronger than just a numerical size bound and justifies the primitive polynomial R_n(bZ)/g_n(b). The formula for the least integer denominator when h_n=1 is also correct; primitive normalization is allowed to use a rational multiplier and is a different operation.

For each fixed integer b>=2, both characteristic-function limits in the exact ratio are positive, so

\[
\frac{|\chi_n(2\pi b)|}{\chi_n(2\pi)}
\longrightarrow L_b:=e^{-\pi^2(b^2-1)/8}>0.
\]

For all sufficiently large n, this ratio is at least L_b/2. Consequently

\[
\frac{|R_{n,b}(2\pi i)|}{|R_n(2\pi i)|}
\ge \frac{L_b}{2}b^{n^2-\delta_n}
\ge \frac{L_b}{2}b^{n^2/2}\longrightarrow\infty.
\]

Both inequalities have the correct direction: dividing by g_n(b)<=b^delta_n gives the first, and the independently proved degree bound delta_n<=n^2/2 gives the second. Even the elementary degree bound delta_n<=n(n-1) would prove a weaker divergence proportional to b^n. The ratio diverges while both its numerator and denominator tend to zero, as follows from the all-dimension upper estimate below. The threshold may depend on b; the theorem provides no uniform comparison for a sequence b_n tending to infinity.

## Other scope checks

The unconditional upper estimate uses g_n(b)>=b^n and |chi_n|<=1, so its exponent n^2-n is correct. At fixed n, the leading term of R_n(b*2*pi*i), combined with g_n(b)<=b^delta_n, gives the stated liminf lower bound. The explicit sufficient fixed-n threshold follows from the usual leading-term triangle inequality with |2*pi*b|>6b. These checks do not promote coefficient-normalized smallness at a transcendental point to an integer contradiction or to a rational approximation theorem.

Combining the upper estimate with the established denominator theorem gives an unconditional dimension statement in both parities. Indeed, \(c_n=d_n/h_n\le d_n\le D'_n\), and the earlier factorial estimate gives

\[
\log(D'_nU_n)=-\tfrac12n^2\log n+O(n^2).
\]

Thus, for every fixed epsilon in (0,1/2), integer frequencies \(b_n\le n^{1/2-\varepsilon}\) satisfy

\[
\log|R_{n,b_n}(2\pi i)|
\le-\varepsilon n^2\log n+O(n^2)\longrightarrow-\infty.
\]

This needs no formula for the odd primitive content and no uniform characteristic-function lower bound for varying b_n. It is only an absolute upper estimate, so it does not conflict with the fixed-b comparison above.
