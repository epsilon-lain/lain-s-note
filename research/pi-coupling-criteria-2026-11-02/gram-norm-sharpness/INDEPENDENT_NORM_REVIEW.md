# Independent review of the normalized endpoint-norm lower bound

2026-11-02. Reviewed the proposed lower bound for the endpoint perturbation matrix after normalization by the positive Hilbert Gram matrix. The normalized monomial test, exponential-tail bound, all-size asymptotic conclusion, and orthonormal-basis invariance are valid. No mathematical blocker was found. This is an independent agent review, not external peer review.

Let omega=2*pi*i, t=2*pi, J=[1/(i+j+1)] for 0<=i,j<n, and S=J^(-1/2). Let V have entries nu_(i+j)(omega), and write B=S V S. The moments are

\[
\nu_k(\omega)=(-1)^k k!\,\omega^{-k-1}
\sum_{j=0}^k\frac{(-\omega)^j}{j!}.
\]

## The monomial test gives the stated norm lower bound

Let e be the last standard coordinate vector and put

\[
c=\sqrt{2n-1}\,J^{1/2}e.
\]

This vector is real, and its Euclidean norm squared is

\[
c^*c=(2n-1)e^TJe=1.
\]

Consequently

\[
c^*Bc=c^TBc=(2n-1)e^TVe
=(2n-1)\nu_{2n-2}(\omega).
\]

The general inequality |c* B c|<=||B||_2 therefore proves

\[
\|B\|_2\ge(2n-1)|\nu_{2n-2}(\omega)|.
\]

No Hermitian symmetry or positivity of B is assumed. The reality of c is what identifies the adjoint quadratic expression with the bilinear moment expression. Equivalently, this tests the endpoint bilinear form on the L^2[0,1]-normalized monomial sqrt(2n-1)*u^(n-1).

## The explicit tail estimate is sufficient from n=9 onward

Because exp(-omega)=1, the difference between the displayed truncated sum and 1 is bounded by the absolute exponential tail. With T=44/7>t, define

\[
E_k(T)=\frac{T^{k+1}}{(k+1)!}
\frac1{1-T/(k+2)},\qquad k+2>T.
\]

Successive terms of the positive tail are bounded by the geometric ratio T/(k+2), so the true tail is strictly smaller than E_k(T). The exact base value is

\[
E_{16}(44/7)
=\frac{24090891970712521146368}
{149394938057999106188625}<\frac16<\frac12.
\]

Moreover,

\[
\frac{E_{k+1}(T)}{E_k(T)}
=\frac{T}{k+2}
\frac{1-T/(k+2)}{1-T/(k+3)}<1
\qquad(k\ge16).
\]

Both factors are positive and smaller than one. Thus for k=2n-2 and every n>=9, the truncated sum has magnitude greater than 1/2. Substitution in the monomial test yields

\[
\|B\|_2>
\frac{(2n-1)!}{2(2\pi)^{2n-1}}
>\frac{(2n-1)!}{2(44/7)^{2n-1}}.
\]

Strictness is retained at each stated step. The second inequality gives a fully rational lower certificate.

## Asymptotic interpretation and its precise limit

The preceding local-disk argument supplies the explicit upper bound

\[
\|B\|_2\le
2\cdot3^8n^3 2^{8n}(2n-2)!.
\]

Stirling's estimate applied to the lower and upper bounds proves

\[
\log\|B\|_2=2n\log n+O(n).
\]

In particular the leading n*log(n) scale in this operator norm is not merely an artifact of the earlier loose inverse-Hilbert estimate. The result does not identify the exact linear-in-n correction.

Changing a real orthonormal polynomial basis replaces B by O B O^T with O orthogonal, preserving its Euclidean operator norm. A complex orthonormal change produces U B U^T with U unitary; U^T is also unitary, so the same invariance holds. This conclusion concerns the fixed L^2[0,1] norm. A nonunitary change must account for the changed domain and range norms and their condition numbers.

Finally, a large perturbation norm does not force a proportionally large determinant change. Cancellation, rank, or special spectral structure can make a determinant far less sensitive than a norm bound suggests. The result therefore does **not** prove that the previous certified disk radius is optimal, that the actual derivative has this scale, or that another determinant argument cannot yield a better rational-approximation estimate. It supplies no bad-approximability conclusion for pi.

## Independent finite checks

Exact rational arithmetic separately checked the displayed E_16 value and its comparison with 1/6; 84 consecutive tail-ratio identities and monotonicity checks for k=16,...,99; and 98 normalized-monomial scalar identities for n=2,...,99. These finite checks supplement the general algebraic arguments above rather than establishing them by extrapolation.
