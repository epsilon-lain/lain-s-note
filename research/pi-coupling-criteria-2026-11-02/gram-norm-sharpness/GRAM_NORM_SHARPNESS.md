# The factorial scale of the normalized endpoint perturbation is intrinsic

2026-11-02. This note concerns the endpoint perturbation in the preceding local Gram disk proof. It proves matching leading logarithmic upper and lower bounds for its operator norm. Thus the factorial size in that norm estimate cannot be removed merely by changing the orthonormal polynomial basis. This is not a lower bound for the determinant, a sharp derivative estimate, or an upper bound for the largest true small-value disk.

## 1. The operator and the result

Let n>=2, omega=2*pi*i, and t=2*pi. Define the real positive Hilbert Gram matrix and the complex symmetric endpoint matrix by

\[
J_n=\left[\frac1{i+j+1}\right]_{i,j=0}^{n-1},\qquad
V_n(Z)=[\nu_{i+j}(Z)]_{i,j=0}^{n-1},
\]
\[
\nu_k(Z)=\sum_{\ell=0}^k
\frac{(-1)^\ell k!}{(k-\ell)!Z^{\ell+1}},\qquad Z\ne0.
\]

Write S_n=J_n^(-1/2), the real positive square-root inverse, and

\[
\mathcal B_n=S_nV_n(\omega)S_n.
\]

All norms below are operator norms induced by the complex Euclidean norm. The main conclusion is

\[
\boxed{\frac{(2n-1)!}{2(2\pi)^{2n-1}}
<\|\mathcal B_n\|_2
\le2\cdot3^8n^3 2^{8n}(2n-2)!
\qquad(n\ge9).}
\tag{1}
\]

In particular,

\[
\boxed{\ln\|\mathcal B_n\|_2=2n\ln n+O(n).}
\tag{2}
\]

The lower bound also has a fully rational weaker version obtained by replacing 2*pi in its denominator by T=44/7. The constants are intentionally nonoptimal.

## 2. One normalized monomial detects the large norm

Let e be the last coordinate vector. Put

\[
x=\frac{J_n^{1/2}e}{\sqrt{e^TJ_ne}}.
\]

This is a real unit vector because J_n is real positive definite. Since e^TJ_ne=1/(2n-1),

\[
x^*\mathcal B_n x
=\frac{e^TV_n(\omega)e}{e^TJ_ne}
=(2n-1)\nu_{2n-2}(\omega).
\]

For any complex matrix B and unit vector x, |x*Bx|<=||B||_2. Thus

\[
\boxed{\|\mathcal B_n\|_2
\ge(2n-1)|\nu_{2n-2}(\omega)|.}
\tag{3}
\]

Equivalently, the test polynomial is the L^2[0,1]-normalized monomial sqrt(2n-1)*u^(n-1). No numerical matrix square root, inverse conditioning estimate, or Hermitian symmetry of V_n is required.

## 3. The endpoint moment is a factorial times a truncated exponential

Reverse the finite sum defining nu_k:

\[
\nu_k(Z)=(-1)^k k!Z^{-k-1}
E_k(-Z),\qquad E_k(w)=\sum_{j=0}^k\frac{w^j}{j!}.
\tag{4}
\]

At omega, exp(-omega)=1. If k+2>T=44/7>t, the exponential series tail obeys

\[
|E_k(-\omega)-1|
\le\sum_{j=k+1}^{\infty}\frac{t^j}{j!}
<\mathcal E_k(T):=
\frac{T^{k+1}}{(k+1)!}\frac1{1-T/(k+2)}.
\tag{5}
\]

The ratio of consecutive tail terms is at most T/(k+2), which gives the geometric majorant. The finite rational base check is

\[
\mathcal E_{16}(44/7)<\frac16<\frac12.
\tag{6}
\]

For k>=16 the first factor T^(k+1)/(k+1)! decreases when k increases, because T/(k+2)<1. The positive reciprocal factor 1/(1-T/(k+2)) decreases as well. Therefore the bound in (6) holds for every k>=16. In particular |E_k(-omega)|>1/2. Taking k=2n-2 in (3)--(4) proves the lower bound in (1) for all n>=9.

More precisely, the tail in (5) tends to zero, so

\[
|\nu_{2n-2}(\omega)|
=\frac{(2n-2)!}{t^{2n-1}}(1+o(1)).
\tag{7}
\]

Equation (7) is an asymptotic for one matrix entry, not an asymptotic equality for the full operator norm.

## 4. The matching upper logarithmic scale

The earlier Gram disk proof gives the upper bound in (1). Its elementary estimates can be recalled here. The inverse-Hilbert diagonal formula gives

\[
(J_n^{-1})_{ii}
=(2i+1)\binom{n+i}{n-i-1}^{\!2}\binom{2i}{i}^{\!2}
\le2n\,2^{8n}.
\]

Since J_n^(-1) is positive, ||J_n^(-1)||_2 is at most its trace, hence at most 2n^2*2^(8n). At 5<=|Z|<=8, (4) and the absolute exponential series give

\[
|\nu_k(Z)|\le k!|Z|^{-k-1}e^{|Z|}
\le3^8(2n-2)!\quad(0\le k\le2n-2).
\]

The norm of an n-by-n matrix whose entries have magnitude at most M is at most nM. Submultiplicativity therefore gives exactly the upper bound in (1). Stirling's estimate applied to the two sides of (1) now proves (2): each logarithm is 2n*ln n+O(n).

## 5. What basis changes can and cannot improve

The functions in the vector S_n(1,u,...,u^(n-1))^T are an orthonormal real polynomial basis for L^2[0,1] restricted to degree below n. Another real orthonormal basis is obtained by a real orthogonal matrix O, and the endpoint matrix becomes

\[
O\mathcal B_nO^T.
\]

Its Euclidean operator norm is unchanged. The same conclusion holds for a complex orthonormal basis, with a unitary congruence \(U\mathcal B_nU^T\): both U and U^T are unitary. Thus the factorial scale is a property of this bilinear endpoint functional in the fixed L^2 geometry, not merely of the original monomial coordinates. A nonorthogonal change must also transform the Gram metric; measuring it in an unrelated Euclidean norm would change the comparison being made.

Consequently, any valid scalar upper bound Lambda_n for this normalized operator, even just at omega, must satisfy

\[
\ln\Lambda_n\ge2n\ln n-O(n).
\]

If one keeps the preceding proof's sufficient-radius prescription r_n=1/(8n*Lambda_n), sharpening the norm bound alone cannot remove the leading 2n*ln n in ln(1/r_n). The existing upper bound already matches that leading scale.

This does **not** prove that the largest genuine disk with a small determinant has that radius. A matrix norm can overestimate determinant variation, and directional or exterior-power cancellations could allow a larger disk. It likewise does not prove a lower bound for |P_n'(omega)|, optimality of the preceding rational-approximation bound, or an obstruction to a genuinely different high-dimensional construction. Those questions require information beyond this norm.

## Verification

The standard-library companion checks the exact finite tail premise, its monotonicity on a sample range, the reversed exponential identity in Gaussian rational arithmetic, the inverse-Hilbert formulas, and orthogonality of the shifted Legendre polynomial basis in the same Gram metric. It does not approximate pi or a matrix square root. The all-size tail monotonicity, norm inequality, and asymptotic statement follow from the proof above.


## Preceding comparison

The upper norm and sufficient-radius prescription are in the [normalized Gram disk proof](../../pi-local-gram-2026-11-02/next-local-period-sensitivity/LOCAL_PERIOD_DISC.md). The lower bound in this note concerns that same operator and metric; it does not turn the sufficient disk into an optimal one.
