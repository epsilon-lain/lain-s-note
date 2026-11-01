# Independent review of the small-escape Hardy bound

2026-11-01. Reviewed `../SMALL_ESCAPE_GRAM_BOUND.md` independently of its computational verifier. No mathematical error was found in its stated bounds. This is a review of an absolute upper bound for period minors, not a proof of a pi approximation claim.

## Row exchange and normalization

The retained row set I is a subset of the independent source rows, so its Gram matrix is invertible. Schur complementation bounds the extra determinant by the product of squared norms of the new rows. Jacobi's identity then gives det G_I=det G_L times the complementary principal minor of G_L^{-1}. Bounding that minor by lambda_min(G_L)^{-k} proves the finite row-exchange inequality. This does not require bucket conservation.

The rows throughout are ordinary Taylor coefficient functionals [t^s] f(j omega+t), not unnormalized derivatives or sparse logarithmic jets. This convention is essential for the uniform row-norm estimate and for the cardinal polynomials. The ordinary and logarithmic complete source jet blocks have the same determinant because their full triangular change has diagonal one; no such change is silently applied to a sparse target set.

## The row-norm constant

With R=14K, A=(K-1)|omega| and r=(R-A)/2, one has r>1. The upper bound |omega|<44/7 gives

\[
\frac{A+r}{R}=\frac{R+A}{2R}<\frac{71}{98}.
\]

The Hardy evaluation norm at z is (1-|z|^2/R^2)^(-1/2). Cauchy on the circle |t|=r adds r^(-s), no factorial. Squaring and using r>1 gives 1/(1-(71/98)^2)=9604/4563, as stated. This estimate applies to all nonnegative target time orders at the fixed centers.

## Hermite interpolation and the inverse Gram bound

For the cardinal polynomial p_{j,s}, the product Lambda_j has zeros of order n at all other centers. At its own center, multiplication of Lambda_j by the truncated reciprocal equals one modulo t^(n-s). The leading t^s therefore gives exactly t^s modulo t^n. This proves the required cardinal identities for all n consecutive ordinary normalized jets.

On |t|=|omega|/2, each reciprocal factor is bounded by (1-1/(2d))^(-n). The product is at most [(j+1)(K-j)]^n because 2d/(2d-1)<=(d+1)/d. The coefficient bound is therefore C_j^n(2/|omega|)^l.

On |z|=R, every |z-x_l| is less than L_K=(142K-44)/7. In each summand of p_{j,s}, s+l<=n-1; hence L_K^(s+l)(2/|omega|)^l<=L_K^n. Together with

\[
\prod_{l\ne j}|x_j-x_l|=|omega|^{K-1}j!(K-1-j)!
\]

this gives norm(p_{j,s})<=n Q_K^n. The factorial estimate follows from binom(K-1,j)<=2^(K-1), and |omega|>6 supplies the factor 3 in Q_K. The K=1 formula Q_1=14 also works directly. In all cases Q_K>1, so the bound is monotone in n and may be replaced by n_*.

If Z is the Hardy jet map and P the cardinal interpolation map, ZP=I. Decompose P=Z^*(ZZ^*)^(-1)+W with ZW=0. Orthogonality gives

\[
P^*P=(ZZ^*)^{-1}+W^*W.
\]

Thus the inverse Gram operator norm is no larger than trace(P^*P)<=K n^3 Q_K^(2n). The full source Gram is block diagonal in the transverse type. Taking its inverse operator norm uses the maximum of the block norms, not their sum, so no extra number-of-types factor is needed.

The claimed logarithmic simplification also checks. For K>=2, L_K<21K and C_j<=(K+1)^2/4<=K^2, while (K-1)!>=((K-1)/e)^(K-1). These imply

\[
\log Q_K<\log21+3\log K+(K-1)(\log7+1)
+(K-1)\log(K/(K-1))
<3K+3\log K+2.
\]

## Absolute period-minor conclusion

The coefficient functions f_{a,P}(z) give exactly the actual-period rows because exp(h j omega)=1 for integral h and j. For a raw monomial, its coefficient has size at most exp(hR)(2R)^|alpha| on the radius-R circle, giving the stated weighted-column cost g_R. The number of nonzero types contributes only log N_A/(2H).

Combining that factorization with the source Gram product and the row-exchange bound yields the uniform absolute estimate g_R-S_R+rho log(Q_K)/v_0+o(1). Therefore every k=o(M) selection retains the original Hardy saving to leading MH order at fixed parameters.

This result is compatible with the generic source-relative spectral spike in `ESCAPE_GRAM_ENVELOPE.md`: that example compares with a possibly much smaller chosen finite polynomial determinant B. The present estimate uses the full Hardy evaluation Gram as its baseline and does not divide by that selected B. Neither result can be substituted for the other without tracking this distinction.

The proof does not cover a simultaneous H-dependent growth of dimension or weights; the fixed-parameter order of limits is explicit. A positive escape fraction still pays the displayed linear loss, and all transport and arithmetic scaling costs remain separate.
