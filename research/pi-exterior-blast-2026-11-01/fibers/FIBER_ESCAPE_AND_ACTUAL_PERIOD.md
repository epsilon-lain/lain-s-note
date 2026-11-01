# Fiber factorization and a surviving two-exit term at the actual period

2026-11-01. The identities and finite nonvanishing certificate below do not settle bad approximability of pi. They check precisely how far the new one-coordinate quadratic error-power law extends.

## 1. Exact factorization over fibers

Let B be a finite lower set in the nonnegative integer lattice. Write an index as (x,y), separating its first coordinate. The fiber over y is B_y={0,...,n_y-1}. Choose n_y distinct nonnegative integers J_y, and let A be the union of the sets J_y times {y}. Order the y values by a common linear extension of coordinatewise order, and order x increasingly within each fiber. Then

\[
\det\left[\binom a\beta\right]_{\beta\in B,a\in A}
=\prod_y\frac{\prod_{x<x',\ x,x'\in J_y}(x'-x)}
{\prod_{k=0}^{n_y-1}k!}>0.                         \tag{1}
\]

Proof: the evaluation matrix is block triangular in y, since binom(y,beta') vanishes unless beta'<=y. Each diagonal block is a one-coordinate Pascal evaluation matrix. Its determinant is the displayed Vandermonde quotient.

If J_y contains k_y indices outside B_y, its first-coordinate displacement is at least k_y^2. The bound is sharp: retain 0,...,n_y-k_y-1 and add n_y,...,n_y+k_y-1. Hence the transport power in this fiber-preserving class is at least

\[
\zeta_1^{\sum_y k_y^2}.
\]

With F active fibers and k=sum_y k_y, the exponent is at least k^2/F, not generally k^2. This is a result for this class of target sets, not for arbitrary multivariable transport.

In particular, move only the last point n_y-1 to n_y in each of k selected fibers. The full grouped transport minor is exactly

\[
\boxed{\det T_{B,A}=(-\zeta_1)^k\prod_{\text{selected }y}n_y.} \tag{2}
\]

The integer coefficient is positive, so these assignments do not cancel after grouping by the target row set.

## 2. Compatibility with separated weights

Take m=2, theta=1/2, v_0=1, integer weights w_1=A>=1 and w_2=W>=A. Put c_1=2A-1 and c_2=2W-1. Choose integer h>=A, N=h and H=2h. In the conserved bucket s+|beta|=N, source types satisfy c_1 beta_1+c_2 beta_2<h. Their total-order cap is automatic because c_i>=1.

For every integer y with c_2 y<h, the fiber length is

\[
n_y=\left\lceil\frac{h-c_2y}{c_1}\right\rceil.
\]

Its first outside point (n_y,y) satisfies
\(h\le c_1n_y+c_2y<h+c_1\), and |a|=n_y+y<=h=N. The original column-type budget also holds:

\[
wa=\tfrac12(ca+|a|)<h+c_1/2\le2h=H.
\]

Thus (2) occurs among allowed source and target rows, with A and W as large and as separated as desired. At fixed weights there are only O(H) fibers in a two-coordinate bucket, versus O(H^2) rows. In m coordinates there are O(H^{m-1}) fibers per bucket and O(H^m) across all buckets, versus a full row count of order H^{m+1}. These first-order exits therefore cover a vanishing row fraction; they do not disprove a useful MH-scale bound for positive-fraction replacement.

## 3. A complete finite packet and an actual-period certificate

Take

\[
m=2,\ K=2,\ \theta=1/2,\ w_0=v_0=1,
\quad(w_1,w_2)=(2,5),\quad H=20.
\]

There are 156 complete source rows (j,s,beta), j=0,1, with s+4 beta_1+10 beta_2<20, and 242 allowed raw columns Y^h X_1^{alpha_1}X_2^{alpha_2} with h+2 alpha_1+5 alpha_2<=20.

In the j=0,N=10 bucket, the source types are
\(B_{10}=\{(0,0),(1,0),(2,0),(3,0),(0,1)\}\).
Change only these two row labels:

\[
(0,7,(3,0))\mapsto(0,6,(4,0)),\qquad
(0,9,(0,1))\mapsto(0,8,(1,1)).
\]

Both target rows have cost 22, outside the cutoff 20, but obey the column-type budget. Keep every other row. Formula (2), and identity transport in the other buckets, give

\[
\boxed{\det T_{L,J}=4\zeta_1^2.}                \tag{3}
\]

There are two boundary exits and exactly two error powers, rather than four. This coefficient already combines every assignment to that target set.

The same example also certifies both of the following determinants as nonzero on one explicit set of 156 raw columns:

- The original source matrix at r_1=2i(22/7), r_2=2i(355/113).
- The target-row matrix at the genuine common period omega=2 pi i.

The JSON certificate lists the selected columns and two nonzero residues modulo 65537. They are 21680 for the factorial-scaled period target and 35592 for the factorial-scaled rational source; for normalized ordinary jets the residues are 57507 and 36181. Equation (3) uses normalized jets. If both row lists are scaled by their time factorials, its coefficient becomes 252 zeta_1^2 instead; nonvanishing is unchanged. To explain the certificate, multiply each ordinary normalized time row by s!. The entry for column (h,alpha) becomes

\[
\binom\alpha\beta
\sum_{k\le\alpha-\beta,\ |k|\le s}
\binom{\alpha-\beta}{k}(jr)^{\alpha-\beta-k}
\frac{s!}{(s-|k|)!}h^{s-|k|}.                    \tag{4}
\]

At formal common period Z, these are integer polynomials in Z. Nonzero determinant at Z=1 modulo 65537 proves that the polynomial determinant is not identically zero. Known transcendence of pi then proves its nonvanishing at Z=2 pi i. This uses known transcendence and does not assert exponential periodicity at Z=1: the formal polynomial is formed after reduction at the true period.

For the rational-center matrix, use the ring homomorphism from Z[i,1/7,1/113] to the field modulo 65537 sending i to 256, since 256^2=-1 there. No homomorphism from the entire field Q(i) is asserted. The denominators 7 and 113 are invertible modulo this prime. The second nonzero residue proves nonvanishing of the Gaussian-rational determinant on the same columns. Row factorials do not affect nonvanishing. The source complete logarithmic time jets and ordinary time jets differ by a unitriangular transformation, so the original source format is also covered.

Finally zeta_1=(r_1-omega)/r_1 is nonzero. Thus (3) multiplies a genuinely nonzero period minor for actual rational approximations and a complete source packet. No unknown interpolation threshold is used in this finite certificate. The displayed weights also equal ceil(ln 7) and ceil(ln 113), though the proof only needs the specified weights.

## 4. Scope

The one-coordinate k^2 law is exact. Extending it unchanged to the total number of boundary exits in multiple coordinates would be false, even for a surviving term of a complete finite packet at the actual period.

A nonzero term can still be very small; different target sets can still cancel. No lower bound on the period minor, eventual high-dimensional low-cost family, or conclusion about pi's approximation constant follows. The remaining problem is the total contribution of compatible period minors, with the number of available fibers and their budgets accounted for.

## 5. The fully aggregated degree-two coefficient also survives

The previous nonzero target-set term alone would not rule out cancellation between different target sets. The following additional certificate addresses that distinction on the same 156 raw columns.

Use a formal marker eta multiplying every boundary column in [I,R], and treat zeta_1,zeta_2 as independent formal error variables. Define

\[
\Phi(Z,\eta,\zeta)=\det(B_L^0(Z)+\eta R(\zeta)B_O^0(Z)).
\]

Set zeta_2=0. The coefficient of eta^2 zeta_1^2 uses exactly two first-order boundary rows. With source rows multiplied by their time factorials, let B(Z) be the scaled source period matrix and E_1(Z) its first-order boundary matrix. For rho=(j,s,beta),

\[
(E_1)_\rho=
\begin{cases}
-s(\beta_1+1)\widehat B^0_{j,s-1,\beta+e_1},
 &s\ge1,\ 20-s-4\beta_1-10\beta_2\le3,\\
0,&\text{otherwise.}
\end{cases}
\]

The target row on the right is itself scaled by (s-1)!, explaining the additional factor s. There are 44 active rows. After multiplying by the common nonzero source row-factorial product, the desired coefficient is exactly [u^2]det(B(Z)+u E_1(Z)); higher error-degree terms in R cannot contribute because two boundary factors already supply degree at least two.

At Z=1 modulo 65537, the certificate gives

\[
\det B=21281,\quad
\operatorname{tr}(B^{-1}E_1)=-6,\quad
\operatorname{tr}((B^{-1}E_1)^2)=6.
\]

The determinant coefficient identity therefore yields

\[
[u^2]\det(B+uE_1)
=\tfrac12\det B\bigl((\operatorname{tr}B^{-1}E_1)^2
-\operatorname{tr}((B^{-1}E_1)^2)\bigr)
\equiv57067\pmod{65537}\ne0.                   \tag{5}
\]

An independent implementation forms entries by polynomial convolution, then reduces the determinant to 44 dimensions and uses exact interpolation at 45 field elements. It recovers the same coefficient, also checking two values against the full 156-dimensional determinant. The common trace values do not imply a rank-six projector: the first-order boundary matrix has rank 44 and its normalized square is not its negative.

The coefficient in (5) is an integer polynomial in Z. Thus it is nonzero at the genuine period 2 pi i by known transcendence. Consequently the full eta^2 coefficient, after combining all target sets, is not divisible by total error degree four: its zeta_1^2 term survives. This is stronger than (3).

This statement is about formal error-degree divisibility at the actual period. It does not lower-bound the coefficient's magnitude, rule out cancellation after substituting a particular nonzero zeta_2, or control the total determinant at eta=1. It is a counterexample to a universal multivariable square-power shortcut, not to bad approximability of pi.
