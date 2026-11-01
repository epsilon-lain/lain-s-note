# Exact exterior transport coefficients and rounded lattice bases

2026-11-01. This note groups the boundary transport expansion at the determinant level. It identifies cancellations between row assignments exactly and constructs nonzero integer transport minors. It does not prove nonvanishing of the corresponding period-function minor or an endpoint bound for pi.

## 1. Group by conserved total derivative order

Keep the complete row packet

\[
L=\{(j,s,\beta):0\le j<K,\quad v_0s+\sum_i v_i\beta_i<H\},
\qquad c_i=v_i-v_0\ge0.
\]

For each center j and integer d=s+|beta|, the source rows in that group are indexed by the finite lower set

\[
B_d=\{\beta\in\mathbb Z_{\ge0}^m:
|\beta|\le d,\quad v_0d+c\cdot\beta<H\}.
\]

The word lower means that beta in B_d and 0<=gamma<=beta imply gamma in B_d. Set r_d=|B_d|. At a fixed center and d, a target index a represents the row (j,d-|a|,a), so admissible targets satisfy |a|<=d.

Let T(-zeta) be the shear transport from the preceding boundary note. Its only nonzero entries are

\[
T_{(j,d-|\beta|,\beta),(j,d-|a|,a)}
=\binom a\beta(-\zeta)^{a-\beta},\qquad a\ge\beta.
\]

It preserves j and d separately. Choose any full-size set J of target rows. Then

\[
\boxed{\det T_{L,J}=0
\quad\text{unless J has exactly }r_d\text{ rows in every group }(j,d).}
\tag{1}
\]

This follows from the block decomposition: a nonzero determinant permutation must match within each group. Thus the targets cannot be freely redistributed among total derivative orders or centers.

Order both source and target groups in the same way. In an admissible group write A_{j,d} for the target transverse indices and define

\[
C(B_d,A_{j,d})=
\det\bigl[\binom a\beta\bigr]_{\beta\in B_d,\ a\in A_{j,d}},
\quad
\Gamma_{j,d}=\sum_{a\in A_{j,d}}a-\sum_{\beta\in B_d}\beta.
\]

All nonzero determinant permutations carry the same error monomial. Hence

\[
\boxed{\det T_{L,J}
=\prod_{j,d}C(B_d,A_{j,d})(-\zeta)^{\Gamma_{j,d}}.}
\tag{2}
\]

If C is nonzero, every coordinate of Gamma is nonnegative: one nonzero matching has a>=beta entrywise and its coordinate increments sum to Gamma. Thus (2) has no negative powers. Prove the factorization first when all zeta_i are nonzero, then use its polynomial form to include zero coordinates. Each C is an integer independent of zeta and the period.

For a matrix of period rows B^0, the exact Cauchy--Binet regrouping is consequently

\[
\boxed{
\det(T(-\zeta)_{L,\Gamma_N}B^0)
=\sum_{J}\left[\prod_{j,d}
C(B_d,A_{j,d})(-\zeta)^{\Gamma_{j,d}}\right]\det B^0_J.}
\tag{3}
\]

Only target sets satisfying (1) need be summed. Formula (3) is an exact regrouping of the already proved boundary expansion, not an additional approximation. The matrix [I,R] after eliminating internal transport has the same full-size minors as T, because its left transformation has determinant one. Thus C also groups all boundary assignment cancellations into one coefficient for a given target set.

## 2. One coordinate: every increasing target set is nonzero

For m=1, B_d={0,...,r_d-1}. Let A={a_0<...<a_{r-1}} be any r distinct targets in {0,...,d}. Since binom(x,k) has degree k and leading coefficient 1/k!, the polynomial Vandermonde identity gives

\[
\boxed{C(B_d,A)=
\frac{\prod_{0\le i<j<r}(a_j-a_i)}{\prod_{k=0}^{r-1}k!}>0.}
\tag{4}
\]

After fixing the target order there is no further cancellation inside this transport coefficient. Signs or cancellations between the period minors in (3) remain possible and are not controlled by (4).

## 3. A sharp quadratic error power for partial boundary replacement

Continue with one transverse coordinate and one homogeneous group. Write its source set as B={0,...,b-1}, and let a target set A of size b contain exactly k outside indices (indices at least b). Then

\[
\boxed{\Gamma=\sum_{a\in A}a-\frac{b(b-1)}2\ge k^2.} \tag{5}
\]

To minimize the target sum, retain the b-k smallest source indices and choose the k smallest outside indices. The unique minimizing target set is

\[
A_{\min}=\{0,\ldots,b-k-1\}\cup\{b,\ldots,b+k-1\}.
\]

Its sum difference is k^2. It is admissible whenever k<=b and b+k-1<=d. Formula (4) says that its transport coefficient is nonzero, so the exponent k^2 is sharp at the transport level. For k=b the target set is consecutive, the Vandermonde quotient equals one, and the full-boundary coefficient is exactly (-zeta)^{b^2}.

Equivalently, multiply all boundary columns of [I,R] in this group by a formal variable eta. The coefficient of eta^k in any determinant formed from these compatible period rows is divisible by zeta^{k^2}. This follows directly from Cauchy--Binet: choosing k boundary columns contributes eta^k, and (5) applies to its target set. Further cancellations from the chosen period columns can increase the order; the statement does not exclude them.

For prescribed replacement counts k_{j,d} over the independent groups,

\[
\boxed{\text{the corresponding coefficient is divisible by }
\zeta^{\sum_{j,d}k_{j,d}^2}.}                    \tag{6}
\]

This is stronger than charging each replaced row independently: the outside target derivative rows must be distinct. It makes the exclusion constraint in the exterior determinant explicit. The exponent applies to a prescribed count in each group, before restoring the diagonal raw-column prefactor. For only a total count k over several groups, one must minimize the sum of squares over its possible distributions. It does not yet bound the remaining coefficient or period mixed minor. In several transverse coordinates, distinct fibers can admit first-order boundary moves, so a square of the total replacement count must not be asserted without a separate argument.

## 4. Several coordinates: exact cancellation can occur

Take m=2, v_0=1, (v_1,v_2)=(3,4), H=8, d=4. The source group is

\[
B_d=\{(0,0),(1,0),(0,1)\}.
\]

Choose targets A={(0,0),(1,1),(2,2)}, all satisfying |a|<=4. The integer coefficient matrix, in that order, is

\[
\begin{pmatrix}1&1&1\\0&1&2\\0&1&2\end{pmatrix}.
\]

Its determinant is zero. There are nonzero row assignments, but they cancel exactly because the two coordinate functions agree on these target nodes. In boundary form, keeping target (0,0) unchanged and replacing the other two rows produces

\[
R_{\{(1,0),(0,1)\},\{(1,1),(2,2)\}}
=\begin{pmatrix}
-\zeta_2&-2\zeta_1\zeta_2^2\\
-\zeta_1&-2\zeta_1^2\zeta_2
\end{pmatrix},
\qquad\det R=0.
\tag{7}
\]

A termwise product of nonzero boundary coefficients cannot certify a surviving determinant contribution. The binomial evaluation rank is the exact extra condition.

## 5. Nonzero multidimensional coefficients from rounded grids

Let B be any finite lower set of nonnegative integer multi-indices. For each coordinate choose strictly increasing nonnegative integer nodes

\[
f_i(0)<f_i(1)<\cdots.
\]

Set f(alpha)=(f_1(alpha_1),...,f_m(alpha_m)) and A=f(B), and order rows beta and columns alpha by one common linear extension of the coordinatewise order on B. Then

\[
\boxed{
\det\bigl[\binom{f(\alpha)}\beta\bigr]_{\beta,\alpha\in B}
=\prod_{\beta\in B}\prod_{i=1}^m
\frac{\prod_{k=0}^{\beta_i-1}(f_i(\beta_i)-f_i(k))}{\beta_i!}>0.}
\tag{8}
\]

Empty products are one. Individual displayed factors can be rational; the full product is an integer because it equals an integer determinant.

**Proof.** Replace the polynomial binom(x,beta) by the multivariate Newton polynomial

\[
N_\beta(x)=\prod_i\frac1{\beta_i!}
\prod_{k=0}^{\beta_i-1}(x_i-f_i(k)).
\]

Both polynomials have leading monomial x^beta/beta!, and N_beta is a linear combination of binom(x,gamma) with gamma<=beta. Since B is lower, this change of polynomial basis is triangular with diagonal one. At node f(alpha), N_beta vanishes if some alpha_i<beta_i. The resulting evaluation matrix is triangular, and its diagonal is precisely the product in (8). All differences are positive. This proves the formula.

For nonnegative integer translations t_i and arbitrary real dilation factors tau_i>=1, the choice

\[
f_i(n)=\lfloor\tau_i n\rfloor+t_i
\]

is strictly increasing. It gives an exact integer-node version of a fractional coordinatewise dilation; its node coordinates differ from tau_i n+t_i by less than one. No floating-point rounding argument is needed. If tau_i is an integer d_i, (8) simplifies to

\[
C(B,f(B))=\prod_i d_i^{\sum_{\beta\in B}\beta_i}.
\tag{9}
\]

In either case the transport error exponent is exactly

\[
\Gamma_i=\sum_{\beta\in B}(f_i(\beta_i)-\beta_i).
\tag{10}
\]

To apply (8) within a homogeneous group, one must separately verify |f(beta)|<=d for every beta in B_d. To use the resulting period rows on the weighted column space, the column support and actual period-minor nonvanishing must also be checked. A nonzero transport coefficient alone gives neither conclusion.

## 6. A warning about shell evaluation and footprint bounds

For a lower integer grid A, polynomial leading terms can give support bounds through its standard monomials. That conclusion cannot be transferred unchanged to a shell E=A\\B.

For example, in one coordinate let B={0,1}, A={0,1,2,3}, E={2,3}, and P(x)=x-2. The leading exponent is 1. Both shell nodes lie above it, but P is nonzero at only one shell node. Thus the shell upward-shadow size 2 is not a lower bound for the support size 1.

This does not disprove a shell packing theorem. It excludes one particular shortcut in its proof. The rounded-grid bases in (8) remain valid wherever their target sets fit, including translations into a shell.

## 7. Verification

The standard-library verifier checks the grouped transport factorization, the exact two-dimensional zero coefficient, arbitrary one-dimensional Pascal minors, multivariate rounded-grid determinants, integer-affine specializations and the shell-support counterexample. These finite identities accompany the algebraic proofs. They do not certify a high-dimensional mixed-minor estimate or a Diophantine property of pi.
