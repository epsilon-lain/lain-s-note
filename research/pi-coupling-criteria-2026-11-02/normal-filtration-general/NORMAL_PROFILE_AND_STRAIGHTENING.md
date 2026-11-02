# Normal profiles, exact block cuts, and when straightening removes all nonzero roots

2026-11-02. This note identifies precise conditions under which adding base derivatives or spatial directions does not create new coupling. It concerns the same formal-period raw matrices as the preceding work. It does not assert that all high-dimensional matrices factor, or that the conditions below are sufficient for a nonzero determinant.

## 1. Raw normal support persists when base derivatives are active

Use spatial variables X_1,...,X_d,T with T=X_(d+1). A row has label (j,s,nu,b), where j is a center index, s is a nonnegative time index, nu is a d-dimensional nonnegative normal derivative index, and b is the derivative index in T. A raw column is Y^h T^a X^kappa. The coefficient functional at a common formal period Z is

\[
[z^su^{\nu}u_T^b]\ e^{hz}(jZ+z+u_T)^a
\prod_{i=1}^{d}(jZ+z+u_i)^{\kappa_i}.
\]

Extracting the spatial coefficients first gives the exact entry

\[
\boxed{\binom{\kappa}{\nu}\binom ab
[z^s]e^{hz}(jZ+z)^{a-b+|\kappa-\nu|},}
\tag{1}
\]

with entry zero unless kappa>=nu coordinatewise and a>=b. Here binom(kappa,nu) is the product of the coordinate binomial coefficients. In particular, allowing b>0 does **not** by itself remove the normal upper-triangular support.

Consider any square M by M raw submatrix. Write R_nu for its number of rows with normal label nu, and C_kappa for its number of columns with normal label kappa. A zero count is allowed. If R_nu=C_nu for every normal label, order the labels by total degree, with identical orders within a degree. Equation (1) gives a block upper triangular matrix with square diagonal blocks

\[
\mathcal B_\nu(Z)_{(j,s,b),(h,a)}
=\binom ab[z^s]e^{hz}(jZ+z)^{a-b}.
\]

With these specified row and column orders,

\[
\boxed{\det\mathcal B(Z)=\prod_\nu\det\mathcal B_\nu(Z).}
\tag{2}
\]

The diagonal blocks are one-spatial-variable mixed time/base jet matrices. They need not be the old P_n blocks, whose base index b was always zero. No coordinate change or mixed-column replacement is required for (2).

## 2. Upper-set capacity, tight cuts, and a necessary surplus

Let Lambda be any finite coordinatewise lower set containing all row and column normal labels. An upper set U in Lambda means that nu in U and kappa>=nu in Lambda imply kappa in U. Put R(U)=sum_(nu in U) R_nu and C(U)=sum_(kappa in U) C_kappa. The rows in U vanish on every column outside U. Thus

\[
\boxed{\operatorname{rank}\mathcal B(Z)
\le M-R(U)+C(U).}
\tag{3}
\]

In particular nonzero determinants require R(U)<=C(U) for every upper set. Conversely these inequalities are sufficient for a perfect matching in the **normal-label support graph** kappa>=nu, with each row and column label repeated according to its count. To see this, the neighbor set of any collection of row labels is its upward closure; the corresponding upper-set inequality implies Hall's condition. This is only a support matching statement: the additional a>=b zeros, frequency restrictions, and algebraic cancellations can still make the actual matrix singular.

If R(U)=C(U), this is a tight cut. Ordering the U rows and columns first gives

\[
\mathcal B=\begin{pmatrix}\mathcal B_{U,U}&0\\ *&\mathcal B_{U^c,U^c}\end{pmatrix},
\qquad
\det\mathcal B=\det\mathcal B_{U,U}\det\mathcal B_{U^c,U^c}.
\tag{4}
\]

The sign is one for the displayed ordering; converting a previously fixed ordering to it may introduce only its explicit row/column permutation sign. Thus an apparent coupled matrix may still split at every saturated normal capacity cut. Formula (2) is the fully balanced instance.

There is also a quick scalar diagnostic. Define the total normal surplus

\[
E=\sum_{\rm columns}|\kappa|-\sum_{\rm rows}|\nu|.
\]

Every nonzero determinant permutation has kappa>=nu at each matched pair. Therefore E>=0 is necessary. If E=0, each matched pair must have kappa=nu, so either the profiles agree and (2) applies, or every determinant term is zero. In particular,

\[
\boxed{\text{A nonzero determinant with genuine transport between normal labels
must have }E>0.}
\tag{5}
\]

Positive surplus is necessary, not sufficient. In one normal dimension its exact decomposition is

\[
E=\sum_{r\ge1}\left(C(\{\kappa\ge r\})-R(\{\nu\ge r\})\right).
\]

It measures the total slack in the tail capacity inequalities. Tight tails still factor the matrix into separate bands.

## 3. The row and column coordinate changes have opposite weight conditions

The normal coordinates U_i=X_i-T straighten the diagonal center line. There are two different closure requirements, and they should not be conflated.

For columns, expansion of T^a U^kappa uses raw monomials

\[
T^{a+|\kappa-\lambda|}X^\lambda,\qquad0\le\lambda\le\kappa.
\]

A selected raw monomial space is invariant under this change and its inverse exactly when its index set is closed under these replacements, keeping h fixed. The expansion is integer triangular with diagonal one. For the full weighted column space w_0h+w_Ta+sum w_i kappa_i<=H, a sufficient condition at every H is

\[
w_T\le\min_i w_i.
\tag{6}
\]

For rows, let t=u_T and v_i=u_i-u_T. A coefficient row in the new variables expands as

\[
[t^b v^\nu]F(v+t,t)
=\sum_{|\eta|\le b}\binom{\nu+\eta}{\nu}
[u^{\nu+\eta}u_T^{b-|\eta|}]F(u,u_T).
\tag{7}
\]

Thus the row index set must be closed under moving base derivative units to normal directions:

\[
(s,\nu,b)\longmapsto(s,\nu+\eta,b-|\eta|).
\]

This is again an integer triangular change with diagonal one and an inverse with signed coefficients. The closure is exact as a statement about coefficient functionals on the full formal polynomial ring. For a complete strict weighted source v_0s+v_Tb+sum v_i nu_i<H, a sufficient condition at every H is

\[
v_T\ge\max_i v_i.
\tag{8}
\]

The two inequalities are also necessary if invariance is demanded for every positive height: a degree-one monomial or derivative at a height between the two unequal weights witnesses failure of the wrong ordering. A particular finite packet can satisfy closure accidentally, so the inequalities are not claimed necessary at a single H.

When v_i=w_i/theta in every spatial direction, the ambient column condition (6) and row condition (8) point in opposite directions. Both are automatic for every height only if all spatial weights equal. The previous sparse block construction instead used a specially closed selected column set with an expensive base variable. Consequently arbitrary weighted straightening cannot be treated as a free simultaneous row-and-column reduction.

If both finite closure conditions hold, applying both changes makes the matrix block diagonal by the adapted normal index. A nonzero square matrix must then have equal adapted row and column counts in every block, and its rank is the sum of the block ranks. This is a separate sufficient mechanism for block reduction; the direct raw support factorization (2) is available even without those closure hypotheses.

## 4. Time-to-base closure turns each balanced block into a monomial

Suppose a block in (2) has row labels (j,s,b) closed under

\[
(j,s,b)\longmapsto(j,s-k,b+k),\qquad0\le k\le s.
\tag{9}
\]

For a complete weighted source, condition (9) holds whenever v_T<=v_0. Define D=partial_z-partial_(u_T). On a block column e^(hz)(jZ+z+u_T)^a, D acts as multiplication by h. Replacing a normalized row partial_z^s partial_(u_T)^b/(s!b!) by D^s partial_(u_T)^b/(s!b!) is the integer unitriangular combination

\[
\sum_{k=0}^s(-1)^k\binom{b+k}{k}\,
\ell_{j,s-k,b+k}.
\tag{10}
\]

Entries with a<b remain zero. For a>=b, the transformed entry is exactly

\[
\frac{h^s}{s!}\binom ab(jZ)^{a-b}.
\]

Factoring Z^a from each column and Z^(-b) from each row, initially for Z!=0, gives

\[
\boxed{\det\mathcal B_\nu(Z)
=c_\nu Z^{\sum_{\rm block\ columns}a-\sum_{\rm block\ rows}b},
\qquad c_\nu\in\mathbb Q.}
\tag{11}
\]

If c_nu is nonzero, the exponent is nonnegative because the original determinant is polynomial. If the exponent would be negative, c_nu must be zero. Thus a balanced normal-profile determinant satisfying (9) is either zero or a single monomial in Z. It has **no nonzero root at all**. Activating base derivatives can therefore remove the previous near-period-root geometry rather than generate new coupling.

This does not apply to every source packet: the previous P_n family had no base derivatives, so moving a positive time derivative into the base direction was not allowed. It also does not exclude new roots when profiles are imbalanced or closure (9) fails.

## 5. A complete legal test family with active base derivatives

For an integer N>=2 use d=1, K=2, H=N, theta=1/3, v_0=v_T=1, v_X=2, w_0=1/4, w_T=1/3, w_X=2/3. All standard strict parameter inequalities hold, and the complete source is

\[
(j,s,\nu,b):\quad j=0,1,\quad s+b+2\nu<N.
\]

At each nu put B=N-2nu. Select raw columns X^nu Y^h T^a with 0<=h<B and 0<=a<2(B-h). There are B(B+1) rows and columns in that normal block. Every column satisfies its height budget, since

\[
h/4+a/3+2\nu/3
\le(2N-2\nu-1)/3<N.
\]

The normal profiles match and time-to-base closure holds. Thus its determinant is a monomial or zero by the general theorem. Direct exact checks for N=2,...,5 show nonzero constants in these instances, certify the predicted powers, and compare original raw matrices at additional rational Z. This finite nonvanishing statement is not extrapolated to all N. The general monomial-or-zero conclusion already holds for every N.

The same verifier checks exact row and column closure expansions, upper-set rank cuts, and the support-matching criterion on finite profiles. All checks passed: 1,080 row-coordinate and 1,080 column-coordinate identities; 989 normal support profiles, including 219 zero-surplus profiles; 112 complete weighted closure cases; four legal complete packets with 6, 14, 26 and 44 rows; 12 raw determinant comparisons; and 5,004 time-to-base row-change entries. The finite complete determinants are respectively Z^5, -Z^15, -Z^35 and Z^70 with the verifier's stated ordering. These checks supplement the proofs above. No result here establishes pi's bad approximability, a new rational approximation, or a critical half-density construction.


## Relation to prior checkpoints

The preceding [complete normal-block construction](../../pi-local-gram-2026-11-02/next-highdim-content-transfer/COMPLETE_NORMAL_BLOCK_PRODUCT.md) had no active base derivative. The present support theorem allows those derivatives and identifies when the factorization still holds. The [earlier mixed ideal cuts](../../pi-period-mechanism-2026-11-02/target-filtration-next/MIXED_NORMAL_IDEAL_CUTS.md) provide related necessary rank conditions. A finite imbalanced example and a monomial control are recorded in the [companion coupled certificate](../coupled-packet-search/COUPLED_RAW_CERTIFICATE.md).
