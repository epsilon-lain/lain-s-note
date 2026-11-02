# Independent review of normal profiles and straightening

2026-11-02. Independently reviewed `NORMAL_PROFILE_AND_STRAIGHTENING.md`, including the raw support formula, upper-set capacity criterion, zero-surplus factorization, finite coordinate-change closure, and time-to-base monomial theorem. No mathematical blocker was found. The conclusions are conditional structural statements for the specified matrices, not impossibility statements for all high-dimensional constructions.

## Raw support, capacity, and factorization

Extracting spatial coefficients first gives exactly

\[
\binom\kappa\nu\binom ab
[z^s]e^{hz}(jZ+z)^{a-b+|\kappa-\nu|},
\]

and gives zero when either coordinatewise kappa>=nu or a>=b fails. Consequently active base derivatives do not remove the normal upper-triangular support. With equal row and column counts at every normal label, ordering by total normal degree produces square diagonal blocks. Equal total degrees and coordinatewise inequality force equal labels, so the stated block product is exact with the displayed ordering.

For an upper set U, its rows are zero on columns outside U. Their rank contribution is at most C(U); all other rows contribute at most M-R(U). This proves the rank bound and the necessity of the upper-set inequalities. The support-only converse is also valid: for any collection of individual rows, take the upward closure U of its normal labels. Its neighboring columns are exactly those in U, while its number of rows is at most R(U)<=C(U). This is Hall's condition. Additional coefficient zeros and cancellations remain outside this support criterion, as the note emphasizes.

If an upper-set cut is tight, ordering those rows and columns first gives square diagonal blocks and an actual determinant factorization. The note records the possible permutation sign when changing an existing ordering. For zero total surplus, every nonnegative matched difference |kappa|-|nu| must vanish, forcing kappa=nu in every surviving permutation. Hence either all profiles agree or the determinant is zero. Positive surplus is necessary for transport between labels but not sufficient for nonvanishing or new roots.

## Opposite finite closure conditions

For columns, expanding T^a(X-T)^kappa generates precisely the replacements described in the note, with nonzero integer binomial coefficients and diagonal coefficient one. A raw monomial span is invariant exactly when it contains every resulting monomial. The inverse change uses the same index replacements, so the stated closure suffices in both directions.

For rows, expand a monomial u^(nu+eta) u_T^(b-|eta|) after u_i=v_i+t and u_T=t. Its coefficient of v^nu t^b is binom(nu+eta,nu). This proves formula (7), with no extra multinomial coefficient: the powers of the single variable t add after the independent binomial expansions. The row change is triangular by decreasing base index b, with diagonal coefficient one and an inverse obtained by signed substitution.

The weighted column replacement lowers weight exactly when w_T<=w_i, while the weighted row replacement lowers weight when v_T>=v_i. Degree-one monomials or coefficient functionals at a height strictly between unequal weights prove necessity if invariance is required at every positive height. A finite selected packet may still be closed for other reasons. With v_i=w_i/theta, simultaneous ambient invariance at all heights indeed forces equal spatial weights. This observation does not prevent the specially closed selected spaces constructed previously.

When both finite closures hold, adapted columns contain the factor U^kappa, and adapted rows extract the coefficient v^nu. On the center curve U=v, so the resulting matrix is block diagonal in the normal label. A nonzero square matrix must have equal row and column counts in each such block. This verifies the claimed reduction without treating an unavailable row change as free.

## Time-to-base closure and the monomial conclusion

The normalized differential expansion is

\[
\frac{(\partial_z-\partial_{u_T})^s\partial_{u_T}^b}{s!b!}
=\sum_{k=0}^s(-1)^k\binom{b+k}{k}
\frac{\partial_z^{s-k}\partial_{u_T}^{b+k}}{(s-k)!(b+k)!}.
\]

Thus the integer coefficient in (10), including its factorial normalization, is correct. Under the asserted closure every required row exists, and the transformation is unitriangular by time index, so its determinant is one. The differential difference acts as multiplication by h on the block column, leaving entry h^s/s! times binom(a,b)(jZ)^(a-b), with zero retained for a<b.

For Z!=0, factoring Z^a from each column and Z^(-b) from each row is legitimate. The remaining entries are constant in Z; zero entries stay zero. This proves a single common Z power for the determinant. Polynomiality handles the possible negative formal exponent exactly as stated: its constant coefficient must then vanish. A nonzero balanced-profile determinant with this closure has no nonzero root.

The test family's complete source count is B(B+1) in each normal block. Its selected columns have the same count, and its budget estimate follows from

\[
h/4+a/3+2\nu/3
\le (2N-2\nu-1)/3+(1/4-2/3)h<N.
\]

The all-N conclusion is monomial or zero. The separate finite nonvanishing checks are not promoted to an all-N nonvanishing theorem. This scope distinction, and the fact that the old P_n rows fail time-to-base closure, are both preserved correctly.
