# Independent review of the nearest complex root scale

2026-11-02. This review checks `NEAREST_ROOT_SCALE.md` using the
separately audited logarithmic derivative theorem, the established period
lower bound, and the established uniform complex disk estimate.

## Verdict

Both nearest-root conclusions are valid:

\[
\log\delta_n^Z=-2n\log n+O(n),\qquad
\log\delta_n^Y=-2n\log n+O(n).
\]

The distances refer to two different complex polynomial zero sets. The first
varies `Z` with the endpoint fixed at `Y=1`; the second varies the independent
endpoint `Y` with `Z=2pi i` fixed. No identification between these sets is used.

## Z roots: the degree supplies the missing global information

An entry in column `(h,a)` has degree at most `a`. Every determinant term uses
each column exactly once, so `D_n=deg P_n <= n(n-1)` is a valid bound without
relying on a sharper degree formula. At zero the two row blocks coincide,
whereas `P_n(2pi i)` is nonzero by the period lower bound. Thus the polynomial
is nonconstant and has a finite nonempty multiset of complex roots.

Factoring this polynomial, with every zero repeated according to its
multiplicity, gives

\[
\frac{P_n'(\omega)}{P_n(\omega)}
=\sum_{j=1}^{D_n}\frac1{\omega-\zeta_j}.
\]

If `delta_n^Z` is the least root distance, the triangle inequality gives
`Q_n <= D_n/delta_n^Z`, where `Q_n=|P_n'/P_n|`. Consequently

\[
\delta_n^Z\le\frac{n(n-1)}{Q_n}
\]

whenever `Q_n>0`, which the derivative theorem guarantees eventually. No
absence of cancellation is assumed. This implication needs the finite-degree
bound; it would be invalid for arbitrary nonvanishing analytic functions.

For the lower bound, set `f(w)=P_n(omega+r_n w)/(2U_n)`. The disk theorem gives
`|f|<1` on the closed unit disk and the period theorem gives `|f(0)|>1/12`.
For a zero `a` inside the disk, division by its Blaschke factor gives an
analytic function bounded by one, hence `|f(0)|<=|a|`. This proves

\[
\delta_n^Z>r_n/12.
\]

The argument is valid for a multiple zero as well: removing one factor is
enough. A root on or outside the radius-`r_n` disk trivially satisfies the
same strict lower bound.

Since `log Q_n=2n log n+c n+O(log n)` and
`log r_n=-2n log n+O(n)`, the displayed two-sided root-distance bounds imply
the claimed leading scale. The derivative also supplies the one-sided
linear refinement

\[
\log\delta_n^Z\le-2n\log n-c n+O(\log n),
\]

but the lower bound does not match its linear coefficient. No two-term root
distance asymptotic is justified by these estimates.

## Y roots: spectral radius and trace

The exact identity

\[
\frac{G_n(\omega,1+\eta)}{P_n(\omega)}
=\det(I+\eta M),\qquad M=H^{-1}V,
\]

shows that a nonzero eigenvalue `lambda` gives the endpoint root
`eta=-1/lambda`. In the present family every eigenvalue is nonzero: the
endpoint-pencil identity gives `det V=A_n^2/omega^(n^2) != 0`, and `H` is
invertible. Therefore `delta_n^Y=1/rho(M)` exactly. This statement uses only
the characteristic polynomial, not diagonalizability.

The trace inequality `|tr M| <= n rho(M)` gives the desired upper distance
bound. The lower distance bound follows from the Gram normalization:
with `S=J^(-1/2)`, the matrix `W=SHS` has norm at most one and
`|det W|=chi_n(2pi)>1/6`. Its smallest singular value therefore exceeds
`1/6`, so `||W^(-1)||<6`. Since `||SVS||<=L_n` and

\[
W^{-1}(SVS)=S^{-1}MS,
\]

the spectral radius satisfies `rho(M)<6L_n`. Thus

\[
\frac1{6L_n}<\delta_n^Y\le\frac n{|\operatorname{tr}M|}.
\]

The separately audited endpoint insertion gives the same leading logarithm
for `|tr M|` as for `Q_n`; and `log L_n=2n log n+O(n)`. This proves the
second root-distance assertion.

## Limits of the conclusion

The proof establishes actual nearby complex zeros. It does not establish
that a nearest `Z` zero lies on the imaginary axis, that there is only one,
that it is simple, or that its displacement equals a Newton correction.
Symmetry of the zero set does not resolve those questions. The endpoint
roots satisfy a different polynomial equation, and are not forced onto the
relation `Y=exp(Z)`.

No rational numerator or denominator is produced, and the argument does
not control the degree and height costs needed for an irrationality-measure
claim about pi.

## Executed checks

The reviewer read and reran `verify_nearest_root_scale.py`; it passed.
It constructs the complete polynomials for `n=2,...,6` and checks their degree,
origin zero, and parity, with 15 independent evaluations against the original
raw determinant. The additional exact checks cover nine factored
logarithmic-derivative examples (including multiplicities), three endpoint
trace/root examples, and 15 Blaschke boundary identities. These are audits of
finite algebraic formulas; the analytic argument above proves the asymptotic
claim.
