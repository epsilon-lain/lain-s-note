# Linear multiplier period map and its primitive gcd budget

2026-11-02. This is a construction-accounting tool for the user's current family. It does not resume the paused prime-pool problem. It explains which gcd savings can be arranged directly in a two-dimensional multiplier block, and how to avoid simply encoding an already selected rational target.

Fix m>=1 and

\[
A=\binom{2m}{m},\quad C=\binom{4m}{2m},\quad
S=2m-1,\quad R=6m+1,\quad T=2^{6m+2},\quad D=AC.
\]

The previously proved base periods are

\[
\alpha_0=A/64^m,\quad\beta_0=0,
\quad\alpha_1=(R/S)\alpha_0,\quad\beta_1=-1/(SC).
\]

For H(t)=h_0+h_1t with integral h_0,h_1, put v=S h_0+R h_1. When v!=0, exact period extraction gives

\[
\boxed{\frac pq=\frac{T h_1}{D v},\qquad
g=\gcd(T h_1,Dv),\quad q=\frac{D|v|}{g}.}
\tag{1}
\]

All signs can be normalized after taking the gcd. This formula needs no matrix inversion, numerical pi, or continued fraction.

## Primitive coefficient blocks have a fixed gcd budget

Suppose gcd(h_0,h_1)=1. Then

\[
\gcd(h_1,v)=\gcd(h_1,S).
\]

Consequently

\[
\boxed{\operatorname{odd}(g)\mid
\operatorname{odd}\bigl(D\,\gcd(h_1,S)\bigr)
\mid\operatorname{odd}(DS).}
\tag{2}
\]

Indeed, at an odd prime the factor T is a unit, so
v_p(g)=min(v_p(h_1),v_p(D)+v_p(v))
<=v_p(D)+min(v_p(h_1),v_p(v)).
The final minimum is v_p(gcd(h_1,S)). A sharper valid divisor is
gcd(h_1,D)gcd(h_1,S), with odd parts understood.

The full content is bounded by the determinant of the integer period map:

\[
\boxed{g\mid TDS.}
\tag{3}
\]

For (3), the two output coordinates are produced by
the integral matrix [[0,T],[DS,DR]], whose determinant is -TDS.
Its adjugate shows that g divides TDS times both input coordinates.
Their gcd is one, proving (3).

These statements do not say that a useful forced gcd is impossible.
They say exactly what must be charged when it is built into primitive
linear multiplier coefficients at a given scale.

## Exact transport of the quadratic quality

Let rho=J_1/J_0, with J_j=I_m(t^j). Since J_0=alpha_0*pi/4,

\[
I_m(H)=J_0(h_0+h_1\rho),\qquad
\alpha_H=\alpha_0 v/S.
\]

Substituting (1) into the exact quality identity gives

\[
\boxed{Q=\pi S D^2\frac{|v|}{g^2}
|h_0+h_1\rho|.}
\tag{4}
\]

Thus smallness of the multiplier at the weighted integral mean rho,
the integer size of v, and the actual content g must be considered together.
Formula (4) is useful when comparing coefficient-generation rules which
force different divisibility patterns before any quality evaluation.

## Encoding an arbitrary target is not a selection rule

Every rational p/q can be represented in this linear multiplier space.
For example, take

\[
h_1=DS p,\qquad h_0=Tq-RDp.
\]

Then v=STq and (1) reduces exactly to the chosen p/q.
This verifies that mere surjectivity, or the ability to force a large
gcd using a preselected target, supplies no independent approximation
mechanism. A useful rule must select h_0,h_1 from the kernel or other
structural data before the target rational is known.

## A bounded control, not a proposed enlarged search

Four fixed integer coefficient rules were compared through m=16:
sum(-1)^j binom(m,j)^2 t^j;
sum(-1)^j binom(m,j)binom(2m,j)t^j;
the shifted Legendre coefficients;
and shifted Chebyshev multipliers.
Generation used no pi or target rational; a rational Machin interval
was applied only afterwards. None of these 64 controls had Q<1.
This finite result is not a failure theorem and is not a reason to
enlarge this particular scan. The exact map (1)--(4) is the more reusable
output from this block.

The period identities used above have all-scale proofs in the preceding
checkpoint. The attached finite integer audit separately checks (1)--(3)
and arbitrary-target reconstruction on primitive small coefficient blocks.
