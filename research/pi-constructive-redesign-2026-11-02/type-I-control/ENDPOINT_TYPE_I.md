# Endpoint-coupled type-I Machin construction: derivation and bounded control

2026-11-02. This is a constructive alternative that was implemented and evaluated, not a continuation of the old prime-pool exclusion problem. Its bounded results are poor, so it should remain a quick negative control rather than become the new main research objective.

## Family

For positive integer a, let

\[
F_a(t)=\int_0^{1/a}\frac{dx}{1+tx^2}
=\sum_{k\ge0}\frac{(-1)^kt^k}{(2k+1)a^{2k+1}}.
\]

Choose an exact Machin identity

\[
w_a\arctan(1/a)+w_b\arctan(1/b)=\pi/4,
\quad \gcd(w_a,w_b)=1.
\]

The tested identities are `(a,b;w_a,w_b)=(2,3;1,1)`, `(3,7;2,1)`,
`(2,7;2,-1)`, and `(5,239;4,-1)`. They follow from the tangent addition formula, with the branch fixed by the angles lying in the correct interval.

Choose degrees m0,ma,mb≥0. Seek polynomials P0,Pa,Pb with integer coefficients, of those maximum degrees, satisfying

\[
R(t)=P_0(t)+P_a(t)F_a(t)+P_b(t)F_b(t)=O(t^L),
\]
\[
w_bP_a(1)-w_aP_b(1)=0,
\quad L=m_0+m_a+m_b+1.
\]

The condition at one endpoint is weaker than requiring the two polynomials to be identical. There are M=m0+ma+mb+3 coefficients and M−1 homogeneous rational constraints. Thus a nonzero rational vector always exists, and clearing denominators gives an integer vector. This argument does **not** prove uniqueness or nonzero endpoint values for all degrees. The implementation explicitly tests nullity and rejects Q=0. All tested cases had nullity one and Q≠0.

Because the weights are coprime integers, the endpoint constraint gives integer Q with

\[
P_a(1)=w_aQ,\quad P_b(1)=w_bQ.
\]

Writing A=P0(1), one obtains R(1)=A+Qπ/4. Hence the fully reduced candidate is

\[
r=-4A/Q,\qquad g=\gcd(4A,Q),\quad q=|Q|/g.
\]

## Joint denominator and remainder formula

The exact quality is

\[
\boxed{q^2|\pi-r|=\frac{4|Q|}{g^2}|R(1)|.}
\]

Write P_a(t)=Σ_j c_{a,j}t^j, and likewise for b. Since L exceeds every polynomial degree, cancellation of the first L Taylor coefficients and the geometric series give the exact remainder

\[
\boxed{R(1)=\sum_{s\in\{a,b\}}\int_0^{1/s}
\frac{\sum_j c_{s,j}(-1)^{L-j}x^{2(L-j)}}{1+x^2}\,dx.}
\]

No numerical value of π is needed for this representation or the rigorous upper bound

\[
|R(1)|\le\sum_{s\in\{a,b\}}\sum_j
\frac{|c_{s,j}|}{(2(L-j)+1)s^{2(L-j)+1}}.
\]

Consequently a small formal Taylor remainder is insufficient unless |Q|/g² is accounted for. The code records both quantities, including exact rational bounds. Its descriptive floating values are not used for certificates. The separate exact Machin intervals measure actual quality after reduction.

## Executed finite experiments

`endpoint_typeI.py` tested all degree triples with max degree≤7 and sum≤15: 456 systems, all nullity one and Q≠0. The best fraction was 25/8, with normalized quality about 1.0619298297. The first balanced degree triples give:

| degrees (m0,ma,mb) | reduced candidate | normalized quality |
| --- | --- | ---: |
| (0,0,0) | 10/3 | 1.7256661 |
| (1,1,1) | 8143/2592 | 21.1458075 |
| (2,2,2) | 141446342/45023769 | 3272153.22 |

The additional angle sweep tested 287 degree triples for each identity (max degree≤6, sum≤12):

| identity parameters | best candidate | normalized quality |
| --- | --- | ---: |
| (2,3;1,1) | 25/8 | 1.0619298 |
| (3,7;2,1) | 68/21 | 42.55764 |
| (2,7;2,−1) | 24/7 | 14.06196 |
| (5,239;4,−1) | 3804/1195 | 59507.15 |

Thus merely using smaller arctangent arguments did not produce the required arithmetic improvement in this design. No infinite failure claim is made.

Some endpoint type-I degree choices recover the previous shared-polynomial type-II candidates exactly: degrees (0,1,1) give 311/99, and (1,2,2) give 94408/30051. This is a warning about Padé duality: a different presentation must not automatically be advertised as a genuinely independent source of approximants. No general equivalence theorem is claimed here.

## Where this leaves constructive research

This family supplies a precise test bed and a joint cost formula, but it has not supplied a good new candidate. Keeping more than one free nullspace direction would permit modular coefficient choices, but without an explicit coefficient-cost bound that can simply reparameterize arbitrary rational approximants. The main constructive work should instead require a demonstrated saving in the integral coefficient height, the common divisor, or the rigorous remainder at fixed height before extending a parameter family extensively.

Run `python endpoint_typeI.py` or `python endpoint_angle_sweep.py` in this directory. Python Fraction and SymPy perform the arithmetic. Saved JSON files contain exact rational intervals and coefficient records. These finite experiments neither prove a quality-zero subsequence nor settle π's bad approximability.
