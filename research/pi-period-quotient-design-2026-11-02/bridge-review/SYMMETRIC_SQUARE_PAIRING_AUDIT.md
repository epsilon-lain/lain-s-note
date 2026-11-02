# Independent audit: period normalization and a cheap symmetric-square invariant

2026-11-02. This checks the arithmetic interface of the supplied S6-period bridge. It uses only its explicitly displayed elliptic subsystem, not the full six-sphere topology. It remains a constructive tool audit and does not revive the old prime-pool exclusion project.

## 1. The bridge's local coefficient statements pass

For f=2F1(1/12,5/12;1;z), the displayed equation gives

\[
(n+1)^2a_{n+1}=(n+1/12)(n+5/12)a_n.
\]

The integer normalization u_n=1728^n a_n satisfies

\[
(n+1)^2u_{n+1}=12(12n+1)(12n+5)u_n.
\]

For odd p>=5, each progression 12k+1 and 12k+5 has at least floor(n/p^e) multiples of p^e among n consecutive indices, cancelling the two factorial valuations. For p=2,3 the progression factors are units, giving v2(u_n)=2s2(n) and v3(u_n)=s3(n). The recurrence and all-index integrality argument are correct.

The logarithmic companion g=f ln z+h has h_1=31/72. Its canonical Wronskian is

\[
fg'-f'g=1/(z\sqrt{1-z}).
\]

This follows directly from Abel's equation and the leading ln z term; it is not a numerical normalization guess.

## 2. Determinant one is not a free pi normalization

The physical modular companion is

\[
\widetilde g=(g-\ln1728\,f)/(2\pi i),\qquad \widetilde g/f=\tau.
\]

Consequently its Wronskian is

\[
f\widetilde g'-f'\widetilde g
=1/(2\pi i\,z\sqrt{1-z}).
\]

A gauge making the canonical Frobenius determinant algebraic or equal to one has absorbed the factor 2pi i. A determinant-one integer transfer matrix preserves a bilinear pairing but does not, by itself, eliminate this connection constant. In particular, at a CM point it is tau that is algebraic; the canonical ratio g/f still contains ln1728+2pi i tau. The bridge correctly uses the CM/Legendre identity to eliminate the additional constants.

This is a scope condition for any simultaneous rational approximation: its initial normalization must be retained in the arithmetic ledger.

## 3. A rational differential invariant for the squared period

The rational-coordinate period in the bridge is

\[
H(v)=\sum_{n\ge0}\binom{2n}{n}^3(v/64)^n=F(v)^2,
\quad F={}_2F_1(1/4,1/4;1;v).
\]

The equality can be proved by substituting F^2 in the third-order equation and checking the normalized analytic solution at v=0. It is also the classical Clausen identity.

The second-order equation for F is

\[
v(1-v)F''+(1-3v/2)F'-F/16=0.
\]

For any twice differentiable y define

\[
\boxed{\mathcal D[y]=2v(1-v)yy''-v(1-v)(y')^2
+(2-3v)yy'-y^2/4.}
\]

Direct substitution proves \(\mathcal D[H]=0\). More generally, let f,g be an independent pair of solutions of the same second-order equation and put y=A f^2+2Bfg+C g^2. Then

\[
\mathcal D[y]=4v(1-v)(AC-B^2)(fg'-f'g)^2.
\]

When the canonical Wronskian is 1/(v sqrt(1-v)), this gives

\[
\boxed{v\mathcal D[y]=4(AC-B^2),}
\]

an exactly conserved discriminant. Thus rank-one period squares can be audited without expanding every product or losing the second-order structure.

The analogous invariant in the original level-one G=f^2 coordinate replaces the final coefficient 1/4 by 5/36. This difference matters when transporting an approximation through the pullback.

There is also a convenient integer normalization for the smaller second-order period. Set w_n=64^n[v^n]F. Then

\[
(n+1)^2w_{n+1}=4(4n+1)^2w_n,
\qquad \sum_{j=0}^n w_jw_{n-j}=\binom{2n}{n}^3.
\]

In fact w_n is the square of the integer 2^n product_(j<n)(4j+1)/n!. At v=1/64 its term denominator is therefore 4096^n before reduction. This proves an actual integer pairing coordinate; it does not make the denominator of the product F(5F+84vF') disappear.

## 4. Polynomial approximants retain a one-term differential defect

Let F_N=sum_(n=0)^N b_n v^n be the ordinary truncation of F. The Gauss recurrence gives exactly

\[
v(1-v)F_N''+(1-3v/2)F_N'-F_N/16
=-(N+1/4)^2 b_N v^N.
\]

Therefore the squared approximant has the explicit rank-one defect

\[
\boxed{\mathcal D[F_N^2]
=-4(N+1/4)^2 b_N v^N F_N^3.}
\]

This is a small structural object compared with a generic scalar Padé defect. It is suitable for exact certification or an approximation scheme that keeps period and derivative together. It does not by itself improve the reduced denominator: the Ramanujan expression becomes F(5F+84vF'), so any proposed rational pair must account for the denominator of this product and its actual gcd.

Enforcing D=0 only at one rational point does not impose a useful arithmetic restriction when y'' remains free: for nonzero rational y and arbitrary rational y', the identity simply solves for a rational y''. Also scaling y preserves D=0, while scaling (5+42theta)y changes its Ramanujan constant. Thus normalization and differential matching must accompany this invariant.

## 5. Source and verification scope

The supplied bridge is at repository commit c4d1ddfc235c7f32ebb901de78349e54c6bd22a9, research/pi-s6-period-bridge-2026-11-02. Its local formulas were checked independently above. The Ramanujan identity is explicitly Formula(2) of Jesus Guillera, *A Method for Proving Ramanujan Series for 1/pi*, arXiv1807.07394v4, https://arxiv.org/pdf/1807.07394 . Its Formula(6) explains the Legendre period pairing and the role of the modular special value. It is not a theorem bounding the reduced denominators of new rational approximants.

The companion standard-library script audits the coefficient formulas, integrality valuations, Clausen square, nonlinear differential invariant and one-term polynomial defect using exact rational polynomial arithmetic. All universal statements have the direct algebraic proofs above; bounded checks do not prove an approximation limit.

## 6. Audit of the transfer into the integer primitive main family

The main construction has an integer endpoint kernel and a rank-two period map to (alpha,R(1)). Let Z be an integer basis of the kernel of this period map, and C be two integer lifts of a basis of its image. All allowed integer data then have the form Zz+Cv with z integral and v in Z^2. Completeness of this splitting must be supplied by Smith normal form and unimodular basis transforms; two arbitrary rational lifts would not suffice.

Let G be the exact positive-definite Hilbert Gram matrix on multiplier coefficients, G_ij=1/(i+j+1). With

\[
A=Z^T GZ,\quad B=Z^T GC,\quad D=C^T GC,
\quad S=D-B^TA^{-1}B,
\]

completion of the square gives

\[
\boxed{\|Zz+Cv\|_G^2
=v^TSv+(z-z_*)^TA(z-z_*),\qquad z_*=-A^{-1}Bv.}
\]

The Schur metric S is positive definite. Exact nearest-plane rounding in the zero-period lattice yields an integer lift with

\[
\|Zz+Cv\|_G^2\le v^TSv+\frac14\sum_j\|b_j^*\|_G^2,
\]

where b_j^* are its Gram--Schmidt vectors. The rounding bill is independent of v, but it is not asserted to be independent of the scale m. Integral LLL may precondition the zero-period basis; it must not be confused with solving the Hilbert closest-vector problem exactly.

For the parity-free kernel r=4/121, Cauchy--Schwarz gives

\[
|I|^2<r^{2m}\|H\|_G^2
\int_0^1(1+x^2)^{-2}dx,
\qquad \int_0^1(1+x^2)^{-2}dx=1/4+\pi/8<9/14.
\]

Consequently the charged quality certificate is

\[
\boxed{Q^2<\frac{72}{7}\frac{q^4}{\alpha^2}
r^{2m}\|H\|_G^2.}
\]

The constant72/7 is correct. This certificate retains the actual reduced denominator.

An important distinction: for a fixed quotient vector v, adding any zero-period vector changes neither I, alpha, R(1), the actual gcd, q, nor the true quality Q. Nearest-plane lifting improves the multiplier's representation cost and the available upper certificate. Only selecting a different quotient vector changes the rational candidate itself. It would be incorrect to describe the zero-period lift alone as creating extra gcd or changing achieved approximation error.

The independent standard-library script verify_hilbert_lifts_independent.py reconstructs the saved Hilbert projections, Gram--Schmidt lengths, nearest-plane integer corrections, exact cost decompositions, unchanged period pairs, reduced fractions and rational quality certificates. It verifies all37 retained candidates across the three declared scales. The generator's Smith/unimodular proof supplies integral completeness; the saved-data audit does not replace that part of the argument.

## 7. Shared-denominator elliptic control reviewed

The separate elliptic-pairing/INTEGER_QUARTER_PERIOD_PAIR.md was also checked. Its smaller integer recurrence, bilinear Ramanujan form, two positive tail bounds, propagated C/E error bound, and exact actual-gcd denominator formula are correct. For degree2s shared-Q approximants matching through3s, its compatibility defect

\[
Q\theta P_0-P_0\theta Q-P_1Q
\]

is divisible by t^(3s+1) and has residual degree at most s-1. The warning about exact compatibility is correct: if P0 and Q are coprime and Q(0) is nonzero, demanding P1/Q=theta(P0/Q) identically forces Q to be constant. A nonconstant paired approximation therefore needs a controlled compatibility defect or a higher derivative denominator. This is a precise condition on that parametrization, not a universal failure theorem for elliptic-period methods.

## 8. Forced period divisibility is a genuine indexed sublattice

The follow-up main construction sets D=odd(lcm(1,...,2m)). All integer endpoint vectors already have D dividing alpha. Let z=(z0,z1) be the R(1) row of an unrestricted integer basis of the rank-two period image, and g=gcd(z0,z1). The odd-prime image theorem implies gcd(g,D)=1.

Choose an integer unimodular U with zU=(g,0), and put

\[
T_D=U\begin{pmatrix}D&0\\0&1\end{pmatrix}.
\]

Then the exact kernel of the congruence zv=0 modulo D is T_D Z^2. Indeed write v=Uw; the congruence is gw0=0 modulo D, equivalent to D dividing w0. Thus its index is exactly D, not D^2. Both resulting periods are D-divisible.

This ensures an actual odd factor D in the period gcd, hence in the raw fraction gcd. The restriction also has a cost: the two lifted period columns become C T_D, and the Schur metric becomes T_D^T S T_D. The implementation retains this transformed metric and the same zero-period lattice; it does not count the gcd saving without the index cost.

The independent verifier's --forced mode checks the exact indices1,3,15, cross-dataset transformed columns and metrics, reconstructed integer lifts, actual period gcds, raw fraction gcds and integer-coordinate contents. All40 retained forced candidates have coordinate content one. In particular the m=2 candidate355/113 has alpha=339, period gcd3, raw fraction gcd768 and coordinate content1. This is not obtained by uniformly multiplying every integer coordinate by3.

For that candidate the new rational Hilbert certificate itself proves

\[
\boxed{Q^2<\frac{217375834098688}{4366350858338469}<\frac1{20}.}
\]

This is a benchmark certificate based solely on the generated integer data and the proved kernel bound. The separate Machin interval gives its much smaller achieved quality, approximately0.00340631. The benchmark fraction is classical; the purpose here is to verify that the congruence restriction, integer lift and rational joint bound work together.

The full unrestricted period image already contains this sublattice, so the forced version is a structured selection restriction, not an enlargement of the original mathematical family. The observed improvement over the unrestricted finite Gauss choices does not imply that every future scale improves. The exact sublattice rule and its charged certificate are the all-parameter results.
