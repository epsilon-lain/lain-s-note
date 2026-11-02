# Roots of the raw period polynomials approach every fixed exponential period

2026-11-02. This note proves an infinite-size root-approximation statement for the specific sparse one-dimensional two-frequency family. For every fixed nonzero integer k, a root of its n-th genuine raw period polynomial lies within O_k(n^-2) of 2*pi*i*k. The proof does not assert that the selected root is imaginary or simple, and does not produce a rational approximation to pi.

## 1. The fixed family and a sharper degree bound

For n>=2 retain the complete source rows (j,s,beta)=(j,s,0), j=0,1 and 0<=s<n, and the raw columns Y^h X^a with h=0,1 and 0<=a<n. For n>=3 these obey the original weighted budget, including the strict condition K*theta<1 with K=2, for

\[
H=2n,\quad w_0=1,\quad w_1=2,\quad v_0=2,\quad v_1=2n,
\quad\theta=1/n.
\]

At n=2 this displayed parameter choice has K*theta=1, so it is only an algebraic boundary test. The same n=2 packet is legal instead with H=5, w_0=1, w_1=2, v_0=5/2, v_1=5 and theta=2/5. This edge case does not affect any subsequent algebraic identity or the large-n convergence theorem.

Rows are ordered by center then time, and columns by frequency then spatial degree. Let

\[
P_n(Z)=\det\bigl[[z^s]e^{hz}(jZ+z)^a\bigr].
\tag{1}
\]

At Z=2*pi*i*k with integer k, this is the actual endpoint determinant because e^(k*2*pi*i)=1. Arbitrary formal Z is not asserted to be an exponential period. Put A_n=prod_(a=0)^(n-1) a! and B_n=prod_(a=0)^(2n-1) a!.

The earlier reflection identity is

\[
P_n(Z)=\det(A+B R)\det(A-B R),
\]

where R=diag((-1)^a) and

\[
A_{sa}=[z^s]e^{z/2}(Z/2+z)^a,\qquad
B_{sa}=[z^s]e^{-z/2}(Z/2+z)^a.
\]

It follows by translating the centers to -Z/2,+Z/2, shifting the local frequency by -1/2, and reflecting the negative-center block. The two parity changes have product sign one.

Write n_0=ceil(n/2) and n_1=floor(n/2). In either factor, column a expands as

\[
\sum_{j=0}^{a}\binom aj(Z/2)^{a-j}v_{j,a\bmod2},
\]

where v_(j,p) is the fixed coefficient vector of z^j(e^(z/2) +/- (-1)^p e^(-z/2)), truncated below degree n. It depends on j and the parity p, but not otherwise on a or Z.

In a surviving determinant term the n_0 even columns must select distinct j values; the n_1 odd columns must also select distinct j values. Otherwise a vector is repeated. Thus the sum of the selected j values is at least binom(n_0,2)+binom(n_1,2). Subtracting from sum a=n(n-1)/2 proves that each factor has degree at most n_0 n_1=floor(n^2/4).

For completeness, let E be the n by n truncated multiplication matrix for e^z, and let T(Z) be polynomial translation f(z)->f(z+Z). The original block determinant reduces to P_n(Z)=det(E T(Z)-T(Z) E). Since T(-Z)=T(Z)^(-1) and det T(Z)=1,

\[
E T(-Z)-T(-Z)E=-T(Z)^{-1}(E T(Z)-T(Z)E)T(Z)^{-1}.
\]

Taking determinants proves P_n(-Z)=(-1)^n P_n(Z). For odd n this lowers the even product-degree bound by one. Hence

\[
\boxed{\deg P_n\le d_n^*:=
\begin{cases}
n^2/2,&n\text{ even},\\
(n^2-3)/2,&n\text{ odd}.
\end{cases}}
\tag{2}
\]

This is only a degree upper bound; exact leading-coefficient formulas are not required below.

## 2. The first nonzero coefficient is explicit

The previously proved formal Hankel identity is

\[
P_n(Z)=\frac{Z^{n^2}}{A_n^2}
\det[\widehat\mu_{i+j}(Z)]_{i,j=0}^{n-1},
\quad \widehat\mu_0=0,\quad
\widehat\mu_j=Z^{-1}-jZ^{-1}\widehat\mu_{j-1}.
\tag{3}
\]

For j>=1 its highest negative-power term is (-1)^(j-1)j! Z^(-j). Factoring Z^(-i) from row i and Z^(-j) from column j shows that the first possible power in P_n is Z^n. Its coefficient is obtained from the factorial Hankel matrix [(i+j)!], with its (0,0) entry replaced by zero, together with an overall sign (-1)^n.

For any integer alpha>=0 and size r,

\[
\det[(i+j+\alpha)!]_{i,j=0}^{r-1}
=\prod_{j=0}^{r-1}j!(j+\alpha)!.
\tag{4}
\]

To prove (4), factor (i+alpha)! from row i. Column j then evaluates the monic degree-j polynomial (i+alpha+1)...(i+alpha+j) at i=0,...,r-1. Its determinant is the Vandermonde product prod j!.

The unmodified determinant in (3) is therefore A_n^2, while its (0,0) cofactor is prod_(j=0)^(n-2) j!(j+2)!=n A_n^2. Replacing that entry by zero subtracts this cofactor. Consequently

\[
\boxed{[Z^n]P_n(Z)=(-1)^{n+1}(n-1),
\qquad \operatorname{ord}_0P_n=n.}
\tag{5}
\]

In particular the polynomial is nonzero for every n>=2, independently of any numerical computation. Define

\[
R_n(Z)=\frac{P_n(Z)}{(-1)^{n+1}(n-1)Z^n},
\qquad R_n(0)=1,
\qquad D_n=d_n^*-n.
\tag{6}
\]

Then deg R_n<=D_n, with D_n=n^2/2-n+O(1). For n>=4 this bound is positive.

## 3. An elementary root-distance lemma

Let R be a polynomial with R(0)=1 and deg R<=D, where D>=1. Let z!=0 and suppose |R(z)|<=v<1. Then R is nonconstant. If delta is the distance from z to its closest root, the factorization R(z)=prod_(rho)(1-z/rho) gives

\[
\left|1-\frac z\rho\right|
=\frac{|z-\rho|}{|\rho|}
\ge\frac{|z-\rho|}{|z|+|z-\rho|}
\ge\frac\delta{|z|+\delta}.
\]

All roots are nonzero because R(0)=1. Since the last ratio lies below one, using the actual degree d<=D yields

\[
|R(z)|\ge\left(\frac\delta{|z|+\delta}\right)^d
\ge\left(\frac\delta{|z|+\delta}\right)^D.
\]

Writing r=v^(1/D)<1, rearrangement proves

\[
\boxed{\delta\le\frac{|z|r}{1-r}.}
\tag{7}
\]

Neither a bound for individual polynomial coefficients nor a root-separation hypothesis is used in this lemma.

## 4. Apply the actual-period Hankel upper bound

Fix a nonzero integer k and put omega_k=2*pi*i*k and t=|omega_k|=2*pi*|k|. The Fourier–Hankel identity from the preceding note gives the upper bound

\[
|P_n(\omega_k)|\le t^{n^2}\frac{A_n^2}{B_n}.
\]

Only the absolute upper bound is used; its proof by the nonnegative Vandermonde-square integral works for every real Fourier parameter t. The special positive cosine lower bound at t=2*pi is not being assumed at the higher periods.

From (6),

\[
|R_n(\omega_k)|\le
v_{n,k}:=\frac{t^{n^2-n}A_n^2}{(n-1)B_n}.
\tag{8}
\]

Whenever v_(n,k)<1, let r_(n,k)=v_(n,k)^(1/D_n). Equation (7) gives the explicit finite conclusion

\[
\boxed{\exists\rho_{n,k}\ne0:\ P_n(\rho_{n,k})=0,
\qquad |\rho_{n,k}-\omega_k|
\le\frac{t r_{n,k}}{1-r_{n,k}}.}
\tag{9}
\]

The left side is in fact positive: P_n has rational coefficients and omega_k is transcendental. Positivity is not needed for the distance upper bound.

This finite bound can be certified using rational arithmetic alone: replace t in (8) and (9) by 44|k|/7, then find any rational q<1 with

\[
q^{D_n}\ge
\frac{(44|k|/7)^{n^2-n}A_n^2}{(n-1)B_n}.
\]

It follows that some root lies within (44|k|/7)q/(1-q) of the actual period. This bound is intentionally much weaker than the earlier narrow root brackets, but it applies by formula at arbitrary size.

## 5. Quadratic convergence to every fixed nonzero period

Elementary factorial summation gives

\[
\ln A_n=\tfrac12 n^2\ln n-\tfrac34n^2+O(n\ln n),\qquad B_n=A_{2n},
\]

and therefore

\[
\ln v_{n,k}
=-n^2\ln n+\left[\frac32+\ln(t/4)\right]n^2+O_k(n\ln n).
\]

Consequently v_(n,k)<1 for every sufficiently large n, and

\[
\ln r_{n,k}
=-2\ln n+3+2\ln(t/4)+O_k(\ln n/n).
\]

Substitution into (9) proves

\[
\boxed{\limsup_{n\to\infty}n^2
\operatorname{dist}(2\pi i k,\{\rho:P_n(\rho)=0,\rho\ne0\})
\le\frac{e^3(2\pi|k|)^3}{16}.}
\tag{10}
\]

Thus a genuine nonzero algebraic root can be selected near each fixed nonzero exponential period for all sufficiently large n. For a fixed finite collection of periods, the resulting shrinking disks are eventually disjoint and the roots in them are distinct.

There is also a uniform extension. The same calculation, with its factorial error independent of k, gives

\[
\ln r_{n,k}=2\ln(t/n)+3-\ln16+O(\ln n/n)
\quad\text{uniformly for }2\pi\le t\le n.
\]

For all sufficiently large n and 1<=|k|<=n/(8*pi), it follows that r_(n,k)<=e^3 t^2/(8n^2)<1/2. Equation (9) consequently yields the uniform, nonoptimal bound

\[
\operatorname{dist}(2\pi i k,\{\rho\ne0:P_n(\rho)=0\})
\le\frac{e^3(2\pi|k|)^3}{4n^2}.
\tag{11}
\]

In particular, if positive integers K_n satisfy K_n=o(n^(2/3)), the maximum distance over 1<=|k|<=K_n tends to zero. The disks are eventually disjoint, so P_n has at least 2K_n distinct nonzero roots in these disks. This is a growing finite collection of approximated periods; no infinite-degree polynomial is assumed.

This proves convergence, but the O(n^-2) rate is a general product bound and is not asserted sharp. The earlier n=4,...,10 imaginary-root certificates approach the first period much faster. No theorem here establishes that a chosen root remains on the imaginary axis, is unique in the disk, is simple, or alternates sides of the period.

These are algebraic roots with an n-dependent degree bound, coming from a sparse one-dimensional family with changing theta. No lower bound on the degrees of the selected roots is proved. Equation (10) is not a construction of rational numbers p/q, does not control their denominators, and does not address pi's bad approximability. Nor does it supply the missing high-dimensional half-density argument.

## Verification

The companion standard-library verifier reconstructs small exact P_n from the formal Hankel identity, checks the original raw block determinant at additional points, confirms (2) and (5), and verifies exact rational instances of the finite root-distance bound. The factorial-Hankel products and reflection degree calculation are checked separately. Infinite convergence follows from the proof above, not a finite root-fitting experiment.
