# A complete higher-dimensional raw determinant that factors into the one-dimensional family

2026-11-02. This note constructs a genuine complete multidimensional source packet and a square set of legal raw monomial columns. Its determinant factors exactly into the previously analyzed even one-dimensional period polynomials. Thus the coefficient-content theorem transfers without any new denominator loss. The same identity also identifies a limitation: adding spatial dimensions at fixed depth creates no new root locations and eventually makes the primitive actual-period value large.

The family is sparse and has equal small normal weights. It does not satisfy the original research target's critical half-density or separated logarithmic-weight regime. Nothing here proves or disproves pi's bad approximability.

## 1. Complete source and legal raw columns

Let N>=2 and d>=1 be integers. There are m=d+1 spatial variables X_1,...,X_d,T, where T=X_m, and two centers j=0,1. Use the parameters

\[
H=4N,\qquad w_0=1,\qquad w_i=2/N\ (1\le i\le d),\qquad w_m=2,
\]
\[
\theta=1/(2N),\qquad v_0=2,\qquad
v_i=w_i/\theta=4\ (i\le d),\qquad v_m=w_m/\theta=4N.
\]

In particular v_0>w_0, 2*theta<1 and 2*(w_0/v_0)*theta<1. The complete strict source inequality

\[
v_0s+\sum_{i=1}^{m}v_i\beta_i<H
\]

forces beta_m=0 and becomes s+2|nu|<2N, where nu=(beta_1,...,beta_d) and |nu| is its coordinate sum. Define n_nu=2(N-|nu|). The complete source rows are therefore

\[
L_{N,d}=\{(j,s,(\nu,0)):j=0,1,\ |\nu|<N,\ 0\le s<n_\nu\}.
\tag{1}
\]

Choose the raw columns

\[
\mathcal C_{N,d}
=\{Y^hT^aX_1^{\kappa_1}\cdots X_d^{\kappa_d}:
h=0,1,\ |\kappa|<N,\ 0\le a<2(N-|\kappa|)\}.
\tag{2}
\]

Every column is inside the original weighted column budget:

\[
h+2a+(2/N)|\kappa|
\le1+2(2N-1-2|\kappa|)+(2/N)|\kappa|\le4N-1<H.
\]

Both sets have size

\[
\boxed{M_{N,d}=4\binom{N+d}{d+1}.}
\tag{3}
\]

Indeed, there are binom(t+d-1,d-1) normal multiindices of total degree t, each contributing 4(N-t) rows and columns; summing gives (3).

The ambient column space is much larger. For example, all choices 0<=h<=N, 0<=a<=N and |kappa|<=floor(N^2/2) already lie in its budget. For fixed d these number of order N^(2d+2), whereas M is of order N^(d+1). For fixed N and increasing d, the ambient subspace with h=a=0 and |kappa|<=2N^2 alone has order d^(2N^2), while M has order d^(N-1). Thus the selected-to-ambient density tends to zero in both limits considered below.

## 2. The raw matrix is block upper triangular

For an arbitrary formal or complex Z, let B_(N,d)(Z) have entries

\[
[z^su_1^{\nu_1}\cdots u_d^{\nu_d}u_m^0]
e^{hz}(jZ+z+u_m)^a
\prod_{i=1}^{d}(jZ+z+u_i)^{\kappa_i}.
\tag{4}
\]

At Z=2*pi*i, these formal-period rows agree with the actual exponential endpoint rows. Extracting the normal variables first makes (4) equal to zero unless kappa>=nu coordinatewise. Otherwise it equals

\[
\left(\prod_{i=1}^{d}\binom{\kappa_i}{\nu_i}\right)
[z^s]e^{hz}(jZ+z)^{a+|\kappa-\nu|}.
\tag{5}
\]

Order both row and column normal blocks by increasing total normal degree, with the same ordering inside each degree. A nonzero block then has |kappa|>=|nu|. If their total degrees are equal, coordinatewise inequality forces kappa=nu. The diagonal block is exactly the 2*n_nu by 2*n_nu one-dimensional matrix defining P_(n_nu)(Z), with rows ordered by center then time and columns by frequency then spatial degree.

Consequently the determinant identity is exact, with no unspecified sign or scalar:

\[
\boxed{\mathcal P_{N,d}(Z):=\det B_{N,d}(Z)
=\prod_{|\nu|<N}P_{2(N-|\nu|)}(Z)
=\prod_{k=1}^{N}P_{2k}(Z)^{r_{N,d,k}},}
\tag{6}
\]
\[
r_{N,d,k}=\binom{N-k+d-1}{d-1}.
\]

An equivalent exact column operation replaces each X_i by U_i+T, with U_i=X_i-T. The set (2) is closed under this triangular substitution: a term indexed by lambda<=kappa has new T degree a+|kappa-lambda|, whereas its allowable upper bound increases by 2|kappa-lambda|. The integer transformation has diagonal entries one and inverse of the same form. In adapted columns Y^hT^aU^kappa the matrix is block diagonal. This gives the same identity, but it is not necessary to regard the selected columns as mixed columns: (6) already holds for the actual raw monomials in (2).

## 3. Exact denominator and exact rational-specialization degree

Write A_j=product_(a=0)^(j-1) a!. The preceding even-dimensional denominator theorem gives

\[
d_{2k}=A_{2k}^2/A_k^4,\qquad
R_{2k}=d_{2k}P_{2k}\in\mathbb Z[Z],\qquad
\deg R_{2k}=2k^2,
\]

with leading coefficient (-1)^k and only even powers. Define

\[
\mathfrak d_{N,d}=\prod_{k=1}^{N}d_{2k}^{r_{N,d,k}},\qquad
\mathcal Q_{N,d}=\mathfrak d_{N,d}\mathcal P_{N,d}
=\prod_{k=1}^{N}R_{2k}^{r_{N,d,k}}.
\tag{7}
\]

This is an integer polynomial with leading coefficient plus or minus one. The highest coefficient of mathcal P is exactly its reciprocal up to sign, so mathfrak d is the **least positive integer coefficient denominator**, not merely a sufficient clearing factor. The integer polynomial mathcal Q is primitive.

Its exact degree is

\[
\boxed{D_{N,d}=2\sum_{k=1}^{N}k^2r_{N,d,k}
=4\binom{N+d}{d+2}+2\binom{N+d}{d+1}.}
\tag{8}
\]

For coprime integers a,b with b>=2, the polynomial mathcal Q(ia/b) is real rational. After multiplication by b^D, its numerator modulo every prime dividing b is the unit leading coefficient times i^D*a^D, which is nonzero. Hence

\[
\boxed{\operatorname{den}_{\rm reduced}\mathcal Q_{N,d}(ia/b)=b^{D_{N,d}}.}
\tag{9}
\]

The previous coefficient cancellation therefore transfers exactly, but the value-denominator exponents add rather than disappearing. At fixed d, D is asymptotic to 4*N^(d+2)/(d+2)!.

## 4. Fixed dimension and increasing depth preserve decay

Let omega=2*pi*i. The preceding two-sided Hankel estimate and even least-denominator theorem imply

\[
\log|R_n(\omega)|
=-\tfrac12n^2\log n+
\left[\tfrac34+\log(\pi/\sqrt2)\right]n^2+O(n\log n).
\tag{10}
\]

All factors at omega are nonzero. For fixed d, sum (10) over the multiplicities in (6). Setting n=2k gives

\[
\boxed{\log|\mathcal Q_{N,d}(\omega)|
=-\frac4{(d+2)!}N^{d+2}\log N+O_d(N^{d+2})
\longrightarrow-\infty.}
\tag{11}
\]

For completeness, (8) gives sum 2*k^2*r =4*N^(d+2)/(d+2)!+O_d(N^(d+1)). In the weighted logarithmic term, replacing log k by log N introduces O_d(N^(d+2)): r is at most a constant times N^(d-1), and sum k^2*log(N/k)=O(N^3), for example by comparison with the integrable function x^2*log(1/x). The remainder in (10) sums to O_d(N^(d+1)*log N). These elementary estimates prove (11) without a uniform-in-d asymptotic assumption.

Thus the complete multidimensional construction retains primitive smallness when d is fixed and the time depth grows. This conclusion is about a transcendental argument, not an integer contradiction.

## 5. Fixed depth and increasing dimension do the opposite

Fix N>=2. The highest normal shell |nu|=N-1 contributes P_2(Z)=-Z^2. It has d_2=1, so

\[
|R_2(\omega)|=4\pi^2>1.
\]

Its multiplicity is

\[
r_{N,d,1}=\binom{N+d-2}{N-1}
=\frac{d^{N-1}}{(N-1)!}+O_N(d^{N-2}).
\]

Every other factor has k>=2 and multiplicity O_N(d^(N-2)); its nonzero value at omega is fixed. Equation (6) therefore proves

\[
\boxed{\log|\mathcal Q_{N,d}(\omega)|
=\frac{\log(4\pi^2)}{(N-1)!}d^{N-1}+O_N(d^{N-2})
\longrightarrow+\infty.}
\tag{12}
\]

In fact (8) also gives D_(N,d)=2*d^(N-1)/(N-1)!+O_N(d^(N-2)), so

\[
\frac{\log|\mathcal Q_{N,d}(\omega)|}{D_{N,d}}
\longrightarrow\log(2\pi)>0.
\tag{13}
\]

The mechanism is explicit: increasing the number of normal coordinates at fixed depth predominantly adds copies of the shallowest, non-small block. This limit statement is not interchangeable with (11), whose constants may depend on d.

There is an even stronger geometric limitation. Since all multiplicities r_(N,d,k) are positive for d>=1, (6) gives the equality of root sets

\[
\boxed{\{Z:\mathcal P_{N,d}(Z)=0\}
=\bigcup_{k=1}^{N}\{Z:P_{2k}(Z)=0\},}
\tag{14}
\]

which does not depend on d. Raising the spatial dimension at fixed N changes root multiplicities but creates **no new root location**, and hence cannot bring a root closer to omega than the best one-dimensional factor already did. These are restrictions of this explicit complete block extension, not impossibility results for arbitrary high-dimensional methods.

## Verification and boundary of the result

The companion standard-library verifier checks the complete source, legal raw columns, exact raw matrix determinants against (6), block degrees, minimal coefficient denominators, and the exact rational-value denominator in small cases. The all-parameter identities and the two opposing limits follow from the arguments above, not from finite measurements.

This construction carries the balanced endpoint content into a true multidimensional raw determinant without losing the complete source. It also shows why this particular extension does not reduce the rational specialization cost or exploit new root locations. A successful critical high-dimensional construction would require coupling not present in this block upper triangular family; no assertion about whether such a coupling exists is made here.


## Prior checkpoints

The one-dimensional factors use the [exact even denominator and its asymptotics](../../pi-primitive-period-2026-11-02/denominator-reflection-route/ENDPOINT_CONTENT_AND_DENOMINATORS.md) and the [actual-period two-sided Hankel estimate](../../pi-period-mechanism-2026-11-02/raw-root-mechanism/HANKEL_PERIOD_MAGNITUDE.md). The proof above establishes the new complete multidimensional packet and its raw-matrix factorization explicitly.
