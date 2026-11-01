# A grouped Hardy-space determinant bound at the true period

2026-11-01. This note continues the exact removal of logarithmic tails. It gives an upper bound for the **whole constant-translation sum**, retaining correlations before taking absolute values. The finite Hilbert-space identities and kernel formulas are unconditional. Nonzero high-dimensional original minors still require a separate rank input. No conclusion about π's bad approximability is obtained.

## 1. A factorization that groups every remaining constant shift

Fix the original full row packet

\[
\rho=(j,s,\beta),\qquad 0\le j<K,\quad
v_0s+w\beta/\theta<H,
\]

and a fixed family of $M$ polynomial columns from $P_W(H)$. All weights, centers and the dimension are fixed before $H\to\infty$. Put $\omega=2\pi i$, $\delta_i=r_i-\omega$, and $\epsilon_{j,i}=j\delta_i$. Let $\mathcal A_H=\{a\ge0:wa\le H\}$ and $N_A=|\mathcal A_H|$. One may replace this by any smaller finite type set containing every nonzero column coefficient and every source type used below.

For a column $P$, write

\[
f_{a,P}(z)=[u^a]P(e^z,z+u_1,\ldots,z+u_m)
=\sum_{d\ge0}c_{a,d,P}(z/R)^d,
\qquad R>7K.
\]

The exact tail-free period translation gives a matrix product $A=ZC$, where the intermediate labels are all pairs $(a,d)$ and

\[
Z_{(j,s,\beta),(a,d)}
=\binom a\beta\epsilon_j^{a-\beta}
[t^s]\left(\frac{j\omega+\ln(1+t)}R\right)^d,
\quad C_{(a,d),P}=c_{a,d,P}.                 \tag{1}
\]

An entry is zero when $a\not\ge\beta$. This is the actual-period identity: the substitution $e^{hj\omega}=1$ has been made only at $\omega=2\pi i$. The remaining constant shifts are summed in the matrix multiplication, not estimated individually.

Take the intermediate space to be the Hilbert direct sum of one Hardy coefficient space $\ell^2(\mathbb Z_{\ge0})$ for each $a$. Each row of $Z$ is square summable, because its test points lie strictly inside the disk of radius $R$. Every column of $C$ is square summable by Parseval. Cauchy–Binet followed by Cauchy–Schwarz on exterior powers gives

\[
\boxed{|\det A|^2\le \det(ZZ^*)\det(C^*C).} \tag{2}
\]

For completeness, first restrict $d$ to a finite range. Finite Cauchy–Binet writes the determinant as the scalar product of two vectors of maximal minors; the squared norms of those vectors are the two Gram determinants. Square summability lets the finite matrices and their Gram matrices converge, yielding (2).

If every $f_{a,P}$ has modulus at most $D_H$ on $|z|=R$, Parseval and Hadamard give

\[
\det(C^*C)\le(N_A D_H^2)^M.
\]

Consequently, with $G(\delta)=ZZ^*$,

\[
\boxed{
\frac{\ln|\det A|}{MH}
\le\frac{\ln\det G(\delta)}{2MH}
+\frac{\ln D_H}{H}+\frac{\ln N_A}{2H}.}     \tag{3}
\]

For the original monomial columns one may use

\[
\frac{\ln D_H}{H}
=\max\{R/w_0,\ln(2R)/w_*\},\qquad w_*=\min_iw_i.
\]

The type-count cost is $o(1)$ because $N_A=O(H^m)$ at fixed parameters. Equation (2), using the actual coefficient Gram matrix instead of its norm relaxation, is often much sharper still. The formulas retain the exact dependence on all shifts and the row correlations represented by $Z$. They do not preserve every cancellation in $\det(ZC)$: Cauchy–Schwarz discards relative phases between the two exterior-power vectors.

## 2. An explicit reproducing-kernel matrix

Define the finite transverse polynomial

\[
W_{j\beta,k\gamma}(\delta)
=\sum_{\substack{a\in\mathcal A_H\\a\ge\beta,\ a\ge\gamma}}
\binom a\beta\binom a\gamma
\epsilon_j^{a-\beta}\overline{\epsilon_k}^{a-\gamma}.
\]

Since the coefficients of $\ln(1+t)$ are real, the full Gram entry is

\[
\boxed{
G_{(j,s,\beta),(k,\ell,\gamma)}
=W_{j\beta,k\gamma}(\delta)
[t^su^\ell]
\frac1{1-(j\omega+\ln(1+t))
              (-k\omega+\ln(1+u))/R^2}.}    \tag{4}
\]

This is a finite matrix whose entries are rational functions of the period and polynomial functions of the shifts and their conjugates. The infinite Taylor sum has been performed exactly. For arbitrary complex nodes, replace $j\omega$ and $-k\omega$ by a node and the conjugate of the other node.

The denominator is a Szegő reproducing kernel. Formula (4) is more than a determinant of rowwise norm bounds: its off-diagonal entries retain the geometric overlap of the rows. It is therefore able to detect cancellation and near dependence that Hadamard's product of row norms misses.

## 3. Exact formula when the approximation error is zero

At $\delta=0$, the transverse factor becomes $\mathbf1_{\beta=\gamma}$. Thus $G(0)$ is a direct sum indexed by $\beta$. For each source type define

\[
n_\beta=\left\lceil\frac{H-w\beta/\theta}{v_0}\right\rceil>0.
\]

The corresponding block has the $n_\beta$ consecutive time tests $s=0,\ldots,n_\beta-1$ at each of $K$ centers. Put $z_j=j\omega/R$. If $n=n_\beta$, then

\[
\boxed{
\det G_\beta(0)
=R^{-Kn(n-1)}
\frac{\displaystyle\prod_{0\le j<k<K}|z_j-z_k|^{2n^2}}
{\displaystyle\prod_{j,k=0}^{K-1}(1-z_j\overline z_k)^{n^2}}.} \tag{5}
\]

The product in the denominator is positive. For the purely imaginary equally spaced nodes used here every individual factor is positive as well.

To prove (5), first replace the test $[t^s]f(j\omega+\ln(1+t))$ by the ordinary normalized derivative test $f^{(s)}(j\omega)/s!$. On a complete consecutive jet block the change is triangular with diagonal one, because $\ln(1+t)=t+O(t^2)$. It therefore does not change a Gram determinant. Rescaling from the physical variable to the unit disk contributes $R^{-2K\sum_{s=0}^{n-1}s}=R^{-Kn(n-1)}$.

For distinct unit-disk nodes $x_i$, the Cauchy determinant identity is

\[
\det\left(\frac1{1-x_i\overline x_j}\right)
=\frac{\prod_{i<j}|x_j-x_i|^2}
       {\prod_{i,j}(1-x_i\overline x_j)}.
\]

It follows by multiplying both sides by the denominator: the numerator is alternating separately in the nodes and their conjugates, has exactly the Vandermonde degrees, and its leading coefficient is one. Replace each node $z_j$ by $n$ nearby nodes, divide their evaluation rows and columns by their within-cluster Vandermonde factors, and take the limit. Divided differences converge to normalized derivatives. Within-cluster numerator factors disappear against these normalizations; cross-cluster factors acquire exponent $n^2$, as do the denominator factors. This proves (5), with no extra derivative factorials.

In particular, (5) is strictly positive for every finite full packet. Thus the grouped Gram norm has a nonzero constant term at zero approximation error. This does **not** make it uniformly bounded away from zero as $H$ grows, and it does not say that the actual determinant is nonzero for arbitrary columns.

## 4. A sharper period-determinant exponent

Write $\Omega=|\omega|=2\pi$ and define

\[
\Psi_{K,R}
=\frac{K^2}{2}\ln R
-\sum_{0\le j<k<K}\ln(\Omega(k-j))
+\frac12\sum_{j,k=0}^{K-1}
       \ln\left(1-\frac{\Omega^2jk}{R^2}\right).          \tag{6}
\]

Taking half the logarithm of (5), then summing over types, gives the exact identity

\[
\frac12\ln\det G(0)
=-\Psi_{K,R}\sum_\beta n_\beta^2+\frac M2\ln R.          \tag{7}
\]

The original simplex count gives

\[
\frac{\sum_\beta n_\beta^2}{MH}
\longrightarrow\frac2{K v_0(m+2)}.
\]

Therefore the exact-period full-row matrix $B=A|_{\delta=0}$ satisfies

\[
\boxed{
\limsup_{H\to\infty}\frac{\ln|\det B|}{MH}
\le\max\{R/w_0,\ln(2R)/w_*\}
-\frac{2\Psi_{K,R}}{K v_0(m+2)}.}           \tag{8}
\]

There is no separate $\ln2/v_0$ time-extraction charge here: the complete consecutive time jets were transformed exactly. This refinement concerns the full identity selection at the period, not arbitrary sparse selections of test time orders.

If columns are pre-scaled by at most $e^{\kappa H}$ and then mixed with coefficients in $\{0,1,2\}$, the only change in this bound is the added $\kappa+\ln(2N)/H$ coefficient-norm cost. The Gram matrix itself is unchanged. Its arithmetic counterpart must still be tracked separately.

### The large-center-number coefficient is explicit

For fixed $r>7$, put $R=rK$. Elementary Riemann sums in (6) give

\[
\frac{\Psi_{K,rK}}{K^2}\longrightarrow
\phi(r):=\frac12\ln\frac r\Omega+\frac34
+\frac12\int_0^1\!\int_0^1
\ln\left(1-\frac{\Omega^2}{r^2}xy\right)dx\,dy.
\tag{9}
\]

Indeed, the weighted node differences have limiting logarithmic integral
$\int_0^1(1-x)\ln x\,dx=-3/4$; its integrable endpoint is handled by first restricting to $x\ge\varepsilon$ and then letting $\varepsilon\downarrow0$. The double-kernel integral has a uniformly convergent series. With $a=(\Omega/r)^2<1$,

\[
\phi(r)=\frac12\ln\frac r\Omega+\frac34
-\frac12\sum_{k\ge1}\frac{a^k}{k(k+1)^2}.
\tag{10}
\]

At $r=14$, this coefficient is rigorously greater than $11/10$. One simple rational certificate uses $\pi<22/7$, hence $14/\Omega>49/22$ and $a<484/2401$. The three-term positive series with $t=27/71$ gives

\[
\ln(49/22)>2(t+t^3/3+t^5/5)>4/5.
\]

Also $\sum_{k\ge1}a^k/[k(k+1)^2]\le a/[4(1-a)]$. Therefore

\[
\phi(14)>\frac{23}{20}-\frac{121}{3834}
=\frac{42881}{38340}>\frac{11}{10}.
\]

For comparison, the previous uniform collision coefficient at $R=14K$ is $\ln2/2$. This stronger coefficient and the absent time-extraction charge apply to the exact-period **full consecutive-jet packet**. They cannot be substituted unchanged for nonzero errors or arbitrary sparse test packets.

## 5. A nonzero-shift exact certificate

Take a complete source-format packet with $m=1$, $K=2$, $H=4$, $v_0=2$, $w_0=w_1=1$, and $\theta=1/4$. It consists of $s=0,1$, $\beta=0$ at both centers. Use columns $1,X,X^2,X^3$, and choose radius $R=20$.

For an exact finite kernel test set the base centers to $0,6i$ and their shifts to $0,i/100$. These rational imaginary nodes are a surrogate for checking the formulas, **not** the actual period. The columns have no $Y$ factor, so the certificate never asserts $e^{6i}=1$. The second actual center is $601i/100$ and the exact determinant is

\[
\det A=(601/100)^4=130466162401/10^8.
\]

The type set $a=0,1,2,3$ contains the entire transverse support of these columns. Direct exact calculation gives

\[
\det G=\frac{97054745474616020001}
{109719937600000000000000000000},
\qquad
\det(C^*C)=4240876297445601.
\]

The squared grouped upper bound divided by the squared actual determinant is

\[
\frac{411597169637914786036953189425465601}
     {186758908887541254503214602841760000}
<\frac94.
\]

Thus this fully grouped bound is less than $3/2$ times the actual determinant modulus. In contrast, replacing $\det G$ by the product of its diagonal entries makes the squared bound over ten thousand times larger. This is a concrete improvement of more than a factor one hundred in the unsquared bound, obtained by retaining row correlations at a **nonzero** shift. It is a comparison with the rowwise Hardy norm relaxation, not a claimed asymptotic improvement of the manuscript's final endpoint criterion.

## 6. What remains unresolved

The formula supplies an explicit finite target: estimate $\det G(\delta)$ for the actual error vector, or keep more of $C^*C$ instead of its uniform column norm. The exact zero-error product cannot be substituted for the nonzero-error Gram determinant without a quantitative stability argument that is uniform in the required $H$ range. Small shifts alone do not give that stability.

At fixed $H$, $G(0)$ is positive definite, so no positive power of the approximation error is forced by this grouped norm. At growing $H$, the Gram matrix can be extremely ill conditioned. Moreover, replacing $C^*C$ by a generic norm forgets algebraic relations between the functions indexed by different transverse types. Those limitations must not be turned into a statement that further cancellation is impossible.

Run `python3 verify_grouped_gram.py`. The standard-library verifier checks nine confluent Gram determinants, 196 independently computed kernel coefficients, four full nonzero/zero-shift certificates, Hermitian symmetry and positive principal minors, exact matrix factorization for the polynomial columns, the Gram determinant inequalities, and the rational lower certificate for $\phi(14)$. All finite calculations use Gaussian rational arithmetic. They do not certify source weight thresholds or any approximation property of π.
