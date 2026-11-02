# A small formal-period polynomial on an explicit complex disk

2026-11-02. This note improves the rational-center comparison for the same sparse one-dimensional, two-frequency determinant. The improvement comes from normalizing its true moment matrix by the positive Hilbert Gram matrix before perturbing the endpoint exponential. It gives a bound on a whole disk, avoiding the much larger global derivative estimate. The resulting rational error lower bound has scale exp(-O(q log q)), still much weaker than any fixed power of 1/q.

## 1. Exact endpoint discrepancy

For n>=2 retain the rational polynomial

\[
P_n(Z)=\det\bigl([z^s]e^{hz}(jZ+z)^a\bigr),
\quad j,h\in\{0,1\},\quad 0\le s,a<n,
\qquad A_n=\prod_{a=0}^{n-1}a!.
\]

The parameter legality and the n=2 algebraic boundary case are as in the companion notes. Let omega=2*pi*i, t=2*pi, and define the n-by-n matrices

\[
H_n(Z)_{ij}=\int_0^1u^{i+j}e^{Zu}\,du,
\qquad J_n=H_n(0)=\left[\frac1{i+j+1}\right]_{i,j=0}^{n-1}.
\]

For Z!=0 set

\[
\nu_k(Z)=\sum_{\ell=0}^{k}
\frac{(-1)^\ell k!}{(k-\ell)!Z^{\ell+1}},
\qquad V_n(Z)=[\nu_{i+j}(Z)]_{i,j=0}^{n-1}.
\]

The moments with a formal endpoint variable Y obey mu_0=(Y-1)/Z and mu_k=Y/Z-k*mu_(k-1)/Z. They are affine in Y, and their Y derivative is precisely nu_k. Subtracting their values at Y=e^Z and Y=1, then using the established formal Hankel identity, gives

\[
\boxed{P_n(Z)=\frac{Z^{n^2}}{A_n^2}
\det\bigl(H_n(Z)+(1-e^Z)V_n(Z)\bigr).}
\tag{1}
\]

This identity holds at general nonzero complex Z. At omega its correction term vanishes exactly. It keeps the distinction between a formal-period polynomial and the actual exponential endpoint determinant; differentiating only the latter would omit a derivative of the correction term.

The Hilbert determinant is

\[
\Delta_n:=\det J_n=\frac{A_n^4}{A_{2n}},
\qquad U_n=t^{n^2}\frac{A_n^2}{A_{2n}}.
\tag{2}
\]

## 2. Normalize the positive Gram matrix first

The real symmetric matrix J_n is positive definite, since it is the Gram matrix of 1,u,...,u^(n-1) in L^2[0,1]. Let S_n=J_n^(-1/2) be its real positive square-root inverse. For the operator norm induced by the Euclidean norm, Cauchy--Schwarz in the integral gives

\[
\boxed{\|S_nH_n(Z)S_n\|_2
\le e^{\max(0,\operatorname{Re}Z)}.}
\tag{3}
\]

Indeed, for any unit complex vectors v,w, the associated polynomials v* S_n(1,u,...,u^(n-1))^T and (1,u,...,u^(n-1))S_n w each have L^2 norm one. Their product integrated against e^(Zu) has absolute value at most the right side of (3). No lower bound for det H_n(Z), and no inverse of H_n(Z), is used.

The diagonal entries of the inverse Hilbert matrix are

\[
(J_n^{-1})_{ii}
=(2i+1)\binom{n+i}{n-i-1}^{\!2}\binom{2i}{i}^{\!2}
\le B_n:=2n\,2^{8n}.
\tag{4}
\]

For a direct justification, divide the Cauchy determinant for J_n with row and column i deleted by det J_n. The ratio is
(2i+1)*prod_(k!=i)((i+k+1)/(i-k))^2, which equals the displayed binomial expression. Each binomial coefficient is at most 2 raised to its top index, giving the safe bound B_n. Since J_n^(-1) is positive definite, its norm is at most its trace, hence at most n B_n.

When 5<=|Z|<=8, reversing the finite sum in nu_k gives

\[
|\nu_k(Z)|
\le k!|Z|^{-k-1}\sum_{j=0}^k\frac{|Z|^j}{j!}
\le3^8(2n-2)!\qquad(0\le k\le2n-2).
\]

Here e<3 and k!<=(2n-2)! suffice. A matrix whose entries have magnitude at most M has operator norm at most nM. Consequently

\[
\boxed{\|S_nV_n(Z)S_n\|_2\le L_n,
\quad L_n:=n^2B_n3^8(2n-2)!
=2\cdot3^8n^3 2^{8n}(2n-2)!.}
\tag{5}
\]

The bound is explicit and deliberately generous. In particular L_n>=n.

## 3. A uniform disk bound

Put

\[
\boxed{r_n=\frac1{8nL_n}.}
\tag{6}
\]

For |Z-omega|<=r_n we have 5<|Z|<8, using 6<t<44/7 and r_n<=1/(8n^2)<=1/32. Also Re Z<=r_n and

\[
|1-e^Z|=|1-e^{Z-\omega}|
\le|Z-\omega|e^{|Z-\omega|}\le2r_n.
\]

Factoring J_n^(1/2) from both sides of the matrix in (1), and using |det M|<=||M||_2^n, equations (3)--(5) imply

\[
\frac{|P_n(Z)|}{U_n}
\le\left(\frac{|Z|}{t}\right)^{n^2}
\bigl(e^{r_n}+2r_nL_n\bigr)^n
\le\exp\left(\frac{n^2r_n}{t}+nr_n+2nr_nL_n\right).
\]

The last step uses 1+x<=e^x and e^r+2rL<=e^(r+2rL). Since t>6 and L_n>=n, the three exponents are bounded by 1/48, 1/8, and 1/4, respectively. Their sum is 19/48<1/2, and e^(1/2)<2. Thus

\[
\boxed{|P_n(Z)|<2U_n
\qquad\text{throughout }|Z-2\pi i|\le r_n.}
\tag{7}
\]

The proof uses the positive Gram normalization and elementary matrix norms, not a Taylor extrapolation from the value at the center. In particular it controls all higher perturbation orders at once. The earlier positive cosine lower bound is not needed for (7).

An optional derivative consequence follows from the Cauchy estimate: on the smaller disk |Z-omega|<=r_n/2,

\[
|P_n'(Z)|\le4U_n/r_n.
\]

Since log(1/r_n)=2n log n+O(n), this has size U_n exp(O(n log n)). The disk shrinks at that same scale; the estimate must not be used on a fixed-radius disk.

## 4. Apply the exact rational-value denominator

Now assume n>=4 is even. The endpoint-content theorem gives the exact least coefficient denominator

\[
d_n=\frac{A_n^2}{A_{n/2}^4},
\qquad R_n=d_nP_n\in\mathbb Z[Z],
\qquad N_n=\deg R_n=n^2/2.
\]

The polynomial R_n is even and has leading coefficient (-1)^(n/2). For a reduced rational p/q, write 2p/q=a/b in lowest terms. If b>=2, the rational-specialization theorem gives

\[
\operatorname{den}_{\rm reduced}R_n(ia/b)=b^{N_n},
\qquad |R_n(ia/b)|\ge b^{-N_n}.
\]

Choose an even n satisfying

\[
\boxed{2d_nU_n<b^{-N_n}.}
\tag{8}
\]

If ia/b were in the closed radius-r_n disk around omega, (7) would contradict this rational nonzero-value floor. Therefore, without a prior proximity assumption,

\[
\boxed{\left|\pi-\frac pq\right|>
\frac{r_n}{2}
=\frac1{32\cdot3^8 n^4 2^{8n}(2n-2)!}.}
\tag{9}
\]

Condition (8) is strictly certified by replacing 2*pi in U_n by 44/7 and comparing rational numbers. Such a search terminates for each fixed b, by the established denominator asymptotic. The b=1 cases involve only q=1,2 and may be handled separately.

## 5. The resulting asymptotic scale

Let

\[
\kappa=e^{3/2}\pi^2/2.
\]

The previous asymptotic gives log(d_nU_n)/N_n=-log n+log kappa+o(1). For any fixed c>kappa, take n to be the least even integer at least cb. Then (8) holds for all sufficiently large b. Stirling's estimate in (9) gives

\[
\boxed{\log(r_n/2)=-2cb\log b+O_c(b).}
\tag{10}
\]

Thus this comparison supplies an exp(-O(q log q)) lower bound. It improves the exp(-O(q^2 log q)) scale obtained with the earlier global derivative bound for this construction. It is still much weaker than a fixed multiple of q^(-2), and even than any fixed inverse power of q. It is not claimed as a new best general irrationality bound for pi.

The rational certificate based on 44/7 has its own slightly larger asymptotic threshold e^(3/2)(22/7)^2/2. Equation (10) concerns the analytic condition (8); the finite verification uses the rational sufficient condition and does not identify those two constants.

No polynomial factorization or critical high-dimensional construction is established in this note. It identifies a genuine benefit from treating the whole endpoint perturbation in its natural Gram geometry, while retaining the remaining value-denominator cost.

## Verification

The standard-library companion checks the formal two-parameter raw determinant against the Hankel recurrence at rational Z,Y, the endpoint-derivative formula, the inverse-Hilbert diagonal identity and bounds, and exact rational sufficient heights and radii. The all-dimension norm argument and disk estimate are proved above; finite checks do not extrapolate them.


## Prior checkpoints

The formal Hankel identity and the Hilbert determinant are proved in the [earlier moment analysis](../../pi-period-mechanism-2026-11-02/raw-root-mechanism/HANKEL_PERIOD_MAGNITUDE.md). The exact even coefficient denominator comes from the [endpoint-content theorem](../../pi-primitive-period-2026-11-02/denominator-reflection-route/ENDPOINT_CONTENT_AND_DENOMINATORS.md); the rational-value denominator and previous comparison are in [rational specialization](../../pi-primitive-period-2026-11-02/rational-comparison-next/RATIONAL_SPECIALIZATION_COST.md). The present norm argument does not use the earlier positive cosine lower bound.
