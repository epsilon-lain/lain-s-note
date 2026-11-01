# A fixed transverse displacement can destroy the period Gram saving

2026-11-01. These results concern the grouped Hardy Gram matrix. They do not
assert an instability of every compatible polynomial determinant, and do not
decide a Diophantine property of pi.

## 1. A complete-jet kernel subproblem with an exact lower bound

Let the two base nodes be $0,\omega$, where $0<|\omega|<R$. At each node
take the $n$ consecutive functionals

\[
[t^s]f(j\omega+\log(1+t),j\delta),\qquad
j=0,1,\quad 0\le s<n.
\]

Here the transverse variable is held fixed during differentiation, and its
coefficient space contains the types $a=0,1$. In the Hardy coefficient basis
$(a,d)$ the evaluation matrix has entries

\[
Z_{j,s;a,d}=(j\delta)^a
[t^s]\left(\frac{j\omega+\log(1+t)}R\right)^d.
\]

Take exactly these two transverse types, or any larger type set; let
$G_n(\delta)=ZZ^*$. Replacing all logarithmic jets by ordinary normalized
derivatives at each node is a lower triangular row transformation of
determinant one. It leaves every full-row minor, and hence the Gram
determinant, unchanged.

Choose the $2n$ columns with $a=0,1$ and $0\le d<n$. Put

\[
(T_j)_{s,d}=[t^s]\left(\frac{j\omega+t}{R}\right)^d.
\]

Each $T_j$ is upper triangular with diagonal $R^{-s}$. With rows grouped by
node and columns grouped by transverse type, the selected matrix is

\[
Z_{\mathrm{sel}}=
\begin{pmatrix}T_0&0\\T_1&\delta T_1\end{pmatrix},
\qquad
\det Z_{\mathrm{sel}}=\delta^nR^{-n(n-1)}.
\]

Cauchy--Binet, first for finite Hardy truncations and then by convergence,
gives

\[
\det G_n(\delta)\ge
|\delta|^{2n}R^{-2n(n-1)}.                         \tag{1}
\]

At zero displacement, all positive transverse types vanish. The confluent
Cauchy product is

\[
\det G_n(0)=R^{-2n(n-1)}
\left(\frac{|\omega|^2}{R^2-|\omega|^2}\right)^{n^2}.
\]

Consequently

\[
\boxed{\frac{\det G_n(\delta)}{\det G_n(0)}
\ge |\delta|^{2n}
\left(\frac{R^2-|\omega|^2}{|\omega|^2}\right)^{n^2}.} \tag{2}
\]

If $R>\sqrt2|\omega|$ and $\delta\ne0$ is fixed, the logarithm of this ratio
is bounded below by a positive quadratic in $n$ plus a linear term. An
arbitrarily small fixed displacement therefore cannot justify a determinant
comparison with $o(n^2)$ logarithmic loss. At the true period $\omega=2\pi i$,
the previous radius condition $R>14$ already implies $R>\sqrt2|\omega|$.

This is an exact kernel subproblem with complete consecutive time jets and
only transverse derivative index $\beta=0$. As $n$ grows it is not the full
fixed-weight source simplex, whose transverse derivative indices also grow.
It disproves a general transfer principle based solely on small displacement;
it is not by itself a lower bound on the Gram determinant of that full
simplex. Additional transverse types in this same row problem preserve (1),
since they add a positive semidefinite matrix to its positive definite Gram
matrix, and leave the zero-displacement denominator unchanged.

## 2. Conditional full-packet obstruction from arithmetic

There is also a direct full-packet consequence of the earlier arithmetic
bound. Let $M$ be the number of source rows and use the complete type set
$\mathcal A_H=\{a:wa\le H\}$, of cardinality $N_A$. Suppose the rational-center
source matrix has full row rank, so that some raw polynomial minor is
nonzero. Scale its column $Y^hX^\alpha$ by
$Q_\alpha=\prod_i(q_iL_i)^{\alpha_i}$ and apply the determinant-preserving
common Laurent gauge. The established arithmetic inequality is

\[
|\det A_{\rm rat}|\ge\exp(MH b_{\kappa,H}),\qquad
b_{\kappa,H}=
\frac1{MH}\sum_{\rm rows}\sum_i\beta_i\log(q_iL_i).
\]

The factorization is $A_{\rm rat}=Z(\delta)C$, with the same $Z$ as before.
Write

\[
a_R=R/w_0,\qquad
b_R=\max_i\frac{\log(q_iL_i)+\log(2R)}{w_i},\qquad
\mathcal G_R=\max\{b_R,(a_R+b_R)/2\}.
\]

The optimized Laurent-gauge bound gives each component function norm at most
$\exp(H\mathcal G_R+R)$. Hence

\[
\det(C^*C)\le
\left(N_A\exp(2H\mathcal G_R+2R)\right)^M.
\]

Combining the arithmetic inequality with
$|\det A_{\rm rat}|^2\le\det G(\delta)\det(C^*C)$ gives the finite bound

\[
\frac{\log\det G(\delta)}{2MH}
\ge b_{\kappa,H}-\mathcal G_R
-\frac{R+\tfrac12\log N_A}{H}.                    \tag{3}
\]

No common nonzero period minor is needed for this inference. The source
rank assumption alone supplies the rational minor; $G(0)$ is positive
definite independently of that column choice.

Define $S_H=\Psi_{K,R}\sum_\beta n_\beta^2/(MH)$ using the previous exact
Cauchy product. Since
$\log\det G(0)/(2MH)=-S_H+\log R/(2H)$, (3) yields

\[
\boxed{
\frac{\log(\det G(\delta)/\det G(0))}{2MH}
\ge S_H+b_{\kappa,H}-\mathcal G_R
-\frac{R+\tfrac12\log N_A+\tfrac12\log R}{H}.}     \tag{4}
\]

For fixed dimension, weights, data and radius, followed by $H\to\infty$,

\[
\liminf\frac{\log(\det G(\delta)/\det G(0))}{2MH}
\ge \frac{2\Psi_{K,R}}{Kv_0(m+2)}
+\frac{\theta}{m+2}\sum_i\kappa_i-\mathcal G_R.    \tag{5}
\]

In the previously constructed high-dimensional family, set $R=14K$ and take
the transverse weights sufficiently large that $b_R\to1$. Then
$a_R\to0$, $\mathcal G_R\to1$, and the row arithmetic term tends to $1/2$.
Meanwhile

\[
\frac{2\Psi_{K,14K}}{Kv_0(m+2)}
\sim\frac{\phi(14)}{m+2}\left(\frac65\right)^m
\longrightarrow\infty.
\]

Thus that regime cannot admit a structured Gram comparison with
$\log\det G(\delta)-\log\det G(0)=o(MH)$. This remains conditional on the
manuscript's rational-center interpolation/rank hypotheses and the earlier
arithmetic setup. It does not exclude retaining the actual coefficient Gram
matrix, using the compatibility between $Z$ and $C$, or a different bound
for the determinant itself.

## 3. Verification scope

`verify_gram_displacement.py` checks (1)--(2) for $n=1,2,3,4$, base nodes
$0,6i$, displacement $i/100$, radius $20$, and exactly two transverse types.
It independently checks the selected minor, ordinary/logarithmic Gram
determinants, the zero-displacement Cauchy product, and the lower bound using
Gaussian rational arithmetic. The base $6i$ is a rational surrogate for this
kernel calculation; no assertion of $e^{6i}=1$ is made. The asymptotic claims
follow from the displayed proofs, not numerical extrapolation.
