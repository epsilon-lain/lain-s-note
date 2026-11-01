# Small escape fractions preserve the full-period Hardy saving

2026-11-01. This note gives a uniform upper bound for actual-period minors whose row sets differ from the complete source packet in k rows. It uses the conditioning of a full Hardy jet Gram matrix, not the inverse of a selected polynomial interpolation matrix. It supplies neither a lower bound for a period minor nor a conclusion about pi's bad approximability.

## 1. The finite row-exchange inequality

Let L be M linearly independent vectors in a Hilbert space, let J be another set of M vectors, and suppose their intersection I has M-k vectors. Write G_L and G_J for their Gram matrices. Let D=L\I. The Schur complement and Hadamard's inequality give

\[
\det G_J\le \det G_I\prod_{r\in J\setminus L}\|r\|^2.
\]

Jacobi's complementary-minor identity gives

\[
\det G_I=\det G_L\det((G_L^{-1})_{D,D}).
\]

Consequently

\[
\boxed{\frac{\det G_J}{\det G_L}
\le \lambda_{\min}(G_L)^{-k}
\prod_{r\in J\setminus L}\|r\|^2.}             \tag{1}
\]

The empty Gram determinant is one. If G_J is singular, (1) is immediate; otherwise the Schur complement proof applies directly. This does not require J to preserve bucket counts, so it also applies to every admissible transport selection.

## 2. The Hardy jet space at the actual period

Put omega=2*pi*i, Omega=|omega|, K>=1, and R=14K. For every transverse type a, take an independent copy of the Hardy coefficient space with norm

\[
\left\|\sum_{d\ge0} c_d(z/R)^d\right\|^2=\sum_{d\ge0}|c_d|^2.
\]

The row (j,s,a), 0<=j<K, is the functional extracting the ordinary normalized Taylor coefficient [t^s]f_a(j omega+t). Its Gram kernel is

\[
G_{(j,s,a),(l,t,b)}
=\mathbf1_{a=b}[x^s y^t]
\frac1{1-(j\omega+x)(\overline{l\omega}+y)/R^2}.
\]

Use the complete source packet

\[
L_H=\{(j,s,\beta):0\le j<K,\quad
v_0s+\sum_i v_i\beta_i<H\},
\qquad n_* =\lceil H/v_0\rceil.
\]

For each source type beta, the time orders at each center are precisely 0,...,n_beta-1, where n_beta=ceil((H-v beta)/v_0)<=n_*. Hence G_L is a direct sum of complete K-node jet Gram blocks. It is strictly positive definite independently of any finite polynomial column space.

Every target row has bounded norm uniformly even in its time order. Set A=(K-1)Omega and r=(R-A)/2. The inequalities 6<Omega<44/7 give r>1 and

\[
\frac{A+r}{R}=\frac{R+A}{2R}<\frac{71}{98}.
\]

A unit Hardy function obeys |f(z)|<=(1-|z|^2/R^2)^(-1/2). Cauchy's coefficient estimate on |t|=r therefore gives

\[
\|[t^s]f(j\omega+t)\|^2
\le r^{-2s}\left(1-\frac{(A+r)^2}{R^2}\right)^{-1}
\le \boxed{\frac{9604}{4563}}.                \tag{2}
\]

This uses ordinary Taylor coefficients. Applying a separate logarithmic-jet change on a sparse target set need not preserve its determinant and is not part of the argument.

## 3. An explicit exponential inverse bound for complete jet blocks

Define the rational constant

\[
L_K=\frac{142K-44}{7},\qquad
Q_K=L_K\frac{(K+1)^2}{4}
\frac{(L_K/3)^{K-1}}{(K-1)!}.                 \tag{3}
\]

Then every complete block with n consecutive jets at each of K centers satisfies

\[
\boxed{\lambda_{\min}(G_{K,n})^{-1}
\le K n^3 Q_K^{2n}.}                         \tag{4}
\]

Here is a direct interpolation proof. Write x_j=j omega and

\[
\Lambda_j(z)=\prod_{l\ne j}
\left(\frac{z-x_l}{x_j-x_l}\right)^n.
\]

Let d_{j,l} be the Taylor coefficients of 1/Lambda_j(x_j+t), and set

\[
p_{j,s}(z)=(z-x_j)^s\Lambda_j(z)
\sum_{l=0}^{n-1-s}d_{j,l}(z-x_j)^l.
\]

These are Hermite cardinal polynomials of degree at most Kn-1:

\[
[t^r]p_{j,s}(x_l+t)=\mathbf1_{j=l}\mathbf1_{s=r}
\quad(0\le r<n).
\]

On the circle |t|=Omega/2, the reciprocal factor is bounded by C_j^n, where

\[
C_j=\prod_{d=1}^{j}\left(1-\frac1{2d}\right)^{-1}
\prod_{d=1}^{K-1-j}\left(1-\frac1{2d}\right)^{-1}
\le(j+1)(K-j)\le\frac{(K+1)^2}{4}.
\]

The first inequality follows term by term from 2d/(2d-1)<=(d+1)/d. Thus |d_{j,l}|<=C_j^n(2/Omega)^l. On |z|=R, all |z-x_l| are at most R+A<L_K, while

\[
\prod_{l\ne j}|x_j-x_l|
=\Omega^{K-1}j!(K-1-j)!.
\]

Since L_K>1 and 2/Omega<1, summing at most n terms yields

\[
\|p_{j,s}\|_{H^2}\le\sup_{|z|=R}|p_{j,s}(z)|
\le n\left[
\frac{L_K^K C_j}{\Omega^{K-1}j!(K-1-j)!}
\right]^n
\le n Q_K^n.
\]

The last step uses Omega>6 and j!(K-1-j)!>=(K-1)!/2^(K-1). When K=1 the construction is simply p_{0,s}=z^s; the same loose bound holds with Q_1=14.

Let Z map a Hardy function to its Kn jet coefficients, and let P map prescribed jet coefficients to the above cardinal interpolant. Then ZP=I. The minimum-norm interpolation map is Z^*(ZZ^*)^(-1), so

\[
G_{K,n}^{-1}\preceq P^*P,
\quad\lambda_{\max}(G_{K,n}^{-1})
\le\operatorname{tr}(P^*P)
\le Kn\,n^2Q_K^{2n}.
\]

This proves (4). Applying it to every transverse block gives

\[
\lambda_{\min}(G_{L_H})^{-1}\le K n_*^3Q_K^{2n_*}.
\tag{5}
\]

There is no factor proportional to the number of transverse types in (5), because the inverse is block diagonal and its operator norm is the maximum of the block norms.

## 4. Uniform control of every small escape fraction

Combine (1), (2) and (5). For every M_H-row target selection J with k=|J\L_H|,

\[
\boxed{\frac{\det G_J}{\det G_{L_H}}
\le\left[\frac{9604}{4563}\,K n_*^3Q_K^{2n_*}\right]^k.}
\tag{6}
\]

Equivalently, with rho_H=k/M_H,

\[
\frac{\ln\det G_J}{2M_HH}
\le\frac{\ln\det G_{L_H}}{2M_HH}
+\rho_H\left[
\frac{n_*}{H}\ln Q_K+
\frac{\ln(9604K/4563)+3\ln n_*}{2H}
\right].                                      \tag{7}
\]

All parameters are fixed before H tends to infinity. The bracket tends to ln(Q_K)/v_0, uniformly in J. In particular, **every sequence with k=o(M_H) loses only o(M_H H) in the logarithm of the Hardy determinant**. A fixed tiny displacement can still move a positive fraction of rows in the full expansion; this statement does not contradict the earlier displacement-instability example.

The constant is explicit and grows only linearly in K to leading order. For K>=2,

\[
\ln Q_K\le 3K+3\ln K+2.                     \tag{8}
\]

For example, use N!>=(N/e)^N, L_K<21K, K+1<=2K, ln(7)<2, ln(21)<4, and (K-1)ln(K/(K-1))<1. For K=1, ln Q_1=ln14<3. Formula (8) is a deliberately loose bound; (3) is sharper.

## 5. Consequence for actual-period polynomial minors

Let B_J be the actual-period matrix on a fixed M_H-element raw column family in the weighted polynomial space w_0 h+w alpha<=H. The factorization B_J=Z_J C uses

\[
f_{a,P}(z)=[u^a]P(e^z,z+u_1,\ldots,z+u_m).
\]

At z=j omega, the factor e^(h j omega) is exactly one. Thus the Z_J above is the appropriate actual-period row matrix. The number N_A of permitted types is O(H^m), and for each raw monomial column and every a,

\[
\|f_{a,P}\|_{H^2(R)}\le
\exp\{H g_R\},\qquad
g_R=\max\{R/w_0,\ln(2R)/w_*\},\quad w_* =\min_i w_i.
\]

Exterior Cauchy--Schwarz gives |det B_J|^2<=det G_J det(C^*C), and Hadamard bounds det(C^*C)<=(N_A exp(2H g_R))^M_H. The exact source Gram product already established in [GROUPED_HARDY_GRAM.md](../pi-period-structure-2026-11-01/gram/GROUPED_HARDY_GRAM.md) gives

\[
\frac{\ln\det G_{L_H}}{2M_HH}
\longrightarrow -S_R,\qquad
S_R=\frac{2\Psi_{K,R}}{K v_0(m+2)}.
\]

For completeness, the fixed-node constant is

\[
\Psi_{K,R}=\frac{K^2}{2}\ln R
-\sum_{0\le j<l<K}\ln(\Omega(l-j))
+\frac12\sum_{j,l=0}^{K-1}
\ln\left(1-\frac{\Omega^2jl}{R^2}\right).
\]

All logarithms have positive arguments. Its complete-jet formula follows from the confluent Cauchy determinant; no estimate for a sparse jet determinant is substituted into this identity.

Therefore, uniformly over admissible target selections,

\[
\boxed{\frac{\ln|\det B_J|}{M_HH}
\le g_R-S_R+\rho_H\frac{\ln Q_K}{v_0}+o(1).}
\tag{9}
\]

Zero minors satisfy this with ln0=-infinity. The bound extends to commonly scaled or mixed columns only after adding their actual coefficient-norm cost; none is included silently here. It concerns an absolute upper bound, not the ratio of an individual polynomial minor to the source minor.

In particular, the previously known full-period saving survives unchanged at leading order whenever rho_H tends to zero. If a specified arithmetic comparison requires a threshold L_ar and g_R-S_R<L_ar, then this upper-bound route also covers every fixed rho below

\[
\frac{v_0(L_{\rm ar}-g_R+S_R)}{\ln Q_K},
\]

provided all other column-scaling and transport costs in that comparison are included separately. This is a conditional protected interval, not a claim that the threshold inequality holds in the intended pi problem. Transport coefficients and the remainder at larger escape fractions must still be controlled.

## Verification scope

The companion exact verifier checks Hermite interpolation, the positive-semidefinite difference P^*P-G^(-1), the stated rational constants, and complete-packet Gram minor exchanges at rationally rotated surrogate nodes. The proof above handles the actual period through 6<2*pi<44/7. The surrogate computation never treats exp(6i) as one and is not a direct exact evaluation at pi.
