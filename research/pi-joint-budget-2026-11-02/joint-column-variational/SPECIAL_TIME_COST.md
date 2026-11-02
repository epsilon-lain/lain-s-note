# Exact joint cost of time and transverse degree

2026-11-02. This note optimizes a positive sum of raw-column degree costs. It gives a sharp joint column-selection law, including a uniform estimate when the time-cost coefficient grows with dimension. It is not a lower bound for a determinant and does not certify the rank or arithmetic admissibility of a minimizing column set.

## 1. The joint simplex model

For raw columns Y^h X^alpha with w_0 h+w dot alpha<=H, put

\[
T=\frac{w_0h}{H},\qquad X=\frac{w\cdot\alpha}{H},\qquad
Z=1-T-X.
\]

At each fixed m and positive weights, the uniform column-counting measure tends as H grows to

\[
(T,X,Z)\sim\operatorname{Dirichlet}(1,m,1).
\]

Equivalently, (T,X) has density m(m+1)x^{m-1} on t>=0, x>=0, t+x<=1. The time coordinate and the budget slack are separate coordinates. Their common mean is 1/(m+2), while E[X]=m/(m+2).

Let a>0 and fix 0<lambda<1. Define

\[
C_a=X+aT,
\qquad
Q_m(a,\lambda)=\min_{\mathbb P(S)=\lambda}
\mathbb E[C_a\mid S].
\tag{1}
\]

The minimum is attained by the lower cost quantile S={C_a<=q}, where F_a(q)=lambda. The distribution is continuous and strictly increasing on its support. Thus this optimizer is unique up to null sets. Exchange of two sets on opposite sides of the cutoff proves optimality directly.

If M_H/N_H tends to lambda, the minimum normalized average over M_H distinct raw columns tends to Q_m(a,lambda). This follows from weighted lattice integration and the zero measure of a cost level set. Selecting the cheapest columns gives the matching upper construction. All m-dependent parameter choices below occur after this fixed-m H limit.

## 2. Exact cost distribution and conditional mean

Write n=m+1. For a>1 the cumulative distribution is

\[
\boxed{
F_a(c)=
\begin{cases}
c^n/a,&0\le c\le1,\\
1-\dfrac{(a-c)^n}{a(a-1)^{n-1}},&1\le c\le a.
\end{cases}}
\tag{2}
\]

It is zero below zero and one above a. For a=1, F_1(c)=c^n on [0,1]. For 0<a<1,

\[
\boxed{
F_a(c)=\frac1a\left[c^n-
\frac{(c-a)_+^n}{(1-a)^{n-1}}\right],
\qquad0\le c\le1.}
\tag{3}
\]

To verify these formulas, integrate m(m+1)x^{m-1} over

\[
0\le x\le\min(c,1),\qquad
0\le t\le\min\{1-x,(c-x)/a\}.
\]

The two affine boundaries cross at x=(a-c)/(a-1). This gives (2)--(3) without treating T and X as independent.

For every a, the conditional mean is exactly

\[
\boxed{Q_m(a,\lambda)
=q-\frac1\lambda\int_0^qF_a(c)\,dc,\qquad F_a(q)=\lambda.}
\tag{4}
\]

For a>=1 this becomes an explicit two-branch formula. If a lambda<=1,

\[
q=(a\lambda)^{1/n},\qquad
\boxed{Q_m(a,\lambda)=\frac n{n+1}(a\lambda)^{1/n}.}
\tag{5}
\]

If a>1 and a lambda>=1, put L=1-lambda and

\[
\delta=[a(a-1)^{n-1}L]^{1/n},\qquad q=a-\delta.
\]

Then

\[
\boxed{Q_m(a,\lambda)
=\frac{(m+a)/(m+2)-L\{a-n\delta/(n+1)\}}\lambda.}
\tag{6}
\]

The formulas agree at a lambda=1. In the upper branch, the discarded tail is an affine simplex at the time vertex. Its mean cost is a-n delta/(n+1), which also proves (6).

For 0<a<1, the unique q solving (3) gives the algebraic formula

\[
\boxed{Q_m(a,\lambda)
=q-\frac{q^{n+1}-(q-a)_+^{n+1}/(1-a)^{n-1}}
{a(n+1)\lambda}.}
\tag{7}
\]

No numerical root choice is required to define these values.

The joint coordinate means are also explicit where (5) applies:

\[
\mathbb E[T\mid C_a\le q]=\frac{q}{a(m+2)},\qquad
\mathbb E[X\mid C_a\le q]=\frac{mq}{m+2}.
\tag{8}
\]

For the upper branch (6), put b=a-1 and k=delta/b. The selected means are

\[
\mathbb E[X\mid C_a\le q]
=\frac{m(1-Lk)}{\lambda(m+2)},
\qquad
\mathbb E[Z\mid C_a\le q]
=\frac{1-L\delta/a}{\lambda(m+2)},
\tag{9}
\]

and the selected mean of T is one minus their sum. These formulas retain the exact relation T+X+Z=1.

## 3. A bound uniform in the time coefficient

The smallest possible average time coordinate among a lambda fraction is

\[
\boxed{\tau_m(\lambda)
=\frac{1-(n+1)L+nL^{1+1/n}}{(n+1)\lambda},
\qquad n=m+1,\quad L=1-\lambda.}
\tag{10}
\]

Indeed, T has CDF 1-(1-t)^n. Its lower lambda-quantile is u=1-L^{1/n}; direct integration up to u gives (10). Selecting the smallest T values gives

\[
\mathbb E[X\mid T\le u]
=\frac m{m+1}\{1-\tau_m(\lambda)\},
\]

because conditional on T the other m+1 simplex coordinates remain symmetric.

For every event S of mass lambda,

\[
\mathbb E[X\mid S]
=1-\mathbb E[T+Z\mid S]
\ge1-\frac{2}{\lambda(m+2)}.
\]

Combining this with E[T|S]>=tau_m yields the following **uniform-in-a** estimate for all a>0:

\[
\boxed{
1+a\tau_m-\frac{2}{\lambda(m+2)}
\le Q_m(a,\lambda)
\le\frac m{m+1}(1-\tau_m)+a\tau_m
\le1+a\tau_m.}
\tag{11}
\]

The upper bound is realized by an actual measurable column selection, the cheapest-time set. It is not generally the joint optimizer. For a>=1 a slightly sharper lower bound is

\[
\boxed{Q_m(a,\lambda)
\ge1+(a-1)\tau_m-\frac1{\lambda(m+2)},}
\tag{12}
\]

using C_a=1+(a-1)T-Z. The error constants in (11) do not depend on a, even when a grows exponentially with m.

Set

\[
d_\lambda=1+\frac{1-\lambda}{\lambda}\ln(1-\lambda)>0.
\]

Taylor expansion in (10) gives tau_m=d_lambda/m+O_lambda(m^{-2}). Hence for arbitrary positive a_m,

\[
\boxed{Q_m(a_m,\lambda)
=1+\frac{d_\lambda a_m}{m}
+O_\lambda\left(\frac{a_m}{m^2}+\frac1m\right).}
\tag{13}
\]

In particular:

| Time-cost scale | Joint optimal positive cost |
|---|---|
| a_m=o(m) | Q_m tends to 1 |
| a_m/m tends to kappa in (0,infinity) | Q_m tends to 1+kappa d_lambda |
| a_m/m tends to infinity | m Q_m/a_m tends to d_lambda |
| a_m=A rho^m, rho>1 | Q_m is asymptotic to A d_lambda rho^m/m |

At lambda=1/2, d_lambda=1-ln 2. Thus a_m asymptotic to kappa m forces a limiting joint cost 1+kappa(1-ln 2). If a=R/w_0 is asymptotic to a constant times (CB)^m with CB>1, the minimum of this positive sum itself grows exponentially divided by m. This identifies a cost in the stated norm estimate; it is not a lower bound for the analytic determinant being estimated.

## 4. Fixed a and the transition near a=1

For finer first-order information, the joint scaling is

\[
(mT,mZ)\ \Longrightarrow\ (E_1,E_2),
\]

where E_1,E_2 are independent exponential variables of mean one. This follows directly from the scaled density

\[
\frac{m+1}{m}\left(1-\frac{u+v}{m}\right)^{m-1}
\mathbf1_{u,v\ge0,\ u+v\le m}.
\]

The second moments of each scaled coordinate are bounded by two. This supplies uniform integrability and justifies convergence of conditional first moments below continuous quantiles.

For fixed a>1, put b=a-1. Then m(C_a-1) tends to V=bE_1-E_2. If a lambda<=1, its lower quantile is v=ln(a lambda)<=0 and its selected mean is v-1. If a lambda>=1,

\[
v=b\ln\frac{b}{a(1-\lambda)},\qquad
\mathbb E[V\mid V\le v]
=\frac{a-2-(1-\lambda)(v+b)}\lambda.
\tag{14}
\]

In either branch, Q_m=1+E[V|V<=v]/m+o(m^{-1}). At half density this reads

\[
m(Q_m(a,1/2)-1)\longrightarrow
\begin{cases}
\ln(a/2)-1,&1\le a\le2,\\
a-3-(a-1)\ln\!\left(2(a-1)/a\right),&a\ge2.
\end{cases}
\tag{15}
\]

For 0<a<1, write b=1-a. The nonnegative variable Y=E_2+bE_1 has survival function

\[
S_a(y)=\frac{e^{-y}-b e^{-y/b}}a.
\]

Let y solve S_a(y)=lambda. Then

\[
\boxed{m(1-Q_m(a,\lambda))\longrightarrow
y+\frac{e^{-y}-b^2e^{-y/b}}{a\lambda}.}
\tag{16}
\]

The limiting case a=0 gives the earlier cheapest-transverse-degree Gamma(2) law. For every sequence a_m tending to one from either side, the same density and uniform-integrability argument gives

\[
\boxed{m\{Q_m(a_m,\lambda)-1\}\longrightarrow\ln\lambda-1.}
\tag{17}
\]

No rate condition on a_m-1 is needed for (17).

## 5. What can be made cheap together

At a fixed dimension and a fixed fraction 0<lambda<1, the separate minimum of E[T|S] and the separate minimum of E[X|S] cannot be attained by the same set. Their unique minimizing sets are respectively a time sublevel set and a transverse-degree sublevel set; these differ on a region of positive simplex measure. For example, the latter includes points near the time vertex that the former excludes.

Consequently Q_m(a,lambda) is strictly larger than a times the separate time minimum plus the separate transverse minimum for every finite a>0. However, choosing the cheapest-time columns already leaves mean transverse degree 1-O_lambda(1/m), as does every positive-density column choice. The simultaneous penalty relative to these separate minima is therefore only O_lambda(1/m), uniformly in a. The unavoidable time average tau_m is what becomes expensive when multiplied by a growing coefficient.

No source-degree inequalities, time-grade matching, or nonzero minor condition have been imposed in this optimization. They can only raise the minimum positive cost. A cheapest-column construction in this relaxed problem is not evidence that the same columns support the needed determinant. Likewise this note does not rule out cancellations or a sharper analytic norm that avoids the displayed positive sum.

## Verification

`verify_special_time_cost.py` integrates the joint polynomial density with rational arithmetic independently of the closed CDF formulas, checks the conditional coordinate moments and mean costs, certifies the uniform bounds using rational quantile enclosures, and records exact-enclosure scaling examples. Infinite limits follow from the arguments above, not finite extrapolation.
