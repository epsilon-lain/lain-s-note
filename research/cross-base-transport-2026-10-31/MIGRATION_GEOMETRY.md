# Migration geometry: digit sensitivity, invariants, and arithmetic fibers

Date: 2026-10-31

We study the fixed π-generated digit word d_n in {0,1,2,3}:
\[
1=\sum_{n\ge1}d_n\pi^{-n},
\qquad
F(t)=\sum_{n\ge1}d_nt^{-n},\quad t>1.
\]
The digits are fixed while the positional base t varies. No decimal digits or continued fractions are used.

## 1. Sensitivity is a tilted position distribution

For t>1 define
\[
P_t(n)=\frac{d_nt^{-n}}{F(t)}.
\]
Terms with d_n=0 have probability zero. Since F(t)>0, this is a probability distribution on the support S={n:d_n>0}.

Let x=log t and Z(x)=F(e^x). Then
\[
P_x(n)=\frac{d_ne^{-nx}}{Z(x)}.
\]
All derivatives converge locally uniformly for x>0, since d_n≤3 and polynomial factors n^k are dominated by e^{-nx}. Define
\[
\mu(x)=\mathbb E_x[N]=\sum nP_x(n),\qquad
V(x)=\operatorname{Var}_x(N).
\]
Differentiation gives
\[
\boxed{\frac{d}{dx}\log Z(x)=-\mu(x),\qquad
\mu'(x)=-V(x)<0.}
\tag{1}
\]
The variance is strictly positive because d_1=3 and d_3=1, so the support contains at least two positions for every x.

Equivalently,
\[
-\frac{d\log F(t)}{d\log t}=\mathbb E_t[N],
\qquad
\frac{d}{d\log t}\mathbb E_t[N]=-\operatorname{Var}_t(N).
\tag{2}
\]
Thus relative sensitivity is the mean active digit position, and that mean strictly decreases as the base grows.

For a single position n, its contribution is a_n(t)=d_nt^{-n}; its logarithmic sensitivity is exactly
\[
\frac{d\log a_n}{d\log t}=-n.
\]
The relative contribution of position n to the total logarithmic sensitivity is nP_t(n). This gives a precise answer to “which positions matter most”: at base t, rank positions by nP_t(n), not merely by d_n.

## 2. Pairwise migration law and a geometric invariant

For two active positions m,n,
\[
\boxed{
\frac{P_t(n)}{P_t(m)}
=\frac{d_n}{d_m}t^{m-n}.
}
\tag{3}
\]
Hence
\[
\frac{d}{d\log t}\log\frac{P_t(n)}{P_t(m)}=m-n.
\tag{4}
\]
The relative odds of any two positions move linearly in log-base time, with a slope depending only on their positional separation. The normalization F(t) cancels.

This gives a family of transport invariants:
- the support S={n:d_n>0};
- the positional gaps n−m;
- the pairwise log-odds slopes m−n;
- the full coefficient word d, recoverable as the Taylor coefficients of f(z)=F(1/z)=Σd_nz^n at z=0.

The last point is exact: knowing the entire transport curve F(t) for all sufficiently large t determines every d_n. A single value F(t_0), in contrast, forgets almost all of this information.

## 3. Cumulants and the moment vector field

Let κ_k(x) be the k-th cumulant of the tilted position distribution P_x. The log-partition function is
\[
\log Z(x)=\log\sum d_ne^{-nx}.
\]
Standard differentiation of this convergent exponential family yields
\[
\boxed{
\frac{d^k}{dx^k}\log Z(x)=(-1)^k\kappa_k(x).
}
\tag{5}
\]
In particular, the second derivative is the variance. The third derivative is minus the third cumulant, and so on.

If M_k(x)=Σ n^k d_ne^{-nx}, then
\[
M_k'(x)=-M_{k+1}(x).
\tag{6}
\]
The normalized moment flow is therefore an infinite-dimensional vector field. The measure support and original weights d_n are fixed; only the exponential tilt changes.

The Fisher information of this one-parameter exponential family is
\[
I(x)=V(x).
\]
Thus the intrinsic statistical line element is
\[
ds_{\rm Fisher}^2=V(x)\,dx^2.
\]
This is invariant under smooth reparametrization of the base coordinate in the usual metric sense. It is not constant along the curve; it measures how distinguishable nearby base parameters are through the digit-position distribution.

## 4. Global arithmetic fibers

The earlier transport note proves F is a strictly decreasing bijection (1,∞)→(0,∞). Therefore every positive target y has a unique base
\[
\tau(y)=F^{-1}(y).
\]
For any subset E⊂(0,∞), its parameter fiber is exactly F^{-1}(E).

Consequences:
- positive rational targets have a countable dense set of unique base parameters;
- positive algebraic targets have a countable dense set of unique base parameters;
- transcendental targets correspond to an uncountable, full-Lebesgue-measure, comeager set of base parameters.

The algebraic-valued and transcendental-valued parameter sets are both dense. They are not intervals separated by analytic singularities.

For integer bases b≥4, the fixed digits are valid base-b digits. Since d is not eventually periodic, F(b) is irrational. This includes b=10. The proof does not establish transcendence of F(10), nor classify F(t) at rational noninteger t.

## 5. Arithmetic exceptional points are not analytic singularities

F is real analytic on (1,∞), with
\[
F^{(k)}(t)=(-1)^k\sum_{n\ge1}
(n)_k d_n t^{-n-k},
\]
where (n)_k=n(n+1)\cdots(n+k-1) is the rising factorial. In particular, no point t>1 is a singularity of the real-analytic transport function.

Thus “arithmetic singularity” should mean a parameter where F(t) lands in a specified arithmetic class, not a singularity of F as an analytic function. These arithmetic fibers can be dense even though F is smooth everywhere on its domain.

## 6. Inverse migration vector field

Since F'(t)<0, the inverse τ(y) is differentiable:
\[
\boxed{
\tau'(y)=
-\left(\sum_{n\ge1}n d_n\tau(y)^{-n-1}\right)^{-1}.
}
\]
At y=1, τ(1)=π. This is a reverse flow from a target value to its unique base parameter.

## 7. Open explicit-value questions

The measure/category theorem is nonconstructive with respect to a specified parameter. It does not decide whether F(10) is transcendental. The following remain open in this project:
1. Prove or disprove that F(10) is transcendental.
2. Classify F(p/q) for rational noninteger p/q>1.
3. Obtain effective bounds on the block complexity of d, then test whether known complexity criteria for algebraic expansions apply.
4. Study the Fisher length and cumulant asymptotics as t↓1 and t→∞.
5. Compare the fixed-word transport F(t) with the regenerated greedy expansion of 1 in base t; the latter changes the digit word itself.

No proof or disproof of P2848 is claimed.
