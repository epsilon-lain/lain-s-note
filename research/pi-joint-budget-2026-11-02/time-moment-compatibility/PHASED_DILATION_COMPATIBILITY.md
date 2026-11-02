# Moment-optimal profiles compatible with every time cut and nonzero transport

2026-11-02. A phase-cycled coordinate dilation upgrades the earlier mean/collision witness. In a separated-weight regime it respects the conserved buckets, coordinatewise transport, distinct target rows, every intrinsic time-cut inequality for any selected raw-column family, and a nonzero multivariate Pascal transport minor. Its type multiplicities approach the sharp cone profile. This does **not** prove that a target period minor is nonzero, control that minor's magnitude, or supply an actual rational approximation tuple.

## 1. Fixed-parameter construction

Fix dimension \(m\), centers \(0\le j<K\), positive weights with

\[
v_i=w_i/\theta>v_0,\qquad 0<\theta<1,
\]

and the complete source

\[
L_H=\{(j,s,\beta):v_0s+v\cdot\beta<H\}.
\]

Choose a rational dilation \(\lambda=L/Q>1\), with integers \(L>Q\ge1\), so that

\[
A=\lambda\theta<1.
\]

Write \(w_*\!:=\min_iw_i\), \(S=m(L-1)\), and choose a periodic enumeration

\[
r(N)\in\{0,\ldots,L-1\}^m,
\]

which visits each phase vector exactly once in a period of length \(L^m\). Define

\[
N_0=\left\lceil AH/w_*\right\rceil+S.
\]

For a source row put \(N=s+|\beta|\). If \(N<N_0\), keep the row unchanged. Otherwise set

\[
\boxed{a_i=\lfloor\lambda\beta_i\rfloor+r_i(N),
\qquad t=N-|a|.}
\tag{1}
\]

The resulting target set is denoted \(J_H\).

Every source type satisfies \(|\beta|<\theta H/w_*\). Thus \(|a|\le\lambda|\beta|+S<N_0\) on a changed bucket, and \(t\ge0\). Also \(a_i\ge\beta_i\), and \(t+|a|=s+|\beta|\) exactly. Within each fixed \((j,N)\), the coordinate maps in (1) are strictly increasing and hence injective. Different buckets or centers retain distinct labels, so all target rows are distinct.

For sufficiently large \(H\), all target types satisfy \(w\cdot a\le H\): their degree is at most \(AH+\sum_iw_i(L-1)\), and \(A<1\). Thus the type universe used in the moment relaxation is respected.

## 2. The full Pascal transport minor is nonzero

At a conserved order \(N\), the source type set is

\[
B_N=\{\beta\ge0:|\beta|\le N,
\ \sum_i(v_i-v_0)\beta_i<H-v_0N\}.
\]

Because \(v_i>v_0\), this is a coordinatewise lower set. On a changed bucket let

\[
\xi_i(k)=\lfloor\lambda k\rfloor+r_i(N).
\]

For \(\beta\in B_N\), consider the multivariate Newton polynomial

\[
\mathcal N_\beta(x)=
\prod_i\frac{\prod_{k=0}^{\beta_i-1}(x_i-\xi_i(k))}{\beta_i!}.
\]

Its expansion in the binomial basis \(\binom{x}{\gamma}\), \(\gamma\le\beta\), is unit triangular. Evaluating at the nodes \(\xi(\kappa)\), \(\kappa\in B_N\), gives zero unless \(\kappa\ge\beta\) coordinatewise. With a common linear extension of that order, the determinant is the product of its nonzero diagonal entries. Hence

\[
\boxed{
\det\left[\binom{\xi(\kappa)}\beta\right]_{\beta,\kappa\in B_N}
=\prod_{\beta\in B_N}\prod_i
\frac{\prod_{k=0}^{\beta_i-1}
 (\lfloor\lambda\beta_i\rfloor-\lfloor\lambda k\rfloor)}{\beta_i!}>0.}
\tag{2}
\]

Translations by the phase cancel from the formula. Identity buckets have determinant one. Therefore the complete ordinary-coefficient transport minor is a nonzero integer coefficient times \((-\zeta)^\Gamma\), with the consistent bucket ordering, whenever every relevant \(\zeta_i\ne0\). This conclusion is proved by interpolation; it is not inferred merely from an entrywise matching.

The extra Pascal coefficient \(\kappa_H\) is also inexpensive in separated weights. Since

\[
\beta_i-k\le\lfloor\lambda\beta_i\rfloor-\lfloor\lambda k\rfloor
\le\lambda(\beta_i-k)+1,
\]

each diagonal factor in (2) is at least one, while its logarithm is at most

\[
\beta_i\log\lambda+\lambda^{-1}\sum_{d=1}^{\beta_i}\frac1d.
\]

Consequently, with the parameters fixed before \(H\to\infty\),

\[
0\le\limsup_{H\to\infty}\frac{\log\kappa_H}{M_HH}
\le\frac{\theta\log\lambda}{m+2}\sum_i\frac1{w_i}
\le\frac{\theta\log\lambda}{w_*}.
\tag{3}
\]

The harmonic-sum remainder is \(O(M_Hm\log H)\), hence disappears on this normalization.

## 3. Every intrinsic raw-column time cut is satisfied

Fix any \(M_H\) distinct allowed raw columns \(Y^hX^\alpha\). Let \(k_\alpha\) be their multiplicity at type \(\alpha\), and

\[
d_r=\dim(W\cap(Y-1)^r\mathbb C[Y,X])
=\sum_\alpha(k_\alpha-r)_+.
\]

Since \(k_\alpha\le T_H:=\lfloor H/w_0\rfloor+1\), for \(1\le r\le T_H\),

\[
M_H-d_r=\sum_\alpha\min(r,k_\alpha)
\ge\frac{M_Hr}{T_H}.
\tag{4}
\]

This estimate does not require the selected raw subset to be downward closed.

Let \(Q_H^{\rm src}\) count the source types \(\beta\) with \(v\cdot\beta<H\). For a fixed \(j,\beta\), unchanged rows with target time \(t<r\) have source time in \([0,r-1]\). Changed rows satisfy

\[
t=s-d_0(\beta)-\sigma(N),\quad
d_0(\beta)=\sum_i\lfloor\lambda\beta_i\rfloor-|\beta|,
\quad 0\le\sigma(N)\le S.
\]

If \(0\le t<r\), their source time lies in an integer interval of length \(r+S\). Therefore

\[
\#\{J_H:t<r\}\le KQ_H^{\rm src}(2r+S)
\le KQ_H^{\rm src}(S+2)r.
\tag{5}
\]

Combining (4) and (5) proves all cuts whenever

\[
KQ_H^{\rm src}(S+2)T_H\le M_H.
\tag{6}
\]

For cuts above \(T_H\), \(d_r=0\) and there is nothing to check. The fixed-weight simplex asymptotics give

\[
\frac{KQ_H^{\rm src}T_H}{M_H}
\longrightarrow\frac{(m+1)v_0}{w_0}.
\]

Thus the simple strict condition

\[
\boxed{w_0>(S+2)(m+1)v_0}
\tag{7}
\]

implies, for every sufficiently large height, simultaneously for every selected allowed raw-column family,

\[
\#\{J_H:t\ge r\}\ge d_r\quad\text{for every }r.
\]

Coordinatewise transport already gives nonnegative crossing counts, so this is exactly the entire family of necessary cut inequalities. It remains a time-support statement and does not assert actual period rank.

## 4. Why cycling the phases gives the sharp cone profile

Temporarily apply the changed-bucket rule even at \(N<N_0\), allowing negative target times only in this auxiliary comparison. Its type multiplicities are denoted \(\widetilde n_a\). For a fixed type \(a\) away from the finitely thin coordinate boundary, exactly \(Q\) phases in each coordinate have a preimage under

\[
a_i=\lfloor L\beta_i/Q\rfloor+r_i.
\]

There are therefore exactly \(Q^m\) compatible phase vectors. Each has one preimage \(\beta=a/\lambda+O_{m,L,Q}(1)\), and its conserved orders lie in one residue class modulo \(L^m\). The admissible source-order interval has length

\[
\frac{(H-v\cdot\beta)_+}{v_0}
=\frac H{v_0}\left(1-\frac{w\cdot a}{AH}\right)_++O(1).
\]

All constants here are fixed before the height limit. Counting each residue class yields

\[
\frac{\widetilde n_a}{H}
=\frac K{v_0\lambda^m}
\left(1-\frac{x_a}{A}\right)_++o(1),
\qquad x_a=\frac{w\cdot a}{H},
\tag{8}
\]

uniformly away from that coordinate boundary. The omitted types are \(O(H^{m-1})\), and all occupancies are \(O(H)\), so they do not change the normalized moments or square cost.

Since

\[
\frac{M_H}{Q_HH}\longrightarrow
c=\frac{K\theta^m}{(m+1)v_0},
\]

equation (8) is exactly \(\widetilde n_a/H\sim c g_A(x_a)\), with the moment-optimal density

\[
g_A(x)=\frac{m+1}{A^{m+1}}(A-x)_+.
\]

In particular

\[
\widetilde t_H\longrightarrow\frac{mA}{m+2},
\qquad
\frac{\sum_a\widetilde n_a^2}{M_HH}
\longrightarrow\frac{2K}{v_0(m+2)\lambda^m}.
\tag{9}
\]

A single fixed dilation without cycling phases would leave lattice holes and would not attain this square cost. Phase cycling is needed for the sharp profile.

## 5. Returning to valid nonnegative-time targets

Only source rows with \(N<N_0\) were changed in passing from the auxiliary rule to the valid construction. Their number \(M_H^{\rm low}\) satisfies

\[
M_H^{\rm low}\le K N_0Q_H^{\rm src},\qquad
\limsup_{H\to\infty}\frac{M_H^{\rm low}}{M_H}
\le\delta_*:=\frac{(m+1)v_0 A}{w_*}.
\tag{10}
\]

The valid mean is at most the auxiliary mean, since the replacement is by \(\beta\le a\). Their difference is at most \((A-\theta)M_H^{\rm low}/M_H+o(1)\). For the square cost, both valid and auxiliary occupancies are at most \(K\lceil H/v_0\rceil\): each \((j,N)\) contributes at most one copy of a given type, and \(N<H/v_0\) because \(v_i>v_0\). Also \(\sum_a|n_a-\widetilde n_a|\le2M_H^{\rm low}\). Hence

\[
\limsup_{H\to\infty}
\frac{\left|\sum_an_a^2-\sum_a\widetilde n_a^2\right|}{M_HH}
\le\frac{4K}{v_0}\delta_*
=\frac{4K(m+1)A}{w_*}.
\tag{11}
\]

By choosing the transverse weights sufficiently large at each fixed dimension, both errors can be made arbitrarily small. Thus the actual valid target packets attain the moment-optimal mean and square cost in the separated-weight limit, while retaining all cuts and the nonzero transport minor. This means: choose weights for a prescribed tolerance, then take the height limit. No simultaneous limit with an uncontrolled phase period or weight size is asserted.

The same construction also permits nearly all rows to escape in the subsequent dimension limit. To see the fixed-dimensional calculation, use \(\tau=v_0s/H\), \(y_i=v_i\beta_i/H\), and put

\[
\lambda_i=\lambda-(\lambda-1)v_0/v_i>1.
\]

Away from the low-order correction, a target remains in the source packet precisely when \(\tau+\sum_i\lambda_i y_i<1\), up to a boundary error tending to zero with height. Inside the source simplex \(\tau+\sum_i y_i<1\), this region has relative volume \(\prod_i\lambda_i^{-1}\). Low-order correction changes membership for at most \(M_H^{\rm low}\) rows. Therefore the outside fraction differs from \(1-\prod_i\lambda_i^{-1}\) by at most \(\delta_*+o(1)\). After transverse separation it tends to \(1-\lambda^{-m}\). Negative auxiliary times cause no additional error: they occur only in the corrected low-order buckets.

## 6. Consequence for the critical fixed-exponential regime

Take the preceding normalization

\[
w_0=B^{-m},\qquad K=\lfloor C^m\rfloor,
\qquad v_0=2K\theta^m w_0,
\]

with rational \(\theta\), \(C\theta<1\), and \(\theta<B<1\). Fix a rational \(t\in(B,1)\), and for sufficiently large \(m\) choose

\[
A_m=\frac{m+2}{m}t<1,\qquad \lambda_m=A_m/\theta>1.
\]

For rational fixed \(t,\theta\), the numerator \(L_m\) can be chosen \(O(m)\). Thus \((S+2)(m+1)=O(m^3)\), whereas \(w_0/v_0=1/(2K\theta^m)\) grows exponentially. Condition (7) eventually holds. Choose the transverse weights large enough that the errors in (10)–(11) and the Pascal cost (3) tend to zero along the dimension sequence, and then take sufficiently large heights at each fixed parameter member.

The valid targets have mean tending to \(t\), and their normalized collision cost tends to zero, because its cone value is

\[
\frac{1}{(m+2)w_0A_m^m}
=\frac1{m+2}\left(\frac B{A_m}\right)^m\longrightarrow0.
\]

Their outside fraction tends to one as well, since \(\lambda_m\to t/\theta>1\) and the separation errors can be chosen to tend to zero.

The coordinatewise target means also approach \(t/m\) after rescaling by the weights. The radial cheap-transverse column selection from the earlier witness has coordinate means \(u_m/m\), with \(u_m\to1\). Taking the separation error smaller than each coordinate margin preserves \(\Gamma_i\ge0\) and \(E_i\ge0\) on those same columns.

Consequently, adding all time cuts, exact bucket conservation, coordinatewise matching, and nonvanishing of the full transport minor does not remove the moment relaxation's positive endpoint witness in this regime. The possible remaining improvements include actual-period rank on the same columns, determinant size, further spatial/period structure, and cancellation in the complete aggregate. This is not an actual rational approximation construction or a counterexample to the underlying determinant endpoint.

## 7. Finite verification

The companion verifier constructs a complete two-dimensional integer packet, its phase-cycled matching, and a selected raw-column family. It checks target distinctness, every conserved bucket and coordinate inequality, all intrinsic time cuts, the exact sufficient counting inequality, and all distinct bucket Pascal determinants by rational arithmetic. It does not infer period-minor nonvanishing from those checks.
