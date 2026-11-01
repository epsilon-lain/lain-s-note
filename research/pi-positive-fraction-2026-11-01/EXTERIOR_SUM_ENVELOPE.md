# Reducing the complete exterior sum to a period-minor envelope

2026-11-01. This observation makes the remaining quantitative target explicit. It is an upper bound by the triangle inequality, not a claim that cancellation is absent or that a largest term dominates from below.

Fix m, K, positive rational weights, 0<theta<=1, and errors with 0<|zeta_i|<1 before H tends to infinity through the admissible weight lattice. Use the source packet in `POSITIVE_FRACTION_LAW.md`, with v_i=w_i/theta>v_0, c_i=v_i-v_0 and M_H source rows. Select any M_H raw columns obeying the polynomial budget. No interpolation theorem is needed to state the following bound; zero minors may occur.

All potentially relevant target rows belong to the finite set

\[
\Omega_H=\{(j,t,a):0\le j<K,\ t\ge0,\
t+|a|<H/v_0,\ w\cdot a\le H\}.
\]

The conserved total order of a source row is less than H/v_0. Targets violating the type cap have zero entries on the selected raw columns and can be omitted. Put D_H=|Omega_H|. A sufficient count is

\[
D_H\le K\left\lceil\frac H{v_0}\right\rceil
\prod_i\left(\left\lfloor\frac H{w_i}\right\rfloor+1\right)
=O(H^{m+1}).
\]

Since M_H is asymptotic to a positive constant times H^(m+1), the number of target subsets contributes no leading MH cost:

\[
\frac{\log\binom{D_H}{M_H}}{M_HH}
\le\frac{D_H\log2}{M_HH}\longrightarrow0.
\tag{1}
\]

If D_H<M_H there are no full target subsets and the determinant is zero. In this setting Omega_H in fact contains the source packet, so D_H>=M_H.

Let T denote the rectangular transport matrix and B the period matrix on Omega_H and the chosen columns. Cauchy--Binet gives

\[
\det(TB)=\sum_{|J|=M_H}\det T_{L,J}\det B_J.
\]

For a nonzero transport minor put rho_J=|J\setminus L|/M_H, and let

\[
\mathcal A=\prod_i\frac{v_i}{c_i},\qquad
\tau=\min_i\frac{-\log|\zeta_i|}{c_i},\qquad
a_0=\frac{\tau}{(m+1)\mathcal A},\qquad
w_* = \min_i w_i.
\]

The proof of the positive-fraction law is uniform in every bucket exit count: its error term is O((h_N+1)^m), independently of the selected target sets. Summing those errors and using the same lattice counts yields a sequence epsilon_H tending to zero, independent of J, such that

\[
\frac{\log|\det T_{L,J}|}{M_HH}
\le\frac{\log2}{w_*}-a_0\rho_J^2+\epsilon_H.
\tag{2}
\]

The errors and weights are fixed here; no uniformity when they vary with H is asserted. Combining (1), (2) and the triangle inequality proves

\[
\boxed{
\frac{\log|\det(TB)|}{M_HH}
\le\frac{\log2}{w_*}
+\max_{J:\det T_{L,J}\det B_J\ne0}
\left\{\frac{\log|\det B_J|}{M_HH}-a_0\rho_J^2\right\}
+o(1).}
\tag{3}
\]

Use log 0=-infinity. If every summand vanishes, the left side is also -infinity. This bound treats the actual period values of B_J, not only their formal-polynomial degrees or finite-field residues.

One can separate J=L before taking a maximum: its transport determinant is exactly one. The original period determinant then keeps its own estimate without the coarse log(2)/w_* charge. The maximum in the remaining bound ranges over J other than L. This refinement still leaves selections with a vanishing but nonzero exit fraction, so a result only for rho bounded away from zero is insufficient.

Thus the complete absolute-value approach has a precise missing input: a uniform upper envelope for the compatible period minors as a function of their escape fraction, including fractions approaching zero. The number of exterior summands itself costs o(MH) at each fixed parameter choice. Rank certificates show that some summands survive; they do not bound this envelope. Conversely, a stronger proof may exploit cancellations and improve on (3), but cannot infer such cancellation from individual nonzero or zero terms.
