# Diagonal quotient capacity for the unshifted test selection

Date: 2026-11-01. This is an unconditional necessary rank condition for the specific selection \(a=\beta,\ell=s\). It does not apply to arbitrary shifted types, and it does not obstruct the independent shell construction.

Let the polynomial domain and original row packet be

\[
P_W(H)=\operatorname{span}\{Y^hX^\alpha:w_0h+w\alpha\le H\},
\qquad
\mathcal R_H=\{(j,s,\beta):0\le j<K,\ v_0s+w\beta/\theta<H\}.
\]

The identity-selected test row at a common diagonal center \(jW\) is

\[
T_{j,s,\beta}(P)
=[t^su^\beta]P(1+t,jW+\log(1+t)+u_1,\ldots,
                    jW+\log(1+t)+u_m).
\]

Put \(w_{\min}=\min_iw_i\), \(w_{\max}=\max_iw_i\), and define the finite integers

\[
L_H=\#\{(s,k)\in\mathbb Z_{\ge0}^2:
             v_0s+(w_{\max}/\theta)k<H\},
\]

\[
D_H=\#\{(h,d)\in\mathbb Z_{\ge0}^2:
             w_0h+w_{\min}d\le H\}.
\]

**Proposition.** If \(KL_H>D_H\), the full identity-selected test matrix is rank deficient for every \(W\), including the actual period \(W=2\pi i\). If the anchor block at \(j=0\) is independent, contraction of that block leaves at most \(D_H-L_H\) dimensions for the \((K-1)L_H\) diagonal row combinations defined below.

**Proof.** For every pair counted by \(L_H\), all multi-indices \(\beta\) of total degree \(k\) occur in the original row packet, because \(w\beta\le w_{\max}k\). Form

\[
U_{j,s,k}=\sum_{|\beta|=k}T_{j,s,\beta}.
\]

The linear combinations \(U_{j,s,k}\) have disjoint, nonempty supports in the original row index set. Hence their coefficient vectors, before applying the evaluation map, are linearly independent.

Let

\[
Q(Y,Z)=P(Y,Z,\ldots,Z).
\]

Setting every transverse variable equal to \(u\) gives the exact identity

\[
U_{j,s,k}(P)=[t^su^k]Q(1+t,jW+\log(1+t)+u).               \tag{1}
\]

The image of the diagonal restriction \(P\mapsto Q\) is exactly

\[
\operatorname{span}\{Y^hZ^d:w_0h+w_{\min}d\le H\}.
\]

One inclusion follows from \(w\alpha\ge w_{\min}|\alpha|\); the reverse inclusion follows by putting all \(d\) powers on a variable of minimum weight. The dimension is therefore exactly \(D_H\). All \(KL_H\) functionals in (1) factor through that space, so they are dependent if \(KL_H>D_H\). Their disjoint source supports turn this into a nontrivial dependence of the original test rows.

If the anchor rows are independent, the \(L_H\) anchor combinations \(U_{0,s,k}\) are independent. The span of all diagonal combinations has dimension at most \(D_H\), and its intersection with the anchor span has dimension at least \(L_H\). Its image after contracting the full anchor span therefore has dimension at most \(D_H-L_H\). This proves the quotient statement. ∎

## Asymptotic form and what it does not prove

For fixed positive weights,

\[
L_H=\frac{\theta H^2}{2v_0w_{\max}}+O(H),
\qquad
D_H=\frac{H^2}{2w_0w_{\min}}+O(H).
\]

Thus

\[
K\theta w_0w_{\min}>v_0w_{\max}                            \tag{2}
\]

forces identity-selection rank deficiency for all sufficiently large \(H\). A necessary asymptotic condition for full rank of this particular selection is the reverse weak inequality.

On the core parameter path \(v_0=2K\theta^mw_0\), condition (2) becomes

\[
\frac{w_{\max}}{w_{\min}}<\frac1{2\theta^{m-1}}.
\]

This shows why a growing separation of transverse weights is necessary even for a simple identity-selection certificate. The manuscript can choose much more strongly separated weights, and a shifted sparse selection no longer has the shell-complete sums in (1). Consequently this proposition neither contradicts the separated-weight interpolation theorem nor rules out the shell construction.

Compared with the earlier \(s=0\) diagonal-slice test, this bound uses all time orders. Its asymptotic threshold is stronger by the factor \(w_0/v_0\).

## A finite example beyond individual-type capacities

Take

\[
m=2,\quad K=2,\quad\theta=\tfrac14,\quad
w_0=1,\quad v_0=\tfrac15,\quad(w_1,w_2)=(1,2),\quad H=100.
\]

Then

\[
L_H=3380,\quad D_H=5151,\quad KL_H=6760>D_H.
\]

The identity test matrix is therefore rank deficient. Yet the \(s=0\) diagonal-slice count is only \(26\), below its capacity \(101\). The fixed-type \(a=0\) count is \(1000\), below its capacity \(5151\), and the verifier checks that **every** occurring individual type satisfies its own capacity bound \(n_a\le D_a(H)\).

The numerical volume conditions \(K\theta^m<1\) and \(Kw_0\theta^m/v_0<1\) both hold. The weights here are intentionally not asserted to meet Theorem 2.1's successive separation thresholds. The example illustrates a genuine collective quotient obstruction that individual-type capacities and the single-time diagonal test both miss.

The companion verifier uses exact integers and rational numbers. It checks these finite counts, every individual-type capacity in this example, and the combinatorial row-sum identity behind (1). It does not rely on a floating-point matrix rank calculation.
