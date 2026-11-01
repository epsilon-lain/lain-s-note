# Fractional-dilation packing of nonzero transport minors

2026-11-01. This note proves a local packing theorem for the exact transport expansion. It does not prove a nonzero complete period determinant or close the endpoint contradiction. In particular, the low-total-order rows cannot be discarded during a rank proof.

## 1. A lower-grid evaluation packing theorem

A finite set B in the nonnegative integer lattice is a lower set when beta in B and 0<=alpha<=beta coordinatewise imply alpha in B. Let B and A be finite nonempty lower sets in Z_{≥0}^m. Let t>=1 be real and K>=1 integer, with t^m>=K. Assume

\[
\lfloor tx\rfloor\in A
\quad\text{for every }x\in\bigcup_{\beta\in B}(\beta+[0,1)^m).
\tag{1}
\]

Then A contains K pairwise disjoint subsets A_1,...,A_K, each of cardinality |B|, such that

\[
\boxed{\det\left[\binom a\beta\right]_{\beta\in B,\ a\in A_j}\ne0
\qquad(1\le j\le K).}
\tag{2}
\]

All determinants are nonzero integers. This applies to fractional t; no integer dilation is required.

### Proof: volume, polynomial support, then base packing

For S subset B, let its upward shadow in A be
\(\operatorname{Sh}_A(S)=\{a\in A:a\ge\beta\text{ for some }\beta\in S\}\).
Since t>=1 and x in beta+[0,1)^m, floor(tx)>=beta. Thus the image of the disjoint union of the cubes indexed by S under x->tx lies in the union of unit cubes indexed by Sh_A(S). Comparing volume gives

\[
|\operatorname{Sh}_A(S)|\ge t^m|S|\ge K|S|.       \tag{3}
\]

Let P_B be the polynomial space spanned by the monomials x^beta, beta in B, over Q. It is also spanned by the binomial polynomials binom(x,beta), because the change of basis is triangular with nonzero diagonal. Evaluation on B is invertible, so evaluation on A has rank |B|. The inclusion B subset A follows from (1), by taking x=beta and using floor(t beta)>=beta together with the lower-set property of A.

For a subspace V subset P_B of dimension r, choose a row-echelon polynomial basis with r distinct leading monomials S subset B, using any monomial order. We claim that the union of the nonzero evaluation supports of V in A has at least |Sh_A(S)| points.

Here is the footprint argument in this special setting. The vanishing ideal I(A) has exactly A as its standard monomials. For every beta outside A, the polynomial binom(x,beta) vanishes on A: no point a in the lower set A can satisfy a>=beta. Its leading monomial is x^beta. Hence the standard monomials of I(A) are contained in A; their number is |A|, because evaluation gives Q[x]/I(A) isomorphic to Q^A, so equality follows. Add the polynomials in V to I(A). The new leading ideal contains the complement of A and the upward shadow of S. Its quotient dimension is therefore at most |A|-|Sh_A(S)|. In the product algebra Q^A this quotient dimension is exactly the number of points where every polynomial in V vanishes. This proves the claimed support bound, and (3) yields

\[
|\operatorname{supp}_A V|\ge K\dim V.             \tag{4}
\]

Now consider the vector matroid whose elements are the evaluation columns indexed by A, with rows indexed by B. For any E subset A, the space V of polynomials vanishing on E has dimension |B|-rank(E). Its support is contained in A\E. Consequently

\[
|A\setminus E|\ge K(|B|-\operatorname{rank}E).     \tag{5}
\]

The matroid base-packing theorem gives K disjoint bases, proving (2). The exact standard theorem used is Theorem 4 in Michel Goemans's MIT 18.438 Lecture 13 notes, October 29, 2009:
https://math.mit.edu/~goemans/18438F09/lec13.pdf . Its rank condition is precisely (5). The volume and footprint steps above are proved here rather than assumed as an extension of that theorem.

The assertion is for the full lower target set A. It does not assert that one prescribed basis B can be held fixed or that the other bases lie entirely in A\B. A direct shell version of the footprint argument is false: with B={0,1}, A={0,1,2,3}, the polynomial x-2 has only one nonzero value on {2,3}, whereas the upward shadow of its leading monomial there has two points.

## 2. Application to a conserved total-order bucket

Use source costs v_0>0, v_i=w_i/theta>v_0, and c_i=v_i-v_0>0. Transport preserves

\[
N=s+|\beta|.
\]

For a fixed integer N with L=H-v_0N>0, the source types are the lower set

\[
B_N=\{\beta\ge0:c\beta<L,\ |\beta|\le N\}.
\]

Write c_* = min_i c_i and C_sum=sum_i c_i. Define

\[
A_N=\{a\ge0:ca<t(L+C_{\rm sum}),\ |a|\le N\}.
\tag{6}
\]

If

\[
\boxed{N\ge\frac{t(L+C_{\rm sum})}{c_*},}         \tag{7}
\]

then A_N satisfies the dilation containment (1). Indeed, for x in a source cube, cx<L+C_sum, hence c floor(tx)<t(L+C_sum), and
\(|\lfloor tx\rfloor|\le t cx/c_*<t(L+C_sum)/c_*\le N\).

Thus K disjoint nonzero Pascal minors exist inside A_N. Associate one target set A_{N,j} to each center j and replace its bucket rows by

\[
(s',a)=(N-|a|,a),\qquad a\in A_{N,j}.
\]

The exact transport coefficient of that bucket minor is

\[
(-\zeta)^{\sum_{a\in A_{N,j}}a-\sum_{\beta\in B_N}\beta}
\det\left[\binom a\beta\right]_{\beta\in B_N,a\in A_{N,j}}.
\tag{8}
\]

It is nonzero whenever every zeta_i is nonzero. The exponent vector is nonnegative: a nonzero Pascal determinant contains a permutation with every matched a>=beta. Its matched rows satisfy s'=N-|a|<=N-|beta|=s. No individual boundary coefficient need be asserted nonzero; (8) is the coefficient after the exact bucketwise determinant regrouping.

At the formal period zero, disjoint target sets across centers give distinct (s',a) in this one bucket. If these rows are contained in a complete auxiliary one-center interpolation packet, the manuscript's conditional interpolation theorem certifies their independence. This is a local bucket conclusion and does not address other buckets at the same time.

## 3. Cost bounds and the gluing obstruction

The target support obeys

\[
wa=\theta(ca+v_0|a|)
\le\theta(1+v_0/c_*)t(L+C_{\rm sum}).             \tag{9}
\]

Hence it retains the fractional radial factor t, with an error that tends to zero as v_0/c_* tends to zero, rather than paying an integer factor two.

If this construction is used for every bucket satisfying (7), then for each type a there is at most one selected row across all K centers in that bucket. Condition (7) is equivalent to

\[
N\ge N_0(H):=\frac{t(H+C_{\rm sum})}{c_*+tv_0}.
\]

One can choose identity transport in the remaining buckets; its coefficient is one. The resulting global **transport** minor is nonzero, because it is a product of nonzero bucket minors. It uses distinct rows within each center and respects every total-order budget.

For this formal selection, type occupancy satisfies the crude bound

\[
n_a\le \lceil H/v_0\rceil+K\lceil N_0(H)\rceil,
\]

so a collision score c_collision sum_a n_a^2/(MH) is at most

\[
c_{\rm collision}\left(\frac1{v_0}+
\frac{Kt}{c_*+tv_0}\right)+o(1).                 \tag{10}
\]

If the quality ratio is a common nu, the corresponding weighted displacement cost is at most

\[
\limsup_{H\to\infty}
\frac{\nu\sum w(a-\beta)}{MH}
\le\frac{\nu\theta}{m+2}
\left[t(1+v_0/c_*)(m+1)-m\right].                \tag{11}
\]

To check (11), use (9), subtract w beta, and note that N=s+|beta|. The resulting expression is bounded above by theta t(1+v_0/c_*)(H-v_0s)-w beta plus an O(1) cube correction per row. Summing uses the original simplex means (H-v_0s)/H ->(m+1)/(m+2) and w beta/H ->m theta/(m+2). On identity buckets the displacement is zero and the same majorant is nonnegative, so extending the sum to all source rows is valid. The O(1) per-row corrections disappear after dividing by MH.

These estimates show that the transport constraints and Pascal nonvanishing alone do not force a large scalar-plus-occupancy cost. For the earlier family theta=1/2, t->13/10, nu->2, v_0->infinity and c_* chosen much larger than Kv_0, the displayed upper costs tend to 0.3, below the old uniform arithmetic threshold 0.5.

**The missing step is substantial:** a nonzero transport minor need not have a nonzero accompanying period minor on the weighted polynomial column space. In the low-N buckets, B_N can fill the entire homogeneous simplex. K copies then cannot be made disjoint at formal period zero. Retaining those rows preserves the nonzero transport coefficient but does not supply a common zero-period rank certificate. They may form only a small fraction of all rows; that does not permit dropping them from a determinant nonvanishing proof. No global nonzero low-cost determinant term is claimed here.

## 4. Verification

The companion script uses exact rational arithmetic and integer determinants. It checks all shadow subsets and all matroid rank inequalities for a small two-dimensional example, constructs explicit disjoint Pascal bases for anisotropic conserved-order buckets, certifies the nonzero exact transport coefficient, and checks the nonnegative permutation matching and column-type support. These are finite checks of the theorem and its interface, not numerical evidence for pi's bad approximability.

## 5. Exact formulation of the remaining rank-selection problem

Fix a finite ambient target row set that contains the identity rows and all the packed rows. On it define two matroids:

1. The direct sum, over each center and conserved total order, of the Pascal evaluation matroids. A basis of this matroid is exactly a size-M target row selection with nonzero full transport coefficient.
2. The linear matroid of the actual-period row functionals on the full weighted polynomial column space P_W(H), truncated to rank M if necessary. A size-M independent set is exactly a target selection with full row rank on that space.

Thus a nonzero transport coefficient together with a full-row-rank period test requires a common base. The identity target selection is such a base when the period source interpolation theorem applies. The packing theorem constructs a low-cost base of the first matroid; it does not prove independence in the second. A suitable low-cost common-base result would close this particular rank gap.

This formulation is on the full column space. A subsequent common bounded-column-mixing step can make the rational source and the selected period test simultaneously nonsingular. It must be charged using the appropriate mixed-column arithmetic and analytic estimates. It does not establish that both are nonzero on an arbitrary preselected raw minor, nor preserve the raw-minor degree prefactor without retaining the full column scaling.

A small proportion of unresolved low-order rows does not automatically give a small-cost repair of a common base: exchange paths can involve many other rows. No such repair bound is claimed.
