# Mixed time and normal-ideal capacity at the true period

2026-11-02. This note adds spatial conditions to the intrinsic time cuts. They apply directly to arbitrary target row sets at common-period centers, including the phase-cycled dilation targets. A finite phase-dilation certificate satisfies all previous time and transport tests but is forced to have zero period determinant by the new spatial capacity test. Its equal-weight source is itself rank deficient; it is not a counterexample in the original separated-weight interpolation regime.

## 1. A family of exact ideal cuts

Let m>=2, let w_1<=w_i for every i, and let

\[
V_H=\operatorname{span}\{Y^hX^\alpha:
w_0h+\textstyle\sum_iw_i\alpha_i\le H\},\qquad N_H=\dim V_H.
\]

Let W be any selected polynomial column subspace of dimension M, contained in V_H. It need not be a monomial space for the abstract cut below. Take any M distinct target labels J={(j,t,a)}, with t>=0 and a in Z_{≥0}^m. Their common-period functionals are

\[
B^Z_{j,t,a}(P)=[z^tu^a]P(e^z,jZ\mathbf1+z\mathbf1+u).
\tag{1}
\]

At Z=2πi these are the actual normalized period rows. All ideal assertions also hold at arbitrary formal or complex Z; such Z need not be an exponential period.

Choose a nonempty subset S of {2,...,m} and define

\[
I_S=(X_i-X_1:i\in S),\qquad
J_{k,r}=(Y-1)^kI_S^r\quad(k,r\ge0).
\]

For a target derivative a set

\[
\eta_S(a)=a_1+\sum_{i\in S}a_i.
\]

Under substitution (1), Y-1 becomes e^z-1, with order one in z. Each X_i-X_1 becomes u_i-u_1 and is independent of z and of the center. Consequently every polynomial in J_{k,r} is killed by the target row unless simultaneously t>=k and η_S(a)>=r. Put

\[
R_{k,r}=\#\{(j,t,a)\in J:t\ge k,\ \eta_S(a)\ge r\},
\qquad d_{k,r}(W)=\dim(W\cap J_{k,r}).
\]

Restrict the matrix map first to W∩J_{k,r}; its image has dimension at most R_{k,r}. The quotient has dimension M-d_{k,r}(W). This proves the exact finite bound

\[
\boxed{\operatorname{rank}B^Z_{J,W}
\le M-d_{k,r}(W)+R_{k,r}.}
\tag{2}
\]

In particular, nonzero square target determinants require d_{k,r}(W)<=R_{k,r} for every k,r,S. The r=0 case is the familiar intrinsic time cut. Positive r adds information about derivatives normal to the diagonal containing every center and the complete time curve.

The test uses a necessary condition for a row to survive; it does not assert that a surviving row is nonzero or independent. It can therefore be weak but cannot overestimate an obstruction.

## 2. The full-space ideal dimension is an exact lattice count

Set y=Y-1, T=X_1 and U_i=X_i-X_1 for i>=2. Give y,T,U_i weights w_0,w_1,w_i respectively. This affine coordinate change preserves the weighted column space in both directions: replacing X_i by U_i+T never raises weight because w_1<=w_i, and the inverse replacement U_i=X_i-X_1 has the same property. Translating Y by one also only lowers weights.

It follows that V_H has the adapted monomial basis y^h T^b U^ν with

\[
w_0h+w_1b+\sum_{i=2}^m w_i\nu_i\le H.
\]

The ideal J_{k,r} is monomial in these coordinates. Therefore

\[
\boxed{D_{k,r}:=\dim(V_H\cap J_{k,r})
=\#\left\{(h,b,\nu):
w_0h+w_1b+\sum_{i=2}^m w_i\nu_i\le H,
h\ge k,\ \sum_{i\in S}\nu_i\ge r\right\}.}
\tag{3}
\]

The elementary subspace intersection inequality gives

\[
d_{k,r}(W)\ge M-N_H+D_{k,r}.
\]

Combining it with (2) produces a column-choice-independent rank bound

\[
\boxed{\operatorname{rank}B^Z_{J,W}
\le \min\{M,\ N_H-D_{k,r}+R_{k,r}\}.}
\tag{4}
\]

Equivalently, the target rows outside the surviving set must fit into the quotient V_H/(V_H∩J_{k,r}). This statement is also valid on all ambient columns at once, with M replaced by N_H in the domain dimension argument and the row count kept separately.

For an actual chosen raw column space, computing d_{k,r}(W) can strengthen (4). One expands the chosen raw columns in the adapted monomial basis and retains the coefficient rows outside the ideal; if that quotient matrix has rank q, then d_{k,r}(W)=M-q exactly.

The coordinate change is used only to compute dimensions and ranks. No Euclidean norm, determinant magnitude or arithmetic denominator claim is transferred through it.

## 3. An exact phase-dilation example with zero target period determinant

Use the full source packet with parameters

\[
m=2,\quad K=2,\quad H=2000,\quad
v_0=1,\quad v_1=v_2=400,
\]

\[
w_0=24,\quad w_1=w_2=40,\quad
\theta=1/10,\quad \lambda=3/2.
\]

Apply the preceding phase-cycled construction: N=s+a+b, retain the source row for N<12, and otherwise map its transverse indices to

\[
\left(\lfloor3a/2\rfloor+(N\bmod3),
\lfloor3b/2\rfloor+(\lfloor N/3\rfloor\bmod3)\right).
\]

The target time is N minus the sum of these two new indices. The exact counts are

| Quantity | Value |
|---|---:|
| Complete source rows | 28,000 |
| Distinct target rows | 28,000 |
| All allowed raw columns | 37,723 |
| Selected raw columns | 28,000 |
| Nontrivial intrinsic time cuts checked | 84 |
| Minimum time-cut slack | 560 |
| Different positive Pascal blocks checked | 45 |
| Largest target total transverse order | 10 |

The selected columns are the 28,000 with smallest total weighted transverse degree, breaking ties by their transverse indices and then h. The two coordinate transport totals are (36,370,36,438), and the selected-column minus target totals are (207,653,212,296), all nonnegative.

The all-time-cut sufficient counting inequality is 15,120<=28,000. The stronger strict parameter condition \(w_0>[2(3-1)+2](2+1)v_0\) is also satisfied, since 24>18. Every bucket Pascal determinant equals its positive Newton product. Thus the example passes both the finite transport and the intrinsic time-support checks.

Now take U=X_2-X_1 and the ideal I^9=(U^9), with k=0. Only 16,548 adapted ambient columns have U-degree less than 9. Exactly 958 target rows have total transverse order at least 9. Formula (4) proves

\[
\boxed{\operatorname{rank}B^{2\pi i}_{J,V_H}
\le16548+958=17506<28000.}
\tag{5}
\]

Thus every 28,000-column target period determinant is zero, with rank deficiency at least 10,494. This is an exact dimension certificate, not a numerical computation of a 28,000-square determinant. The simpler cutoff U-degree>10 already gives rank<=19,405, but (5) also uses the fact that very few rows reach orders 9 and 10.

**Scope limitation:** the source in this example has total transverse order at most 4. It consequently sees at most 9,972 ambient adapted columns and is itself rank deficient at the period. Equal spatial weights are known to be outside the original critical separated-weight interpolation regime. The example proves that the time/transport criteria alone do not ensure period rank; it does not invalidate an argument assuming a nonzero source minor or supply a counterexample in the original separated-weight parameter range.

## 4. The phase targets retain an equal-weight high-dimensional rank obstruction

Suppose all spatial weights equal w and take S={2,...,m}. Then η_S(a)=|a|. If every target has |a|<=AH/w+O(1), with fixed 0<A<1, all adapted columns of normal total degree greater than this bound are invisible. Counting the visible columns gives, at fixed m and weights before H tends to infinity,

\[
\frac{\#\text{visible ambient columns}}{N_H}
\longrightarrow I_A(m-1,3)
:=\frac{\int_0^A x^{m-2}(1-x)^2\,dx}
{\int_0^1 x^{m-2}(1-x)^2\,dx}.
\tag{6}
\]

Here the m-1 normal coordinates share total normalized budget x, while time, tangent and unused budget are the other three simplex coordinates. Direct integration gives

\[
I_A(m-1,3)=A^{m-1}\left[
\frac{m(m+1)}2-(m^2-1)A+\frac{m(m-1)}2A^2\right].
\]

If M_H/N_H tends to a fixed density δ>0, then for every chosen column family,

\[
\boxed{\limsup_{H\to\infty}
\frac{\operatorname{rank}B^Z_{J_H,V_H}}{M_H}
\le \min\{1,\delta^{-1}I_A(m-1,3)\}.}
\tag{7}
\]

The phase-cycled targets obey the needed bound with A=λθ: the additive phase contributes only O(1) at fixed parameters, and the low-order correction decreases the transverse indices. At the critical half-density, the right side is min{1,2I_A(m-1,3)}. For fixed A<1 it tends to zero in the subsequent dimension limit. The same conclusion holds for A_m tending to any fixed value below one, with the height limit taken separately at each m.

This extends the visible-normal-profile obstruction to arbitrary, possibly non-lower and phase-dependent target packets; it does not use an invertible normal-coordinate change on their row labels. It remains an equal-weight obstruction. Arbitrarily separated weights do not preserve total-order row cutoffs on the same normalized scale, so (7) cannot be silently applied to the original separated-weight path.

## 5. Verification and remaining issue

`verify_mixed_normal_cuts.py` uses only the Python standard library. It independently expands raw columns into the adapted coordinates, verifies 144 exact selected-column ideal/rank inequalities on small matrices, and checks the full 28,000-row phase certificate, every intrinsic time cut, all different Pascal blocks, coordinate margins and 4,420 mixed ideal bounds. The best tested ambient bound is (5).

The general ideal inequalities are proved by annihilation and dimension, not extrapolated from finite experiments. The finite large example never asserts the rank equals its upper bound. What remains open here is whether mixed or finer spatial ideal tests force a useful deficit when the source is genuinely full rank in the original strongly separated regime.
