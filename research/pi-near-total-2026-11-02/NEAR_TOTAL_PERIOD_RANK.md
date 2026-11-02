# Near-total escape with actual-period full rank in one transverse dimension

2026-11-02. This gives an unconditional complete-packet construction for every sufficiently large height in a specified one-dimensional parameter region. It closes the global period-rank gap for that region, including every low-total-order row. It does **not** reach the critical high-dimensional parameter path, prove a common raw-column theorem for all heights, bound a nonzero period determinant from below, or establish noncancellation of the full exterior sum. No conclusion about the bad approximability of π follows.

## 1. Exact weighted space and complete source packet

Fix an integer K≥1, positive weights w₀,w₁,v₀,v₁ with v₁>v₀, and put c=v₁−v₀ and θ=w₁/v₁. Use exactly

\[
V_H=\operatorname{span}\{Y^hX^\alpha:h,\alpha\in\mathbf Z_{\ge0},\ w_0h+w_1\alpha\le H\}
\]

and the complete source packet

\[
L_H=\{(j,s,\beta):0\le j<K,\ s,\beta\ge0,\ v_0s+v_1\beta<H\}.
\]

For a center parameter Z define the coefficient-normalized formal-period rows by

\[
\mathcal R_{j,s,\beta}(Z)P=[z^su^\beta]P(e^z,jZ+z+u).
\tag{1}
\]

At Z=ω=2πi these are the actual-period rows after using e^{jω}=1. Z is a formal polynomial parameter elsewhere; no assertion e^Z=1 is made for arbitrary Z. The same arguments for complete packets work if the inner z in the X coordinate is replaced by any formal G(z) with G(0)=0. Thus the complete-source rank argument is unaffected by the original local logarithm-coordinate convention.

For a conserved total order N=s+β, the source types at each center form

\[
B_N=\{0,\ldots,b_N-1\},\qquad
b_N=\min\left(N+1,\left\lceil\frac{H-v_0N}{c}\right\rceil\right),
\quad 0\le N<H/v_0.
\tag{2}
\]

No low-N rows will be removed.

## 2. A full adjacent-interval transport selection

If 2b_N−1≤N, call N packed and set A_N={b_N,…,2b_N−1}. Otherwise set A_N=B_N. At every center j replace the rows in bucket N by

\[
\{(j,N-a,a):a\in A_N\}.
\tag{3}
\]

Call the resulting full target packet J_H. All time orders are nonnegative. Distinct buckets have distinct conserved orders, and center labels are retained, so all target rows are distinct. A packed bucket lies entirely outside L_H; the other buckets are unchanged.

For any integer b≥1 the adjacent Pascal block satisfies

\[
\det\left[\binom a\beta\right]_{0\le\beta<b,\ b\le a<2b}=1.
\tag{4}
\]

Indeed, the polynomials binom(x,β) have degree β and leading coefficient 1/β!. The determinant is the Vandermonde of the consecutive integers b,…,2b−1 divided by ∏_{β=0}^{b−1}β!, hence one. Identity blocks also have determinant one. The full normalized transport coefficient is consequently

\[
\boxed{\det T_{L_H,J_H}=(-\zeta)^{\Gamma_H},\qquad
\Gamma_H=K\sum_{N\text{ packed}}b_N^2.}
\tag{5}
\]

This is a nonzero coefficient whenever ζ≠0. It is the determinant after regrouping a complete bucket, not an assertion about an arbitrary individual matching. With factorial-scaled time rows the coefficient must additionally be multiplied by ∏_{L_H}s!/∏_{J_H}t!.

The number outside is exactly

\[
k_H=K\sum_{N\text{ packed}}b_N,\qquad M_H=|L_H|=K\sum_Nb_N.
\tag{6}
\]

## 3. An unconditional full-period-rank theorem

Assume auxiliary weights u₀,u₁ exist such that

\[
\boxed{w_0<u_0<v_0,\quad Kw_1<u_1<v_1,\quad
0<\Delta:=u_1-u_0<c/2.}
\tag{7}
\]

Let

\[
\kappa=\max\{u_0/v_0,2\Delta/c\}<1,
\qquad \eta=\max\{w_0/u_0,Kw_1/u_1\}<1.
\]

For

\[
H>\max\left\{\frac{\Delta}{1-\kappa},
\frac{w_1(K-1)}{1-\eta}\right\},
\tag{8}
\]

both L_H and J_H have full row rank on the original V_H at every nonzero center parameter Z, in particular Z=2πi.

### 3.1 Embed every selected row in one complete auxiliary packet

Consider

\[
U_H=\{(j,t,a):0\le j<K,\ t,a\ge0,\ u_0t+u_1a<H\}.
\]

The identity rows lie in U_H because u₀<v₀ and u₁<v₁. In a packed bucket the cap b_N=N+1 is impossible, so b_N=ceil((H−v₀N)/c). Every selected a therefore obeys

\[
a\le2b_N-1<2(H-v_0N)/c+1.
\]

Consequently

\[
\begin{aligned}
u_0(N-a)+u_1a
&=u_0N+\Delta a\\
&<\frac{u_0}{v_0}(v_0N)+\frac{2\Delta}{c}(H-v_0N)+\Delta\\
&\le\kappa H+\Delta<H.
\end{aligned}
\tag{9}
\]

Thus the **entire** target packet is a subset of U_H. This handles the low-N identity rows at the same time as the escaped rows.

### 3.2 A complete auxiliary packet has an explicit Hermite staircase

For each integer t≥0 with u₀t<H, let

\[
C_t=\left\lceil\frac{H-u_0t}{u_1}\right\rceil.
\]

Choose the K C_t columns

\[
Q_{t,\alpha}=(Y-1)^tX^\alpha,
\qquad 0\le\alpha<KC_t.
\tag{10}
\]

Their entire raw-monomial support obeys the original weighted budget, because

\[
w_0t+w_1\alpha
<H\max\{w_0/u_0,Kw_1/u_1\}+w_1(K-1)<H.
\tag{11}
\]

The row time r and column time t block is zero for r<t. For r=t it is

\[
\left[\binom\alpha a(jZ)^{\alpha-a}\right]
_{0\le j<K,\ 0\le a<C_t;\ 0\le\alpha<KC_t},
\]

the ordinary confluent Vandermonde for K distinct nodes jZ. Its determinant, in the usual center/type order, is

\[
\prod_{0\le i<j<K}((j-i)Z)^{C_t^2}\ne0.
\tag{12}
\]

Thus the entire square U_H matrix is block lower triangular with nonzero diagonal. Every subset of its rows is independent on V_H. In particular both L_H and J_H are independent. This uses only elementary one-dimensional Hermite interpolation; no unspecified large-weight threshold or general-position argument enters the proof. It also avoids a transcendence argument for the general full-row-rank assertion.

The auxiliary column set is downward closed in t at each α, since C_t is nonincreasing. Replacing Q by its supporting raw monomials is a determinant-one triangular column change for this full auxiliary square matrix. For a smaller subset of its rows, row independence implies some raw minor is nonzero, but does not identify the same minor for two different row packets and center parameters.

### 3.3 What common columns do and do not follow

The source matrix at a nonzero rational center r=2ip/q and the target matrix at ω each have rank M_H on the full V_H. Hence common bounded **mixed** columns exist: if X is a variable matrix with M_H columns, det(AX)det(BX) is a nonzero polynomial of degree at most two in each entry of X. A nonzero polynomial with that coordinate degree cannot vanish on the whole grid {0,1,2}^d, so some integer mixing matrix with entries 0,1,2 makes both determinants nonzero.

This argument does not give a common raw-column base and does not preserve a raw-minor arithmetic prefactor without redoing the corresponding mixed-column estimates. The finite certificates below do give common raw columns. A theorem giving them for every sufficiently large H remains open here.

## 4. Escape fraction and the precise displacement cost

Take H→∞ with all parameters fixed. In the scaled variable x=N/H, the limiting bucket profile is

\[
b(x)=\min\{x,(1-v_0x)/c\},\qquad 0\le x\le1/v_0.
\]

Packing holds beyond x*=2/(v₁+v₀), apart from an O(1) rounding shift in N. Riemann sums give

\[
M_H\sim\frac{KH^2}{2v_0v_1},\qquad
\Gamma_H\sim\frac{KcH^3}{3v_0(v_1+v_0)^3}.
\]

In particular

\[
\boxed{\lim_{H\to\infty}\frac{k_H}{M_H}
=\frac{v_1(v_1-v_0)}{(v_1+v_0)^2},}
\tag{13}
\]

\[
\boxed{\lim_{H\to\infty}\frac{w_1\Gamma_H}{M_HH}
=\frac{2\theta v_1^2(v_1-v_0)}{3(v_1+v_0)^3}.}
\tag{14}
\]

Similarly cΓ_H/(M_HH) tends to 2v₁(v₁−v₀)²/[3(v₁+v₀)³]. In the subsequent limit v₁/v₀→∞,

\[
k_H/M_H\longrightarrow1,\qquad
w_1\Gamma_H/(M_HH)\longrightarrow2\theta/3,
\qquad c\Gamma_H/(M_HH)\longrightarrow2/3.
\tag{15}
\]

If ν=(−log|ζ|)/w₁ is the scalar error-quality ratio, the normalized negative logarithm of (5), divided by νθ, tends to 2/3 in this successive limit. In particular, it tends to 2νθ/3 if ν and θ are held fixed; when they vary but have limits, their limiting values enter that expression. This is a statement about the **transport scalar**, not a lower bound on its accompanying period minor or on an aggregate coefficient.

At a fixed type a, at most one row occurs for each pair (j,N), so n_a≤K ceil(H/v₀). Therefore the collision score satisfies

\[
\frac{c_{\mathrm{col}}\sum_an_a^2}{M_HH}
\le\frac{c_{\mathrm{col}}K}{v_0}+o(1).
\tag{16}
\]

Provided ν and c_col stay bounded along the parameter family, the scalar-plus-occupancy cost can thus be small in this parameter region: choose θ small, v₀ large compared with K, and then v₁ much larger than v₀, always taking H→∞ first. This does not verify any claimed total-determinant lower bound.

## 5. Parameter scope and the time-grade obstruction

Conditions (7) are nonempty. For fixed K and θ<1/(2K), choose w₀<v₀ and then v₁=w₁/θ sufficiently large compared with v₀. One may choose u₀ between w₀ and v₀ and u₁ strictly between Kw₁ and c/2+u₀, also ensuring u₁<v₁ and u₁>u₀.

For the certificates take K=2, θ=1/8, w₀=1, v₀=2, v₁=8w₁ and

\[
u_0=3/2,\qquad u_1=3w_1.
\]

For w₁≥2 one has κ=3/4 and η=2/3, so the sufficient inequalities are H>12w₁−6 and H>3w₁. These are explicit, finite conditions.

The critical high-dimensional family in the main investigation instead has v₀/w₀→0, with θ=1/2 and K growing exponentially with dimension. It violates w₀<v₀ and θ<1/(2K). This construction cannot be cited as a solution in that regime.

There is also an exact reason the source's own minimal Hermite staircase does not solve the common-column issue. In shifted columns (Y−1)^hX^α, a target row of time order t vanishes whenever h>t. Therefore a nonzero target determinant requires

\[
\sum h\le\sum_{J_H}t=\sum_{L_H}s-\Gamma_H.
\tag{17}
\]

The source's own staircase has ∑h=∑_{L_H}s and is therefore incompatible with every nontrivial transport selection. The larger auxiliary packet provides sufficient additional low-time columns to certify row rank, but selecting a common raw square minor remains a separate task. The finite shifted certificates satisfy the necessary deficit Δ_time=∑s−∑h≥Γ_H exactly as required.

## 6. Finite common raw-column certificates

Run `python verify_near_total_period_rank.py`. The standard-library program retains every source row, checks the auxiliary embedding and every distinct adjacent Pascal block exactly over the rationals, and records the complete selected column lists in `near_total_verification.json`.

All determinants in this table are factorial-scaled time-row determinants modulo 65537, with i↦256. Period columns are specialized at formal Z=1; source columns use r=2ip/q.

| H | w₁ | p/q | Rows M | Outside | Γ | Raw target determinant | Raw source determinant |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 24 | 2 | 22/7 | 32 | 22 | 30 | 6974 | 5673 |
| 40 | 2 | 22/7 | 72 | 48 | 88 | 16233 | 54224 |
| 60 | 5 | 355/113 | 80 | 70 | 102 | 60022 | 10964 |

The same lexicographic pivot index sets also work in shifted columns. Their source residues are respectively 38901, 54224, 28917; the period residues are unchanged. Their time deficits are 44,144,254, respectively, all at least Γ.

For each listed raw family the scaled formal-period determinant is a polynomial P(Z)∈Z[Z]. Nonzero reduction at Z=1 proves P is not identically zero. The known transcendence of 2πi then proves P(2πi)≠0 on that **particular same raw-column family**. For the rational source, reduction is from the localization Z[i,1/q] to the finite field, not from the entire field Q(i). Its nonzero residue proves characteristic-zero nonvanishing.

The entries 22/7 and 355/113 are established approximations used to give genuine rational-center certificates. No new rational approximation or record is claimed.

## 7. What changed

The earlier local packing theorem gave a nonzero transport minor without guaranteeing global period rank; an explicit ten-row counterexample showed that implication could fail. The construction here supplies a full global rank certificate for a nontrivial infinite-height family by putting **all** selected rows inside one explicit auxiliary Hermite packet. Nearly all rows can escape in a successive large-transverse-weight limit, while the transport and occupancy costs stay small. The unresolved parts are the critical high-dimensional regime, a general common raw-column theorem, quantitative period-minor size, and cancellation among different target selections.
