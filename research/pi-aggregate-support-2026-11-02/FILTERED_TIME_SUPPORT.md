# Positive time surplus: cut constraints, intrinsic column filtration, and marker degree

2026-11-02. The total inequality \(|\Gamma|\le\Delta_{\rm time}\) is the sum of a stronger family of necessary inequalities. Those inequalities also give an intrinsic formulation for arbitrary selected raw columns, a rank bound on upward feedback, and a sharper degree bound for a marker counting outside rows. Two genuine period-jet examples below distinguish the statements. These results give support and exact rank information, not a quantitative lower bound on a determinant or a Diophantine conclusion.

## 1. The time-cut constraints

Let \(L\) be the complete source packet, with source row time \(s\), and fix a square family of shifted columns with time orders \(h\). For each positive integer \(r\), set

\[
S_r=\#\{\ell\in L:s(\ell)\ge r\},\quad
d_r=\#\{\text{columns}:h\ge r\},\quad C_r=S_r-d_r.
\]

For a target set \(J\), put

\[
E_r(J)=S_r-\#\{j\in J:t(j)\ge r\}.
\]

If both its complete transport minor and its period minor on the fixed columns are nonzero, then

\[
\boxed{0\le E_r(J)\le C_r\quad\text{for every }r\ge1.}
\tag{1}
\]

Indeed a nonzero transport permutation matches each source time \(s\) to a target time \(t\le s\). Its contribution to \(E_r\) is one precisely when \(t<r\le s\). A nonzero period permutation matches column time \(h\) to row time \(t\ge h\), so the number of columns at or above a cut cannot exceed the number of target rows there.

Equivalently, the target times and column times sorted in increasing order must obey \(t_{(i)}\ge h_{(i)}\) for every index. These are the Hall conditions for this nested zero pattern. They guarantee a matching for the time zero pattern alone, not nonvanishing of the actual period determinant.

Summing (1) gives the previous bound, since

\[
\sum_{r\ge1}E_r=\sum_Ls-\sum_Jt=|\Gamma|,
\qquad \sum_{r\ge1}C_r=\Delta_{\rm time}.
\]

A cut with \(C_r=0\) cannot be crossed by even one edge of a nonzero transport matching contributing a nonzero period minor. A positive budget at a different cut cannot pay for such a crossing.

## 2. An intrinsic formulation that applies to selected raw columns

Let \(W\) be the \(M\)-dimensional span of the chosen polynomial columns, and replace the definition of \(d_r\) by

\[
\boxed{d_r=\dim\bigl(W\cap (Y-1)^r\mathbb C[Y,X]\bigr).}
\tag{2}
\]

Every target row with \(t<r\) annihilates that subspace. If the square period evaluation on \(W\) is invertible, its restriction to this \(d_r\)-dimensional subspace is injective and supported on rows with \(t\ge r\). Therefore (1) remains necessary. For shifted monomial columns, (2) reduces to the preceding count exactly.

An adapted basis of the finite filtration has time orders \(h_1,\ldots,h_M\) with \(d_r=\#\{i:h_i\ge r\}\). The effective total surplus is therefore

\[
\Delta_W=\sum_Ls-\sum_{r\ge1}d_r.
\tag{3}
\]

This basis change is used to describe support; it need not be an integral determinant-one transformation and is not an untracked change in an arithmetic lower bound.

There is an explicit formula for a selected **raw monomial** family. For each spatial multi-index \(\alpha\), let \(k_\alpha\) be the number of distinct selected exponents of \(Y\) paired with \(X^\alpha\). Then

\[
\boxed{d_r=\sum_\alpha\max(k_\alpha-r,0),\qquad
\Delta_W=\sum_Ls-\sum_\alpha\binom{k_\alpha}{2}.}
\tag{4}
\]

To prove this, at a fixed \(\alpha\) with selected powers \(h_1,\ldots,h_k\), the first \(k\) coefficient jets at \(Y=1\) have matrix \(\binom{h_j}{i}\), \(0\le i<k\). Its determinant is

\[
\frac{\prod_{a<b}(h_b-h_a)}{\prod_{i=0}^{k-1}i!}\ne0
\]

after putting the distinct powers in increasing order. Thus the vanishing orders of an adapted basis are exactly \(0,1,\ldots,k-1\), regardless of gaps in the selected raw powers. Distinct \(\alpha\) blocks remain independent as polynomials in \(X\). This yields (4), while preserving each block's spatial degree and the original column span. The weighted budget is preserved because these are combinations of already selected columns.

This avoids the invalid replacement of an arbitrary raw subset by the same-index shifted subset. Their spans and their time filtrations can differ.

## 3. A bound on feedback across each time cut

Assume the source evaluation \(B:W\to\mathbb C^L\) is invertible. Let \(C=B_O\), \(X=CB^{-1}\), and let \(R\) be the exact boundary transport, whose entries vanish unless outside time \(t\) is strictly smaller than source time \(s\). Set \(N=RX\).

Write \(O_{<r}\) for outside rows of time less than \(r\), and \(L_{\ge r}\) for source coordinates of time at least \(r\). Then

\[
\boxed{\operatorname{rank}X_{O_{<r},L_{\ge r}}\le C_r.}
\tag{5}
\]

The inverse image under \(B\) of the source-coordinate subspace supported on \(L_{\ge r}\) has dimension \(S_r\). It contains the subspace in (2), of dimension \(d_r\), and every row of \(C\) indexed by \(O_{<r}\) annihilates that subspace. Hence this restricted map factors through a quotient of dimension \(C_r\).

Since rows of \(R\) with source time at most \(r\) use only outside times less than \(r\),

\[
\boxed{\operatorname{rank}N_{L_{\le r},L_{\ge r}}\le C_r.}
\tag{6}
\]

Thus a small surplus at a cut bounds the rank of all feedback across that cut, including the time-\(r\) diagonal block. Zero-surplus cuts give exact triangular separation. When every \(C_r=0\), this recovers strict lower triangularity and nilpotence. For positive surplus, it does not imply nilpotence or a norm bound.

## 4. The boundary marker can have much lower degree than total displacement

Define \(M_s=\#\{\ell\in L:s(\ell)=s\}\), and let \(O_t\) be the number of available outside rows of time \(t\). In the exact boundary-marker determinant

\[
F(\eta)=\det(B+\eta RC),
\]

Cauchy–Binet associates \(\eta^k\) with selecting \(k\) outside rows. For any nonzero such contribution, choose a nonzero permutation in the corresponding \(R\) minor. Let \(e_s\) count its omitted source rows of time \(s\). Each such match strictly lowers time and therefore crosses cut \(s\). It follows that

\[
e_s\le M_s,\qquad e_s\le E_s\le C_s.
\]

Summing proves

\[
\boxed{\deg_\eta F\le\rho_{\rm source}:=
\sum_{s\ge1}\min(M_s,C_s)\le\Delta_W.}
\tag{7}
\]

Similarly, every selected outside row of time \(t\) is reached by a match crossing cut \(t+1\), so

\[
\deg_\eta F\le\rho_{\rm outside}:=
\sum_{t\ge0}\min(O_t,C_{t+1}).
\tag{8}
\]

These are upper bounds; actual coefficients can vanish. The usual rank bounds can further reduce them. In multivariate shear variables a contributing term still obeys \(k\le|\Gamma|\le\Delta_W\). The marker counts escaped rows, whereas the shear degree counts their total displacement, so these two degrees need not agree.

There is a limitation in the critical parameter range. Every polynomial in the allowed weighted space has \(Y\)-degree at most \(\lfloor H/w_0\rfloor\), so \(d_s=0\) for \(s>H/w_0\). At each such source time, \(C_s=S_s\ge M_s\), and consequently

\[
\rho_{\rm source}\ge\#\{\ell\in L:s(\ell)>H/w_0\}.
\tag{9}
\]

For fixed dimension \(m\) and fixed weights, the fraction on the right tends, as \(H\to\infty\), to \((1-v_0/w_0)_+^{m+1}\). This follows by translating the time coordinate in the source simplex; the remaining simplex has height \(H(1-v_0/w_0)\). On the critical family \(v_0/w_0=2K\theta^m\), \(\theta=1/2\), \(K=\lfloor1.3^m\rfloor\), this fraction then tends to one as \(m\to\infty\). Thus this particular marker-degree upper bound is itself asymptotically almost \(M\) there, regardless of the chosen allowed column space. It cannot alone force a small actual marker degree in that regime. This is a statement about the limitation of the bound, not a lower bound on the actual determinant degree. The height limit is taken first for each fixed parameter member.

## 5. Genuine six-row example: the total test passes but a cut kills the target

Take \(K=2\), \(H=3\), \(v_0=1\), \(v_1=4\), \(w_0=1/4\), \(w_1=1\), \(\theta=1/4\). The complete source consists of \((j,s,0)\) with \(j=0,1\), \(s=0,1,2\). Select the valid shifted columns

\[
1,\ X,\ X^2,\ (Y-1),\ (Y-1)^2,\ (Y-1)^2X.
\]

At formal center difference \(Z\), with increasing-time source order, the source determinant is \(-2Z^3\), nonzero in particular at \(Z=2\pi i\). Here

\[
C_1=1,\qquad C_2=0,\qquad \Delta_W=1.
\]

Replace only the source row \((0,2,0)\) by \((0,1,1)\). Its transport determinant is \(-\zeta\), and \(|\Gamma|=1=\Delta_W\), so the total test passes. Nevertheless this replacement crosses cut \(2\), where \(E_2=1>C_2=0\), and its period determinant is zero.

By contrast, replacing \((0,1,0)\) by \((0,0,1)\) crosses only cut \(1\). Its period determinant is \(-Z^3\), so this target survives. Among all 36 transport-supported target selections, five pass the total-surplus test, three pass all cut constraints, and precisely those three have nonzero period determinant in this example.

The relative boundary operator has rank two, yet its exact marker determinant is

\[
\det(I+\eta N)=1-\eta\zeta.
\]

Thus positive surplus permits nontrivial determinant change, while the location of the surplus controls which target can cause it.

## 6. Genuine family where the marker-degree improvement is sharp

Suppose the complete source has only spatial type \(\beta=0\), with times \(s=0,\ldots,S\) at each of \(K\) centers. Choose \(M=K(S+1)\) columns \(X^\alpha\), \(0\le\alpha<M\), all with \(h=0\), assuming their weighted budget fits. The source matrix is the ordinary confluent Vandermonde and is invertible at distinct centers.

For an outside row \((j,t,a)\) with \(t+a=s\), its evaluation on these columns is \(\binom{s}{a}\) times source row \((j,s,0)\). Since the inside transport is the identity in this single-type source, the relative boundary operator is diagonal, with entry

\[
(1-\zeta)^s-1
\]

at each \((j,s,0)\). Consequently

\[
\boxed{\det(I+\eta N)=
\prod_{s=1}^S\bigl[1+\eta((1-\zeta)^s-1)\bigr]^K.}
\tag{10}
\]

Its marker degree is \(KS=\rho_{\rm source}\), while its shear degree and total surplus are \(KS(S+1)/2\). Thus the improvement in (7) is sharp in genuine period matrices and can be a factor \((S+1)/2\).

A finite weighted instance is \(K=2\), \(S=3\), \(H=4\), \(v_0=1\), \(v_1=5\), \(w_1=1/2\), \(\theta=1/10\). The eight columns are \(1,X,\ldots,X^7\), their costs are at most \(7/2<4\), and the complete source indeed has only \(\beta=0\). The source determinant in increasing-time order is \(Z^{16}\); the marker degree is six although the total surplus is twelve.

The finite verification program uses symbolic \(Z\), not an assumption that arbitrary \(Z\) is an exponential period. Substitution \(Z=2\pi i\) yields the genuine period examples. None of these exact support statements gives the missing quantitative size estimate in the critical parameter family.
