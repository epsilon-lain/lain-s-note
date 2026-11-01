# Coalescence order, genuine cancellation, and the small-error distinction

2026-11-01. This is an exact algebraic supplement to the sparse-shell checkpoint. It establishes neither bad approximability of π nor a new rational approximation. The explicit staircase Vandermonde minor was already established in the earlier nonzero-determinant report; the additional points here concern **every** one-dimensional minor, extra cancellation, and the difference between center coalescence and a small approximation error.

## 1. A universal order for every one-dimensional minor

Fix $K\ge2$. At each time order $s$, let the transverse derivatives be $\beta=0,\ldots,B_s-1$, where only finitely many nonnegative integers $B_s$ are nonzero. Use every row $(s,j,\beta)$, $0\le j<K$, and any equally sized family of polynomial columns $P\in\mathbb Q[Y,X]$. Order the rows by $s,j,\beta$. Let $G(t)\in t\mathbb Q[[t]]$, and define

\[
 A_{s,j,\beta;P}(r)=[t^su^\beta]P(1+t,jr+G(t)+u),
 \qquad D(r)=\det A(r).
\]

This is a polynomial in $r$. Put

\[
 W=\frac{K(K-1)}2\sum_s B_s^2,
 \qquad
 V_{s,a;P}=[t^su^a]P(1+t,G(t)+u),\quad 0\le a<KB_s.
\]

Then

\[
 r^W\mid D(r),\qquad
 [r^W]D(r)=
 \left(\prod_s\prod_{0\le i<j<K}(j-i)^{B_s^2}\right)\det V. \tag{1}
\]

In particular the order is exactly $W$ if and only if the square collapsed-row matrix $V$ is nonsingular. This is a criterion involving actual coefficients, rather than the existence of a support matching.

**Proof.** At a fixed $s$, the exact translation matrix from the coefficient rows indexed by $a$ to the rows indexed by $(j,\beta)$ has entries

\[
 C_{(j,\beta),a}(r)=\binom a\beta(jr)^{a-\beta},
\]

with zero entries for $a<\beta$. Expand its exterior product. A contribution with a repeated $a$ vanishes. For a distinct set $A_s$ of $KB_s$ nonnegative integers, its scalar determinant is homogeneous of degree

\[
 \sum_{a\in A_s}a-K\sum_{\beta=0}^{B_s-1}\beta
 \ge \frac{K(K-1)}2B_s^2.
\]

The unique set attaining equality is $A_s=\{0,\ldots,KB_s-1\}$. Its coefficient determinant is the normalized confluent Vandermonde determinant

\[
 \det\left(\binom a\beta j^{a-\beta}\right)
 =\prod_{i<j}(j-i)^{B_s^2}.
\]

For completeness, the last identity follows by considering independent nodes $c_j$. Coalescing any two nodes forces at least $B_s^2$ vanishing orders: after expanding the second node around the first, the new $B_s$ derivative degrees must be distinct and at least $B_s$. Their minimum extra total degree is $B_s^2$. Thus every $(c_j-c_i)^{B_s^2}$ divides the determinant. The total degree equals the sum of these exponents. Comparing the coefficient obtained successively from the highest available powers of the last node gives multiplicative constant one. Dividing derivatives by $\beta!$ is essential for this normalization. Multiplying the fixed-$s$ exterior expansions proves (1). ∎

The same argument with independent centers proves

\[
 \prod_{i<j}(c_j-c_i)^{\sum_sB_s^2}\mid D(c_0,\ldots,c_{K-1}).
\]

For the old staircase columns $Y^hX^a$, $h=s,0\le a<KB_s$, with nonincreasing $B_s$, the collapsed determinant is one and the total polynomial degree is $W$. This recovers the already-known exact Vandermonde minor; it is not a new high-dimensional construction.

## 2. Individually nonzero leading terms can genuinely cancel

Take a complete source-format packet

\[
 m=1,\ K=2,\ H=6,\ v_0=2,\ w_0=w_1=1,\ \theta=1/6.
\]

The strict row inequality $2s+6\beta<6$ gives $\beta=0$ and $s=0,1,2$ at each center. Take $G(t)=\ln(1+t)$ (only coefficients through degree two matter) and the six columns

\[
 X,\quad Y,\quad Y^2X,\quad Y^2X^2,\quad Y^2X^4,\quad Y^3X^2.
\]

Each column satisfies $h+\alpha\le H$. Direct exact calculation gives

\[
 \boxed{D(r)=5r^7.} \tag{2}
\]

Keep the center-zero rows fixed and translate the three center-$r$ rows to zero. If their new transverse indices are $(a_0,a_1,a_2)$, their scalar factor is $r^{a_0+a_1+a_2}$. Among all $5^3=125$ assignments $0\le a_s\le4$, exactly eleven coefficient determinants are nonzero. The lowest individually nonzero degree is six, with precisely these terms:

| New types $(a_0,a_1,a_2)$ | Coefficient determinant |
|---|---:|
| $(2,2,2)$ | $-12$ |
| $(2,3,1)$ | $16$ |
| $(4,1,1)$ | $-4$ |

They sum to zero. Grouping all eleven terms by degree gives

\[
 [r^6]D=-12+16-4=0,\qquad [r^7]D=5,\qquad [r^8]D=-2+8-6=0.
\]

Thus the lowest nonzero **individual** term degree is six, while the order of the full determinant is seven. This is an explicit demonstration that certifying individually nonzero low-cost terms does not rule out a stronger estimate after summation. The dimension and constants here are finite; no unknown separated-weight threshold is asserted, and this example is not an asymptotic endpoint estimate.

## 3. Coalescence does not supply a power of the approximation error

The preceding variable $r$ is the actual separation between centers. It is small when the centers merge. In the π problem the corresponding center is $r=2ip/q\approx\omega=2\pi i$, so the small quantity is instead $\varepsilon=r-\omega$.

For **any** nonzero one-dimensional raw polynomial minor $D(r)\in\mathbb Q[r]$, known transcendence of π gives

\[
 D(\omega)\ne0.
\]

Consequently $D(\omega+\varepsilon)$ has a nonzero constant term. For the exact example (2),

\[
 D(\omega+\varepsilon)=5\omega^7+35\omega^6\varepsilon+\cdots.
\]

There is no positive power of $\varepsilon$ dividing this expansion. The same statement holds over $\mathbb Q(i)$. This is an obstruction to a **universal exact vanishing factor** in the approximation error; it is not a quantitative lower bound for $D(\omega)$, and does not exclude rapid decay of that nonzero constant as the cutoff grows.

For several coordinates a nonzero polynomial $D(r_1,\ldots,r_m)$ may vanish identically on the diagonal. The analogous conclusion requires the additional certificate $D(Z,\ldots,Z)\not\equiv0$. Neither the one-variable argument nor the shell construction grants that certificate for an arbitrary fixed raw minor.

## 4. A general exact rank test for extra order

This diagnostic applies to any square polynomial matrix $A(r)=\sum_{k\ge0}A_kr^k$ over a characteristic-zero field, with nonzero determinant. Let $n$ be its size. For an integer $t\ge0$, form the block lower-triangular matrix

\[
 \mathcal T_t=(A_{i-j})_{0\le i,j\le t},
 \qquad A_k=0\quad(k<0).
\]

It represents multiplication by (A(r)) modulo $r^{t+1}$. There are nonnegative integers $e_1,\ldots,e_n$ such that

\[
 \operatorname{rank}\mathcal T_t=\sum_{i=1}^n\max(0,t+1-e_i),
 \qquad \operatorname{ord}_{r=0}\det A=\sum_i e_i. \tag{3}
\]

To prove this, perform Gaussian elimination over the formal power-series ring. Choose an entry with least order, divide its unit factor, and clear its row and column; all needed divisions remain in the ring because every other entry has at least this order. Repeat on the remaining block. This yields diagonal entries $r^{e_i}$ under invertible row and column operations. The same operations are invertible modulo every $r^{t+1}$. Multiplication by $r^{e_i}$ there has rank $\max(0,t+1-e_i)$, proving (3).

Thus the integer deficiencies

\[
 n(t+1)-\operatorname{rank}\mathcal T_t
 =\sum_i\min(t+1,e_i)
\]

eventually stabilize at the exact determinant order. This avoids treating a support matching or a single rank at $r=0$ as a complete valuation certificate.

An exact four-row example uses $s=0,1$, $\beta=0$, $K=2$, columns $1,X,X^3,YX^2$. It is a complete source-format packet with $H=4,v_0=2,w_0=w_1=1,\theta=1/4$. Its matrix and determinant are

\[
 A(r)=\begin{pmatrix}
 1&0&0&0\\
 1&r&r^3&r^2\\
 0&1&0&0\\
 0&1&3r^2&r^2+2r
 \end{pmatrix},\qquad D(r)=r^4(1-r).
\]

The universal order (1) is two, but the actual order is four. The successive ranks for $t=0,1,2,3,4$ are $2,5,8,12,16$; the deficiencies are $2,3,4,4,4$. The local invariant orders are $0,0,1,3$. This certifies the additional degeneracy exactly.

## Reproduction and scope

Run `python3 verify.py`. The standard-library verifier reconstructs determinant polynomials by exact rational interpolation using proved degree bounds, checks them at additional rational values, checks the leading coefficient formula for four named examples and 32 additional minors, exhausts all 125 assignments in (2), and computes the block ranks in (3). Results are in `verification.json`.

The new useful conclusion is two-sided: real cancellation between individually nonzero terms is present, but coalescence order alone does not produce a small-error factor near $2\pi i$. Any endpoint use still needs a quantitative estimate at the actual period, uniformly through the required asymptotic range.
