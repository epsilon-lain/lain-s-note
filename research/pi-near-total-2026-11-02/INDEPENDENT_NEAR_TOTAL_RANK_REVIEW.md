# Independent audit of the near-total period-rank construction

The complete-packet construction, auxiliary Hermite rank proof, escape fraction, and displacement-cost formulas are correct. This review independently checks the points that can otherwise conceal a low-order or limit-order gap. The companion note has incorporated the qualification that the error-quality and occupancy coefficients must be controlled when its parameter family varies.

## 1. Exact buckets and transport

For fixed total order \(N=s+\beta\), the source condition is

\[
0\le\beta\le N,\qquad (v_1-v_0)\beta<H-v_0N.
\]

Thus the stated count

\[
b_N=\min\{N+1,\lceil(H-v_0N)/(v_1-v_0)\rceil\}
\]

correctly handles both strict inequalities and integer endpoints. Replacing the interval \(0,\ldots,b_N-1\) by \(b_N,\ldots,2b_N-1\) is possible exactly under the stated nonnegative-time condition \(2b_N-1\le N\). If that condition fails, retaining the source bucket preserves every low-order row.

The adjacent Pascal determinant is one, including its sign with both index intervals increasing. Its shear exponent is

\[
\sum_{a=b}^{2b-1}a-\sum_{\beta=0}^{b-1}\beta=b^2.
\]

Consequently \(\Gamma_H=K\sum_{\mathrm{packed}}b_N^2\) and \(k_H=K\sum_{\mathrm{packed}}b_N\) are exact. A packed target is genuinely outside the source packet because its spatial type is at least \(b_N\) at the same conserved order and center. The source and target row counts agree and there are no collisions across buckets or centers.

## 2. Continuous profile and the order of limits

Write \(c=v_1-v_0>0\). With all weights and \(K\) fixed, the limiting profile is

\[
b(x)=\begin{cases}
x,&0\le x\le1/v_1,\\
(1-v_0x)/c,&1/v_1\le x\le1/v_0.
\end{cases}
\]

The nonzero packed part begins at \(x_*=2/(v_1+v_0)\). This is strictly beyond the corner \(1/v_1\). The discrete ceiling and parity affect the threshold by a bounded number of indices for fixed weights. They contribute only lower-order terms to the three sums.

Direct integration gives

\[
\int_0^{1/v_0}b(x)\,dx=\frac1{2v_0v_1},
\]

\[
\int_{x_*}^{1/v_0}b(x)\,dx
=\frac{c}{2v_0(v_1+v_0)^2},
\qquad
\int_{x_*}^{1/v_0}b(x)^2\,dx
=\frac{c}{3v_0(v_1+v_0)^3}.
\]

These establish all three leading sums, including the stated coefficient of \(\Gamma_H\). Dividing by the source mass yields

\[
\lim_{H\to\infty}\frac{k_H}{M_H}
=\frac{v_1c}{(v_1+v_0)^2},
\]

\[
\lim_{H\to\infty}\frac{w_1\Gamma_H}{M_HH}
=\frac{2\theta v_1^2c}{3(v_1+v_0)^3},
\qquad
\lim_{H\to\infty}\frac{c\Gamma_H}{M_HH}
=\frac{2v_1c^2}{3(v_1+v_0)^3}.
\]

Only after this fixed-parameter height limit may one let \(v_1/v_0\to\infty\), obtaining \(1\), \(2\theta/3\) (with \(\theta\) fixed), and \(2/3\), respectively. A simultaneous limit with varying weights is not proved by these Riemann sums. For example, taking \(v_1\) larger than the height changes the bucket scale and cannot be justified by the displayed fixed-weight asymptotics.

For scalar transport cost, the invariant formulation is

\[
\frac{1}{\nu\theta}
\lim_{H\to\infty}
\frac{-\log|\det T_{L_H,J_H}|}{M_HH}
=\frac{2v_1^2c}{3(v_1+v_0)^3}
\longrightarrow\frac23,
\]

where division assumes \(\nu\theta\ne0\). A literal limit \(2\nu\theta/3\) requires those factors to be held fixed or to have the specified limiting values. The occupancy bound follows from

\[
\sum_an_a^2\le(\max_an_a)\sum_an_a
\le K\lceil H/v_0\rceil M_H.
\]

The resulting small-cost parameter choice requires bounded or otherwise controlled \(\nu\) and \(c_{\mathrm{col}}\). The revised companion note makes this explicit. Neither scalar nor occupancy cost controls the accompanying period-minor magnitude or cancellation in the complete exterior sum.

## 3. Auxiliary containment includes all low-order rows

The conditions

\[
w_0<u_0<v_0,\quad Kw_1<u_1<v_1,\quad
0<\Delta=u_1-u_0<c/2
\]

are sufficient. Unchanged rows belong to the auxiliary packet by coordinatewise weight comparison. In a packed bucket, the cap \(b_N=N+1\) is impossible, and

\[
a\le2b_N-1<2(H-v_0N)/c+1.
\]

Thus

\[
u_0(N-a)+u_1a
<\max(u_0/v_0,2\Delta/c)H+\Delta<H
\]

under the stated height condition. This is a strict inequality, so no omitted boundary case remains.

For auxiliary time \(t\), exactly \(C_t=\lceil(H-u_0t)/u_1\rceil\) spatial types occur at each center. The proposed \(KC_t\) shifted columns satisfy

\[
w_0t+w_1\alpha
<H\max(w_0/u_0,Kw_1/u_1)+w_1(K-1)<H.
\]

Their raw expansions also satisfy the budget because lowering the power of \(Y\) lowers its cost. Ordinary-coefficient jets give a lower time-block triangular matrix with the stated confluent Vandermonde diagonal blocks. The coefficient \([z^t](e^z-1)^t\) is exactly one, so no missing time factorial occurs. The diagonal blocks are invertible for nonzero \(Z\) (for \(K=1\) even the nonzero-center condition is unnecessary). This proves independence of the entire auxiliary packet and every row subset on the original weighted space.

The proof is valid with \(z\) inside the spatial argument replaced by any formal series \(G(z)\) satisfying \(G(0)=0\): the time support and diagonal blocks are unchanged. This observation justifies the stated rank invariance; it is not by itself a derivation of a different shear transport formula in those coordinates.

## 4. Common-column scope

Separate full row rank does not imply a common raw-column basis. The bounded mixed-column argument is valid: the product of the two determinant polynomials is nonzero over the complex coefficient field, has degree at most two in each mixing variable, and therefore is nonzero somewhere on the Cartesian grid \(\{0,1,2\}^d\).

That argument proves existence of common integer mixed columns. It does not retain a particular raw-column arithmetic prefactor without a new estimate. The source's minimal shifted Hermite staircase has zero time surplus and is correctly excluded for the nontrivial target by the time-grade condition.

The parameter region is nonempty when \(K\theta<1/2\), \(w_0<v_0\), and \(v_1/v_0\) is sufficiently large. It remains separate from the investigation's critical high-dimensional regime. No fixed-dimension proof here supplies a uniform estimate as dimension grows.

## 5. Independent finite verification

The three saved certificates were recomputed independently modulo \(65537\), using direct binomial expansion of \((Y-1)^h\) rather than the script's Stirling-number recurrence. Both the raw and shifted determinants matched every saved residue:

| Height | Basis | Target residue | Rational-source residue |
|---:|---|---:|---:|
| 24 | Raw | 6974 | 5673 |
| 24 | Shifted | 6974 | 38901 |
| 40 | Raw | 16233 | 54224 |
| 40 | Shifted | 16233 | 54224 |
| 60 | Raw | 60022 | 10964 |
| 60 | Shifted | 60022 | 28917 |

The use of \(256^2\equiv-1\pmod{65537}\) and denominator localization is valid for the displayed rational centers. The formal-period determinant after factorial time scaling is an integer polynomial in \(Z\). Its nonzero residue at \(Z=1\) proves it is a nonzero polynomial, and hence it cannot vanish at the transcendental actual center \(2\pi i\). This supports only the specific listed common raw and shifted families; it does not replace the remaining general common-raw-column argument.

**Audit outcome:** the construction and asymptotic formulas pass, with the fixed-parameter limit order and coefficient-control qualifications above. The rank theorem is unconditional in its stated one-dimensional region; quantitative determinant size, total-sum noncancellation, and the critical regime remain open in this work.
