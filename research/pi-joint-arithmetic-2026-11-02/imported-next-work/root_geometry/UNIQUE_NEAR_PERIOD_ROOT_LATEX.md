# A unique simple imaginary zero at the Newton scale

Research derivation, 2026-11-02. This strengthens the nearest-root result for
the pinned determinant family. It gives no rational approximation theorem.

Let \(t=2\pi\), \(\omega=it\),
\[
A_j=\prod_{a=0}^{j-1}a!,\qquad
\Delta_j=\frac{A_j^4}{A_{2j}},\qquad A_0=\Delta_0=1,
\]
and use the repository's moment matrices
\[
H_n(Z)_{ij}=\int_0^1u^{i+j}e^{Zu}\,du,\qquad
V_n(Z)_{ij}=\frac1Z\int_0^\infty e^{-v}(1-v/Z)^{i+j}\,dv.
\]
Put \(D_n(Z,\eta)=\det(H_n(Z)+\eta V_n(Z))\). The exact endpoint identity is
\[
P_n(Z)=\frac{Z^{n^2}}{A_n^2}D_n(Z,1-e^Z).
\tag{1}
\]
The established inputs are
\[
|\det H_n(\omega)|>\Delta_n/6,
\qquad
T_n:=\operatorname{tr}(H_n(\omega)^{-1}V_n(\omega))
=i(-1)^nK_n\frac{\chi_{n-1}(t)}{\chi_n(t)}
 (1+O(n^{-1/2})),
\tag{2}
\]
\[
K_n=\frac{[(2n-2)!]^2(2n-1)!}{[(n-1)!]^4t^{2n-1}},
\qquad 1/6<\chi_j(t)\le1\quad(j\ge2).
\tag{3}
\]
In particular \(K_n/12\le|T_n|\le12K_n\) eventually, and \(K_n\) grows
faster than every fixed power of \(n\).

**Theorem.** For all sufficiently large integers \(n\), exactly one zero of
\(P_n\), counted with algebraic multiplicity, lies in
\(\lvert Z-\omega\rvert\le48/K_n\). This zero \(\zeta_n\) is simple,
purely imaginary, and is the nearest zero to \(\omega\). Moreover
\[
\boxed{\zeta_n-\omega
=-\frac{P_n(\omega)}{P_n'(\omega)}(1+O(n^{-2}))}
\tag{4}
\]
and
\[
\boxed{\log|\zeta_n-\omega|
=-2n\log n-\bigl[\log(16/\pi^2)-2\bigr]n+O(\log n).}
\tag{5}
\]
It lies below \(\omega\) for even \(n\), and above it for odd \(n\).

## 1. A uniform majorant for every endpoint order

For \(0\le k\le n\), set \(m=n-k\) and
\[
B_{n,k}(s)=\frac{\Delta_m}{\Delta_n}A_k
 \prod_{j=0}^{k-1}(2m+j)!\;s^{-k(2n-k)}.
\tag{6}
\]
Write
\[
\frac{D_n(Z,\eta)}{\det H_n(\omega)}
=\sum_{k=0}^n C_{n,k}(Z)\eta^k.
\]
Then, for \(Z\ne0\),
\[
|C_{n,k}(Z)|\le
6e^{m\max(\Re Z,0)+k|Z|}B_{n,k}(|Z|).
\tag{7}
\]
To prove this, expand the Vandermonde-square determinant integral with
\(m\) ordinary variables \(u_j\in[0,1]\) and \(k\) endpoint variables
\(x_a=1-v_a/Z\), with factor \(1/(m!k!)\). The physical integral is bounded
by \(m!\Delta_m e^{m\max(\Re Z,0)}\). Each cross term obeys
\[
|x_a-u_j|\le(v_a+|Z|)/|Z|.
\]
The endpoint Vandermonde contributes
\(\prod_{a<b}(v_b-v_a)^2/|Z|^{k(k-1)}\); the endpoint measures contribute
\(|Z|^{-k}\). Thus the total denominator exponent is
\(2mk+k(k-1)+k=k(2n-k)\).
Shift \(y_a=v_a+|Z|\), then enlarge the integration domain to the positive
orthant. The resulting nonnegative integral is
\[
\frac1{k!}\int_{[0,\infty)^k}e^{-\sum y_a}
 \prod_a y_a^{2m}\prod_{a<b}(y_b-y_a)^2\,dy
=A_k\prod_{j=0}^{k-1}(2m+j)!.
\tag{8}
\]
For completeness, the monic Rodrigues polynomial
\[
L_j(y)=(-1)^je^yy^{-2m}\frac{d^j}{dy^j}
 (e^{-y}y^{j+2m})
\]
is orthogonal to smaller degrees under \(y^{2m}e^{-y}dy\), by integration
by parts, and has squared norm \(j!(j+2m)!\). The factorial moment Hankel
determinant is their norm product, proving (8). This argument includes
\(m=0\); empty products handle \(k=0\).

The exact ratios are
\[
B_{n,0}(t)=1,\quad B_{n,1}(t)=K_n,\qquad
\frac{B_{n,k+1}(t)}{B_{n,k}(t)}
=\frac{k!K_{n-k}}{(2n-2k)_k},
\tag{9}
\]
where \((a)_k=\prod_{j=0}^{k-1}(a+j)\). Since \((2m)_k\ge k!\), and
\[
\frac{K_{m+1}}{K_m}
=\frac{8(2m-1)^2(2m+1)}{mt^2}>1\quad(m\ge2),
\]
we have \(\max_{1\le m<n}K_m=K_{n-1}\) for all sufficiently large \(n\).
Consequently
\[
B_{n,k}(t)\le K_nK_{n-1}^{k-1}\quad(1\le k\le n),
\qquad
\frac{K_{n-1}}{K_n}
=\frac{(n-1)t^2}{8(2n-3)^2(2n-1)}=O(n^{-2}).
\tag{10}
\]

## 2. A locally uniform linear limit

On \(|Z-\omega|\le1/n\), (7) implies
\[
|C_{n,k}(Z)|\le6e c^kB_{n,k}(t),\qquad c=e^{t+1}.
\tag{11}
\]
Indeed \(m\max(\Re Z,0)\le1\), \(|Z|\le t+1/n\), and
\(\log(t/|Z|)\le1/(5n)\). Multiplying the last bound by
\(k(2n-k)\le2nk\) costs at most \(2k/5\), while \(k/n\le k/2\).

Fix any constant \(L>0\), and put \(Z=\omega+w\),
\(|w|\le L/K_n\). This disk is eventually contained in the radius
\(1/(2n)\) disk. The physical integral and the Cauchy derivative bound give
\[
|C_{n,0}(\omega+w)-1|\le6n|w|e^{n|w|},
\quad
|C_{n,1}(\omega+w)-T_n|\le12ec\,nK_n|w|.
\tag{12}
\]
With \(\eta=1-e^w\), one has \(|\eta|\le2|w|\) and
\(|\eta+w|\le|w|^2\) eventually. Equations (10)--(11) bound the whole
higher-order tail by a geometric series:
\[
\sum_{k=2}^n|C_{n,k}(\omega+w)\eta^k|
\le\frac{6ec^2K_nK_{n-1}|\eta|^2}
 {1-cK_{n-1}|\eta|}
=O_L(n^{-2}).
\tag{13}
\]
The first two errors are \(O_L(n/K_n)=o(n^{-2})\). Hence the normal form is
\[
\boxed{\frac{D_n(\omega+w,1-e^w)}{\det H_n(\omega)}
=1-wT_n+O_L(n^{-2})}
\tag{14}
\]
uniformly on the entire scaled disk.

## 3. Rouché, symmetry and Newton displacement

Take \(L=48\). On the boundary, \(|1-wT_n|\ge3\), while the error in
(14) tends to zero. Rouché gives exactly one zero, counted with multiplicity;
it is therefore simple. At this zero,
\[
\zeta_n-\omega=T_n^{-1}(1+O(n^{-2})).
\tag{15}
\]
In particular it lies within \(24/K_n\) eventually, so all other zeros,
which lie outside \(48/K_n\), are farther away.

Real coefficients and the exact parity
\(P_n(-Z)=(-1)^nP_n(Z)\) make the root set invariant under
\(Z\mapsto-\overline Z\). This map fixes \(\omega\) and preserves the
unique-root disk. The unique root is therefore purely imaginary.

The repository's chain rule gives
\[
P_n'(\omega)/P_n(\omega)=F_n'(\omega)/F_n(\omega)-T_n,
\qquad F_n'(\omega)/F_n(\omega)=O(n^2).
\]
Together with (15), this proves (4). From (2),
\[
\zeta_n-\omega
=-i(-1)^n\frac{\chi_n(t)}{K_n\chi_{n-1}(t)}
 (1+O(n^{-1/2})),
\tag{16}
\]
which proves the side of approach. Applying the existing Stirling expansion
of \(K_n\) proves (5).

At fixed \(Z=\omega\), the same argument gives exactly one simple endpoint
zero in \(|Y-1|\le48/K_n\), with
\[
Y_n-1=-T_n^{-1}(1+O(n^{-2})).
\]
The reciprocal-conjugate image \(1/\overline{Y_n}\) remains in this disk,
because \(|Y_n-1|\le24/K_n\) eventually. Uniqueness therefore implies
\(Y_n=1/\overline{Y_n}\), so \(|Y_n|=1\).

## Verification and scope

The companion exact checks pass 1377 coefficient-majorant recurrences,
81 final endpoint identities, 1377 geometric bounds, 49 Laguerre factorial
determinants and 27 Gaussian-rational mixed coefficient bounds. The latter
use a synthetic physical measure \(du\) and rational imaginary endpoint
\(Z=6i\); they audit the algebra and majorant, not the actual period value.
Optional 180-digit root examples for \(n=4,6,8,10,12\) are illustrations,
not certified isolations or proofs of the asymptotic assertions.

The theorem determines location, simplicity, uniqueness and a sharper
distance scale for this family. It neither bounds the algebraic degrees or
heights sufficiently for a rational approximation theorem, nor constructs
reduced fractions \(p/q\). In particular it gives no conclusion about
bad approximability of \(\pi\), and claims no result for the original
high-dimensional parameter range.
