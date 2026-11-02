# A sharper first-coordinate filter from the actual Wronskian

2026-11-02. This note concerns the fixed Legendre arctangent family at
repository checkpoint `05e6e140c6a41a5ed62f38d6d583487d0455bc7c`.
This note applies the already established exact Wronskian more sharply.
Neither that Wronskian nor its differential form below is claimed as a new
identity. The project-level advance is a stricter first-coordinate valuation
filter, which removes complete prime factors from the gain-only multiple.

Let
\[
A_0=A_1=1,\quad B_0=0,\ B_1=1,\qquad
X_{n+1}=(2n+1)X_n+n^2X_{n-1}.
\]
Write \(a_n=A_n/n!\), \(b_n=B_n/n!\),
\(g_n=\gcd(A_n,4B_n)\), and
\(s_p(n)=v_p(g_n)-v_p(n!)\) for odd primes.

## 1. A constant witness using the actual Cauchy numerator

Let \(P_n(z)\) be the standard Legendre polynomial, and let
\[
W_n(z)=\frac12\int_{-1}^1
\frac{P_n(z)-P_n(t)}{z-t}\,dt.
\]
These polynomials satisfy the same recurrence
\[
(n+1)Y_{n+1}=(2n+1)zY_n-nY_{n-1},
\]
with \((P_0,P_1)=(1,z)\) and \((W_0,W_1)=(0,1)\).
The fixed integer sequences are exactly
\[
P_n(i)=i^n a_n,\qquad W_n(i)=i^{n-1}b_n.
\]

For all \(n\ge0\), one has the polynomial identity
\[
\boxed{P_n(z)^2+(1-z^2)
\bigl(P_n'(z)W_n(z)-P_n(z)W_n'(z)\bigr)=1.}
\tag{1}
\]
Here primes denote derivatives with respect to \(z\), not the index.

One elementary proof uses the identities, for \(n\ge1\),
\[
(1-z^2)P_n'=n(P_{n-1}-zP_n),
\]
\[
(1-z^2)W_n'=n(W_{n-1}-zW_n)+P_n.
\tag{2}
\]
The first is the Legendre derivative relation. The second follows by
induction from the displayed polynomial recurrence and the initial
values \(W_0=0,W_1=1\); equivalently, differentiating the Cauchy
numerator gives the same relation. The recurrence also yields
\[
P_{n-1}W_n-P_nW_{n-1}=1/n:
\]
its value is one at \(n=1\), and each succeeding determinant is
multiplied by \(n/(n+1)\). Substitution into (2) proves (1).
For \(n=0\), (1) is immediate.

The constant on the right depends on the actual normalization of the
second solution. This is not an identity inferred solely from setting a
free central recurrence variable to zero, so the generic window theorem
does not forbid it.

## 2. Stronger direct discrete proof, including primes dividing the index

The differential identity in Section 1, evaluated at $z=i$ and simplified
using the derivative relations, is exactly the existing discrete Wronskian
\[
n(a_nb_{n-1}-a_{n-1}b_n)=(-1)^n.                 \tag{3}
\]
Indeed the integer recurrence proves
$A_nB_{n-1}-A_{n-1}B_n=(-1)^n((n-1)!)^2$.
Thus the differential form does not produce an independent second witness.

Fix $n\ge2$ and an odd prime $p$. Put
\[
h_p=v_p(n),\qquad e_p=\lfloor\log_p(n-1)\rfloor,
\quad k_p=e_p-h_p,\quad m_p=v_p(a_n),
\]
where $e_p=0$ for $p>n-1$ and $k_p$ may be negative.
The established integrality and convolution are
\[
a_j\in\mathbb Z[1/2],\qquad
b_j=\sum_{r=1}^{j}\frac{a_{r-1}a_{j-r}}r.
\]
Consequently $v_p(a_{n-1})\ge0$ and $v_p(b_{n-1})\ge-e_p$.
If $m_p>k_p$, then
\[
v_p(na_nb_{n-1})\ge h_p+m_p-e_p>0.
\]
Equation (3) forces $na_{n-1}b_n$ to have valuation zero. Hence
\[
\boxed{m_p>k_p\ \Longrightarrow\
v_p(b_n)=-v_p(n)-v_p(a_{n-1})\le0.}              \tag{4}
\]
There is no field extension and no exception for $p\mid n$. In particular,
\[
\boxed{s_p(n)>0\ \Longrightarrow\
1\le v_p(a_n)\le
\lfloor\log_p(n-1)\rfloor-v_p(n).}              \tag{5}
\]
This strengthens the earlier cap on the minimum
$\min(v_p(a_n),v_p(b_n))$ to a cap on the first coordinate alone.
For $n=p^r$, one has $k_p=-1$, so (4) applies even when $m_p=0$;
prime-power indices are part of the exact audit.

For $p>\sqrt n$, (5) implies that positive surplus requires
$v_p(C_n)=1$, with $C_n=\operatorname{odd}(2^na_n)$.
Thus a selected prime with $p^2\mid C_n$ cannot provide positive surplus.

## 3. Refined gain-only and deficit certificates

Define
\[
E_n=\prod_{p\text{ odd}}p^{\max(s_p(n),0)},\qquad
D_n=\prod_{p\text{ odd}}p^{\max(-s_p(n),0)}.
\]
The old gain-only multiple was
\[
G_n^{\rm old}=\prod_{p<n\text{ odd}}
p^{\min(m_p,\max(k_p,0))}.
\]
By (5), every prime with $m_p>k_p$ must be removed completely:
\[
\boxed{E_n\mid G_n^{\rm first},\qquad
G_n^{\rm first}=
\prod_{\substack{p<n\text{ odd}\\1\le m_p\le k_p}}p^{m_p}.} \tag{6}
\]
For these removed primes, (4) gives the exact nonpositive surplus
\[
\boxed{s_p(n)=-v_p(n)-v_p(a_{n-1}).}              \tag{7}
\]
Any strictly negative right side also provides a certified factor of
$D_n$ without constructing $B_n$.

If $M\mid D_n$ is a compatible existing deficit certificate, and
$U_n=2^na_n$, the established strict quality floor is strengthened to
\[
Q_n>
\frac{2\,2^{2(n-v_2(U_n))}M^2}
{(2n+1)(G_n^{\rm first})^2}.                    \tag{8}
\]
At $n=68$, the old multiple contains $37$, but
$v_{37}(C_{68})=2>k_{37}=1$. Equation (6) removes it, raising the
quality floor by $37^2$ at fixed $M$. This is a stricter finite
rejection certificate, not a new rational approximation.

## 4. Refined asymptotic budget

Use the preceding growing-band threshold \(t\asymp n^{2/3}\), and set
\[
W(n)=\sum_{t<p\le n,\ p\mid C_n}\ln p,
\qquad
W_1(n)=\sum_{t<p\le n,\ v_p(C_n)=1}\ln p.
\]
The selected positive-surplus product now divides the product counted
by \(W_1\). All primes below \(t\) contribute only \(o(n)\) to the
positive budget, by the established cap. Hence
\[
\ln E_n\le W_1(n)+o(n).
\]
Retaining the existing multi-band deficit estimate gives
\[
\boxed{J_n\le W(n)+W_1(n)-n+o(n),}
\]
\[
\boxed{\ln Q_n\ge(2+\ln2)n-2W(n)-2W_1(n)+o(n).}
\tag{9}
\]
Compared with the previous floor, the correction is twice the logarithmic
weight of selected prime factors whose multiplicity in \(C_n\) is at
least two. Nothing here proves that correction has positive linear
size, nor bounds the remaining simple-prime concentration.

There is nevertheless a further necessary concentration bound. The known
height formula is
\[
\ln C_n=Hn+O(\ln n),\qquad H=\ln(2+\sqrt2).
\]
Every prime counted in \(W-W_1\) occurs at least twice in \(C_n\), so
\[
2W(n)-W_1(n)\le\ln C_n.                         \tag{10}
\]
If \(Q_n\) stays bounded along a subsequence, (9) requires
\(W+W_1\ge Kn-o(n)\), with \(K=1+\tfrac12\ln2\).
Eliminating \(W\) between these two inequalities gives
\[
\boxed{W_1(n)\ge
\left(\frac{2+\ln2-\ln(2+\sqrt2)}3-o(1)\right)n
=(0.4884000010\ldots-o(1))n.}                   \tag{11}
\]
Thus bounded-quality candidates require at least
\((0.4884000010\ldots-o(1))n/\ln n\) selected distinct large primes
that occur **exactly once** in \(C_n\). These are eligible prime
factors, not a claim that every one supplies actual positive surplus.
Equation (11) is a necessary condition, not a contradiction with the
available height.

## Verification and boundary

`verify_differential_witness.py` checks (1) as an exact rational polynomial
identity for degrees 0 through 25. It independently builds the actual
integer recurrence through the stated finite limit, checks the integer
determinant, tests all odd primes up to and including the index, and
verifies (4), (5), and $E_n\mid G_n^{\rm first}$. The accompanying JSON
records counts and boundary cases. These tests supplement the general
proof above.

The new project-level input is the sharper application of the existing
Wronskian to the first coordinate. Its denominator cost still prevents an
immediate global contradiction;
the proof does not bound \(W_1(n)\) below the necessary threshold and
does not construct a quality-zero subsequence.
