# Exact 2-adic valuation of the Legendre numerator sequence

2026-11-02. Let a₀=a₁=1 and, for integers n≥1,

\[
a_{n+1}=(2n+1)a_n+n^2a_{n-1}.
\]

For a nonzero rational number x, v₂(x) denotes the exponent of 2 in its numerator minus that in its denominator. The following formulas hold for every integer m≥0:

\[
\boxed{\begin{aligned}
v_2(a_{4m})&=2m,\\
v_2(a_{4m+1})&=2m,\\
v_2(a_{4m+2})&=2m+2,\\
v_2(a_{4m+3})&=2m+3+v_2(m+1).
\end{aligned}}
\tag{1}
\]

These are exact all-index statements, not extrapolations from finite values. The argument uses a rational-coefficient duplication identity and a finite polynomial congruence. No assertion about odd-prime gcd growth is made here.

## 1. A binomial expression and its even-index valuation

Put ℓ_n=a_n/n!. Dividing the recurrence by (n+1)! gives

\[
(n+1)\ell_{n+1}=(2n+1)\ell_n+n\ell_{n-1}.
\]

Its generating function is (1−2z−z²)^(−1/2). Equivalently, coefficient expansion gives

\[
\ell_n=\sum_{k=0}^{\lfloor n/2\rfloor}
\frac{n!}{2^k(k!)^2(n-2k)!}.
\tag{2}
\]

The identity can also be checked directly against the displayed recurrence and initial values. In particular all denominators of ℓ_n are powers of 2.

For m≥0, reverse the sum in (2) to write

\[
a_{2m}=C_mS_m,\qquad
C_m=\frac{[(2m)!]^2}{2^m(m!)^2},\qquad
S_m=\sum_{j=0}^m R_j(m),
\]
\[
R_j(m)=\frac{2^j(m)_j^2}{(2j)!},\qquad
(m)_j=m(m-1)\cdots(m-j+1).
\tag{3}
\]

Here (m)₀=1. Since v₂((2m)!)=m+v₂(m!), one has v₂(C_m)=m. Moreover

\[
R_j(m)=\binom mj^2\frac{2^j}{\binom{2j}j},
\quad
v_2(R_j(m))\ge j-s_2(j)=v_2(j!),
\tag{4}
\]

where s₂(j) is the number of ones in the binary expansion of j. The central-binomial valuation used here follows immediately from v₂(N!)=N−s₂(N).

The j=0,1 terms sum to 1+m². If m is even, all j≥2 terms are even, so S_m is odd. If m is odd, terms j≥4 are divisible by 8, while

\[
R_2(m)+R_3(m)
=\frac{[m(m-1)]^2[(m-2)^2+15]}{90}
\]

is divisible by 16: the squared factor contributes at least two powers of 2, the bracket at least three, and the denominator has valuation one. This also holds when m=1, when the falling products vanish. Therefore S_m≡1+m²≡2 modulo 4 for odd m. We have proved

\[
v_2(a_{2m})=m+(m\bmod2).
\tag{5}
\]

## 2. The odd-index sum and a zero of a convergent rational series

For odd indices, the same reversal yields

\[
a_{2m+1}=D_m T(m),\qquad
D_m=\frac{[(2m+1)!]^2}{2^m(m!)^2},\qquad v_2(D_m)=m,
\]
\[
T(m)=\sum_{j=0}^m t_j\binom mj^2,
\qquad t_j=\frac{2^j(j!)^2}{(2j+1)!}.
\tag{6}
\]

One has v₂(t_j)=j−s₂(j)=v₂(j!), which tends to infinity. Thus the series ∑t_j is convergent in the 2-adic sense: for every modulus 2^N, its sufficiently long partial sums have a fixed residue. We claim that this limit is zero.

To prove it without logarithms, let

\[
F(z)=\sum_{j\ge0}c_jz^j,\qquad c_j=\frac{(j!)^2}{(2j+1)!}.
\]

The coefficient ratio proves the formal differential equation

\[
z(4-z)F'(z)+(2-z)F(z)=2.
\tag{7}
\]

Set H(z)=(2−z)F(z(4−z))−2F(z). Direct differentiation using (7), and the identity 4−z(4−z)=(2−z)², gives

\[
z(4-z)H'(z)+(2-z)H(z)=0,\qquad H(0)=0.
\]

Its coefficient recurrence forces every coefficient to vanish. Hence the formal duplication identity is

\[
\boxed{(2-z)F(z(4-z))=2F(z).}
\tag{8}
\]

It is legitimate to substitute z=2 in the 2-adic sense. Indeed v₂(c_j)=−s₂(j), so both F(2) and F(4) converge. In expanding the composition at z=2, every term of c_jz^j(4−z)^j has valuation at least 2j−s₂(j), which tends to infinity; hence regrouping is valid. Substitution in (8) gives 0=2F(2), and therefore

\[
\boxed{\sum_{j\ge0}t_j=F(2)=0.}
\tag{9}
\]

This statement concerns the limit of rational residues modulo increasing powers of 2. The same series converges over the real numbers to a different, nonzero real value. No identification of those two limits is used.

## 3. The odd-index valuation reduces to five terms

When m is even, all terms j≥2 in (6) are even and the j=0,1 terms give T(m)≡1+m²/3≡1 modulo 2. Thus v₂(T(m))=0.

Now let m be odd, put u=m+1, and set s=v₂(u)≥1. Define

\[
P_j(u)=\prod_{r=1}^j(1-u/r)
=(-1)^j\binom{u-1}{j}.
\]

For the positive integer u in question, P_j(u) is an integer. Vandermonde's identity gives

\[
P_j(u)-1=\sum_{r=1}^j(-1)^r\binom ur.
\]

Since binom(u,r)=(u/r)binom(u−1,r−1),

\[
v_2(P_j(u)-1)\ge s-\lfloor\log_2j\rfloor.
\]

The other factor P_j(u)+1 is an integer, so

\[
v_2\left(t_j\frac{P_j(u)^2-1}{u}\right)
\ge v_2(j!)-\lfloor\log_2j\rfloor.
\tag{10}
\]

The right side is at least 2 for every j≥6. For j=6,7 this is direct; for 2^r≤j<2^(r+1) with r≥3, it is at least 2^r−1−r≥4.

Equation (9) and the finite vanishing of binom(m,j) for j>m therefore imply

\[
\frac{T(m)}{u}
=\sum_{j\ge1}t_j\frac{P_j(u)^2-1}{u}
\equiv\sum_{j=1}^5t_j\frac{P_j(u)^2-1}{u}\pmod4.
\tag{11}
\]

For a finite certificate of this last sum, let it be H₅(u). Exact polynomial expansion gives

\[
H_5(2v)-2=\frac4{155925}W(v),
\tag{12}
\]

where the integer polynomial W has coefficients, in ascending order,

\[
(-133890,\ 106091,\ -112150,\ 95560,\ -57235,
\ 25853,\ -8200,\ 1800,\ -240,\ 16).
\]

Its denominator is odd. Since u is even, (11)–(12) give T(m)/u≡2 modulo 4, so

\[
\boxed{v_2(T(m))=1+v_2(m+1)\quad(m\text{ odd}).}
\tag{13}
\]

Combining (6), (13), and the even-m case proves

\[
v_2(a_{2m+1})=
\begin{cases}
m,&m\text{ even},\\
m+1+v_2(m+1),&m\text{ odd}.
\end{cases}
\]

Together with (5), this is precisely (1).

## 4. Consequence for the gcd attack

For g_n=gcd(a_n,4b_n), whatever the second sequence b_n is, one always has v₂(g_n)≤v₂(a_n). Formula (1) therefore supplies an exact upper bound on the 2-part of g_n. Since v₂(n!)=n−s₂(n), it implies

\[
v_2(g_n)-v_2(n!)\le-\frac n2+O(\log n).
\]

Thus the 2-part alone loses an exponential factor relative to n!. To obtain g_n√n/n!→∞, the odd-prime part would have to compensate for this loss. An all-prime compensation bound is a separate problem; the displayed 2-adic result does not establish it.

## Verification

The standard-library verifier independently derives the finite polynomial certificate (12), verifies the coefficient differential and duplication identities to finite order, checks the elementary valuation inequalities, and compares (1) with the integer recurrence through n=2000. The general theorem follows from the proof above. Finite recurrence checks alone would not prove it.
