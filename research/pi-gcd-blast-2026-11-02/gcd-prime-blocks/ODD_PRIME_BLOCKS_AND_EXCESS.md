# Odd-prime blocks, exact gcd excess, and a divisibility counterexample

2026-11-02. This note studies the actual gcd in the same pi Padé family. It gives an all-n upper bound for each odd-prime excess, an exact generic valuation formula, and a counterexample to the stronger conjecture that no odd-prime excess occurs. These results do not bound the full reduced denominator sharply enough to decide the sparse-subsequence criterion.

## 1. The integer sequences and the p-adic question

Use ell_0=ell_1=1, w_0=0,w_1=1 and

\[
(n+1)x_{n+1}=(2n+1)x_n+n x_{n-1}\qquad(n\ge1)
\]

for ell and w. Set

\[
A_n=n!\ell_n,\quad B_n=n!w_n,\quad
g_n=\gcd(A_n,4B_n).
\]

Both integer sequences satisfy

\[
X_{n+1}=(2n+1)X_n+n^2X_{n-1}.
\tag{1}
\]

The reduced pi approximant is 4B_n/A_n=p_n/q_n with q_n=A_n/g_n. For an odd prime p, define the rational p-adic valuation v_p in the usual way (numerator exponent minus denominator exponent), and

\[
e_p(n)=v_p(g_n)-v_p(n!)
=\min(v_p(\ell_n),v_p(w_n)).
\tag{2}
\]

For n>=1 both ell_n and w_n are nonzero. The factor 4 in g_n is a p-adic unit for odd p. A positive e_p(n) is an extra odd-prime cancellation beyond the factorial scale.

## 2. Exact denominators and a uniform excess bound

The explicit Legendre coefficient formula shows ell_n belongs to Z[1/2]. Indeed its coefficients have the form

\[
2^{-n}\binom nk\binom{2n-2k}{n}
\]

up to signs already absorbed when evaluating at i. Hence ell_n is p-integral for every odd p.

The generating functions give a useful exact formula for w. Put

\[
F(t)=\sum_{n\ge0}\ell_n t^n=(1-2t-t^2)^{-1/2},
\qquad G(t)=\sum_{n\ge0}w_n t^n.
\]

Equation (1) before factorial scaling gives

\[
(1-2t-t^2)G'=(1+t)G+1,
\qquad (G/F)'=F.
\]

With G(0)=0 this proves

\[
\boxed{G(t)=F(t)\int_0^tF(s)\,ds,
\qquad
w_n=\sum_{j=1}^{n}\frac{\ell_{j-1}\ell_{n-j}}j.}
\tag{3}
\]

In particular v_p(w_n)>=−floor(log_p n). Here floor(log_p n) means the largest h with p^h<=n; no real logarithm calculation is needed.

The exact rational Wronskian is

\[
\ell_nw_{n-1}-\ell_{n-1}w_n=\frac{(-1)^n}{n}.
\tag{4}
\]

For n>=2 write h=floor(log_p(n−1)). Both terms on the left of (4) have valuation at least e_p(n)−h: use p-integrality of ell_(n−1), the lower bound v_p(w_(n−1))>=−h, and definition (2). The valuation of their difference is at least the minimum of the term valuations. Therefore

\[
\boxed{e_p(n)\le\lfloor\log_p(n-1)\rfloor-v_p(n)
\qquad(n\ge2,\ p\text{ odd}).}
\tag{5}
\]

This is an unconditional all-n bound for the actual gcd. It implies:

- If p>sqrt(n−1), any positive excess is at most one.
- If n=p^k, then e_p(n)<=−1: at least one p factor is missing relative to n!.
- Let M_(n−1)=lcm(1,...,n−1), and let odd(N) denote the odd part of a positive integer N. Then

\[
\boxed{\operatorname{odd}(n g_n)
\ \mid\ \operatorname{odd}(n! M_{n-1}).}
\tag{6}
\]

The last divisibility is a compact global expression of (5). Its lcm factor is still too costly to establish a bounded value of g_n sqrt(n)/n!, so it does not close the desired denominator estimate.

The forward Wronskian gives a second constraint. Apply (4) at n+1 and use v_p(w_(n+1))>=−floor(log_p(n+1)), while ell_(n+1) is p-integral. The same valuation argument proves

\[
\boxed{e_p(n)\le\min\left\{
\lfloor\log_p(n-1)\rfloor-v_p(n),
\lfloor\log_p(n+1)\rfloor-v_p(n+1)\right\}
\quad(n\ge2).}
\tag{6a}
\]

For example, if n+1=bp^k with 1<=b<p, the second expression is zero and there is no positive excess. This statement uses a single nonzero base-p digit for n+1, not arbitrary divisibility of n+1 by p.

There is a stronger necessary condition, useful for a sieve using ell_n alone. Suppose e_p(n)>0 and put u_p=v_p(ell_n)>0. In the backward Wronskian, ell_(n−1)w_n has strictly positive valuation, whereas the right side has valuation −v_p(n)<=0. Thus the other term must have exactly that latter valuation. Applying the same argument forward gives

\[
\boxed{v_p(w_{n-1})=-v_p(n)-u_p,
\qquad v_p(w_{n+1})=-v_p(n+1)-u_p.}
\tag{6b}
\]

The denominator bounds for these neighboring w values imply that **u_p itself**, not just the smaller quantity e_p(n), is at most the right side of (6a). Therefore a prime for which v_p(ell_n) exceeds that cap contributes no positive excess at all. Merely truncating its valuation to the cap would retain a prime which can instead be removed entirely. This refinement was identified during independent review and follows directly from the two Wronskians.

## 3. Lucas digits and the leading p-pole of w_n

The digit factorization here is a classical Schur/Lucas-type congruence, not claimed as a new result. A primary source discussing Legendre and scaled-Legendre versions is Tony D. Noe, *On the Divisibility of Generalized Central Trinomial Coefficients*, Journal of Integer Sequences 9 (2006), Article 06.2.7, Sections 6–7: https://cs.uwaterloo.ca/journals/JIS/VOL9/Noe/noe35.pdf. The following short derivation records the specialization needed for this family.

Work in F_p[[t]] with p odd and D(t)=1−2t−t². Frobenius and the square-root identity give

\[
F(t)=D(t)^{(p-1)/2}F(t^p)\pmod p.
\]

One formal justification is to square both sides and use D(t)^p=D(t^p); the two resulting square roots have the same constant term 1, so they agree in odd characteristic. The polynomial D^((p−1)/2) has degree p−1. Its coefficients of degree less than p agree with ell_0,...,ell_(p−1). Consequently

\[
\boxed{\ell_{bp+r}\equiv\ell_b\ell_r\pmod p,
\quad 0\le r<p,}
\tag{7}
\]

and iterating multiplies the ell values at every base-p digit. In particular

\[
\ell_{p-1}\equiv\chi_p:=(-1)^{(p-1)/2}\pmod p.
\]

Now fix n>=1 and write n=bp^h+r, where h=floor(log_p n), 1<=b<p and 0<=r<p^h. The rational number p^h w_n is p-integral by (3). In that sum, every term with p^h not dividing j vanishes modulo p after multiplication by p^h. Put j=kp^h in the remaining terms, and apply (7):

\[
\ell_{kp^h-1}\equiv\chi_p^h\ell_{k-1},
\qquad
\ell_{n-kp^h}\equiv\ell_r\ell_{b-k}.
\]

Thus

\[
\boxed{p^h w_n\equiv
\chi_p^h\ell_r w_b\pmod p.}
\tag{8}
\]

Since b<p, w_b is p-integral and its residue is well-defined. If ell_r w_b is nonzero modulo p, equation (8) gives the exact valuation, rather than just a bound:

\[
\boxed{v_p(w_n)=-h,
\qquad e_p(n)=-h.}
\tag{9}
\]

The second equality uses p-integrality of ell_n. Thus the generic case has a deficit of floor(log_p n) p factors. Any positive excess requires cancellation of this leading pole and further p-adic lifting. A necessary condition for positive excess is ell_n=0 modulo p, so at least one base-p digit must lie in the finite zero set {0<=r<p: ell_r=0 modulo p}. This condition alone is not sufficient.

Digit zeros give only a modulo-p criterion, not additive valuations: for n=72,p=17 the base-17 digits are (4,4), with ell_4=17/2, but v_17(ell_72)=1, not two.

Repeating that zero digit also does not sustain excess. Exact recurrence evaluation at n=4,72,1228,20880 (base-17 strings consisting entirely of 4) gives respectively (v_17 ell_n,v_17 w_n,e_17(n)) equal to (1,0,0), (1,1,1), (2,−1,−1), (2,0,0). This rejects that particular extrapolation; it is not an infinite-family theorem.

## 4. Exact prime-power transfer blocks

For either integer solution define v_n=(X_n,X_(n−1))^T, n>=1. Then

\[
v_{n+1}=T_n v_n,\qquad
T_n=\begin{pmatrix}2n+1&n^2\\1&0\end{pmatrix}.
\]

For every modulus m, T_(n+m)=T_n modulo m. Define S_m=T_m T_(m−1)...T_1 and the prefix C_s=T_s...T_1, C_0=I. Writing n−1=km+s with 0<=s<m gives the exact finite computation

\[
\boxed{v_n\equiv C_s S_m^k v_1\pmod m,}
\tag{10}
\]

with v_1=(1,1)^T for A and v_1=(1,0)^T for B. This holds in particular for every prime power m=p^j. Matrix exponentiation computes the block power without iterating through the full index n. It does not assert that every block is invertible.

For an odd prime, the full prime block has a particularly simple nilpotent reduction:

\[
\boxed{S_p\equiv-\chi_p
\begin{pmatrix}1&-1\\1&-1\end{pmatrix}\pmod p,
\qquad S_p^2\equiv0\pmod p.}
\tag{11}
\]

To prove this, S_p sends the A initial vector to (A_(p+1),A_p), zero modulo p, because ell is p-integral. Hence its column sum is zero. Formula (8) at n=p and p+1 gives pw_p=pw_(p+1)=chi_p modulo p. Wilson's identity (p−1)!=−1 modulo p then gives B_p=B_(p+1)=−chi_p modulo p, fixing its first column. The displayed matrix follows. Wilson's identity itself follows by pairing every nonzero residue with its inverse; only 1 and −1 remain unpaired.

Higher powers must use (10) modulo p^j; the modulo-p nilpotent block alone does not determine higher valuations or a large gcd.

## 5. Counterexample to no odd excess

The conjectural bound e_p(n)<=0 for every odd p is false. A complete exact check through n=500 has its first counterexample at

\[
\boxed{n=72,\quad p=17,\quad
v_{17}(A_{72})=v_{17}(B_{72})=5,
\quad v_{17}(72!)=4.}
\]

Therefore v_17(g_72)=5 and e_17(72)=1. Equivalently the odd part of g_n need not divide the odd part of n!. Other exact examples include (n,p)=(107,73),(166,131),(170,41),(176,23). For n=480,p=31, the A and B valuations are respectively 16 and 17, while v_31(480!)=15, again giving excess one.

The verifier records the first nonzero residues after removal of these p powers and separately checks them via modular recurrence. It also verifies the Lucas relation, the leading-pole formula, prime and prime-power transfer blocks, and the uniform bound (5). These finite checks certify the counterexamples and audit the formulas; the general statements are proved above.

The remaining issue is quantitative simultaneous control of exceptional primes and their lifts. Neither a modulo-p digit criterion nor a finite absence of composite excess can be promoted into an all-n upper bound for their product.
