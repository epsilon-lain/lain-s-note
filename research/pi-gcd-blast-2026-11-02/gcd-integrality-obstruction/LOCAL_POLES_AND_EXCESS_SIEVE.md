# Exact neighboring poles and a stronger arithmetic sieve

2026-11-02. This is an independent arithmetic attack on the fixed Legendre--arctangent family. It gives a stronger necessary condition for exceptional odd-prime cancellation and a factorization-free refinement of the reduced-denominator sieve. It does not prove a global bound strong enough to exclude all quality-zero subsequences.

## 1. Setup and the arithmetic target

Use the rational polynomials ell_n(a), w_n(a) with initial values
ell_0=1, ell_1=a, w_0=0, w_1=1 and

\[
(n+1)x_{n+1}=(2n+1)a x_n+n x_{n-1}.
\]

Unless an argument is displayed, evaluation below is at a=1. Put

\[
A_n=n!\ell_n,\quad B_n=n!w_n,\quad
g_n=\gcd(A_n,4B_n),\quad
\frac{p_n}{q_n}=\frac{4B_n/g_n}{A_n/g_n}.
\]

The preceding exact integral and Legendre error estimates give

\[
Q_n=q_n^2|\pi-p_n/q_n|
>\frac{2}{2n+1}\left(\frac{n!}{g_n}\right)^2,
\tag{1}
\]

and the sharper asymptotic

\[
Q_n\sim\frac{\sqrt2}{n}\left(\frac{n!}{g_n}\right)^2.
\tag{2}
\]

The open issue for this family is whether g_n can exceed n!/sqrt(n) by an unbounded factor along a sparse subsequence.

For every odd prime p, ell_n(a) has p-integral coefficients. The convolution identity

\[
w_n(a)=\sum_{j=1}^{n}\frac{\ell_{j-1}(a)\ell_{n-j}(a)}j
\tag{3}
\]

therefore shows that lcm(1,...,n)*w_n(a) has p-integral coefficients. In particular, at a=1,

\[
v_p(w_n)\ge-\lfloor\log_p n\rfloor.
\tag{4}
\]

These are coefficient-denominator statements, not guesses about cancellation after evaluation.

## 2. A differential identity independently restricts common divisors

Differentiating the polynomial recurrences, or proving the following by induction, gives

\[
(1+a^2)\ell_n'=n(a\ell_n+\ell_{n-1}),
\]
\[
(1+a^2)w_n'=n(a w_n+w_{n-1})-\ell_n.
\]

The polynomial Wronskian is
ell_n*w_(n-1)-ell_(n-1)*w_n=(-1)^n/n. Subtracting the two derivative identities after multiplying by ell_n and w_n proves

\[
\boxed{(1+a^2)(\ell_n w_n'-\ell_n'w_n)
=(-1)^n-\ell_n^2.}
\tag{5}
\]

At a=1, suppose v_p(ell_n)=alpha>0 and v_p(w_n)>0. The right-hand side of (5) is a p-adic unit. Since ell_n' is p-integral and w_n' has valuation at least -floor(log_p n), alpha cannot exceed floor(log_p n). This independently explains why positive common divisors are restricted to primes and powers already present in the coefficient denominators of w_n. It is not a bound saying that no odd excess exists.

The neighboring Wronskians sharpen (5) substantially.

## 3. Positive excess forces exact neighboring poles

For n>=2 and p odd write

\[
\alpha=v_p(\ell_n),\quad\beta=v_p(w_n),\quad
e_p(n)=v_p(g_n)-v_p(n!)=\min(\alpha,\beta).
\]

Suppose e_p(n)>0. In the backward Wronskian

\[
\ell_n w_{n-1}-\ell_{n-1}w_n=(-1)^n/n,
\]

the second term has positive valuation because ell_(n-1) is p-integral and beta>0. The right-hand side has valuation -v_p(n)<=0. Consequently the first term must have exactly that nonpositive valuation; cancellation with the positive-valuation term is impossible. The forward identity gives the analogous conclusion. Thus

\[
\boxed{v_p(w_{n-1})=-v_p(n)-\alpha,
\qquad
v_p(w_{n+1})=-v_p(n+1)-\alpha.}
\tag{6}
\]

Combining (6) with (4) proves the stronger necessary condition

\[
\boxed{e_p(n)>0\quad\Longrightarrow\quad
0<\alpha\le c_p(n),}
\tag{7}
\]
\[
c_p(n)=\min\{\lfloor\log_p(n-1)\rfloor-v_p(n),
\lfloor\log_p(n+1)\rfloor-v_p(n+1)\}.
\]

This bounds the **entire** p-adic valuation of ell_n when excess occurs, rather than only the smaller value e_p(n). In particular, if v_p(ell_n)>c_p(n), p contributes no positive excess at all. Merely retaining p^c_p(n) in an upper bound would lose this information.

The companion [odd-prime note](../gcd-prime-blocks/ODD_PRIME_BLOCKS_AND_EXCESS.md) proves the unconditional bound e_p(n)<=c_p(n). Its proof and the additional Lucas and transfer-block identities were independently checked and found sound. Equation (7) is a useful strict-case refinement of that argument.

## 4. A sieve using only the first integer solution

Define

\[
U_n=2^n\ell_n\in\mathbb Z_{>0},\qquad
s_n=v_2(U_n),\qquad O_n=U_n/2^{s_n}.
\]

These integers can be generated without factorials or the second solution:

\[
U_0=1,\quad U_1=2,\qquad
(n+1)U_{n+1}=2(2n+1)U_n+4nU_{n-1}.
\tag{8}
\]

Let L_j=lcm(1,...,j), and for n>=2 set

\[
E_n=\gcd\left(
\frac{L_{n-1}}{\gcd(L_{n-1},n)},
\frac{L_{n+1}}{n+1}\right),\qquad
G_n=\gcd(O_n,\operatorname{odd}(E_n)).
\tag{9}
\]

Here odd denotes removal of all factors of two. The odd-prime exponents in E_n are max(0,c_p(n)). Consequently the numerator of the odd part of the reduced rational number g_n/n! divides G_n. This is the two-neighbor version of the root investigation's basic sieve.

Equation (7) improves the divisor G_n. Form

\[
H_0=\gcd(O_n/G_n,G_n),\qquad
H_{j+1}=\gcd(G_n,H_j^2).
\tag{10}
\]

This iteration stabilizes at H_infinity. For a prime p with alpha=v_p(O_n)>0 and c=max(0,c_p(n)), its exponent in G_n is min(alpha,c). Its exponent in H_0 is positive exactly when alpha>c>0. Repeated squaring and gcd therefore fill out precisely the full prime-power part of G_n at such primes. Define

\[
\widehat G_n=G_n/H_\infty.
\tag{11}
\]

Every prime removed in (11) is forbidden from contributing positive excess by (7). Thus the odd numerator of g_n/n! divides the smaller integer Ghat_n as well. This construction requires gcds and integer arithmetic but no prime factorization.

Since g_n divides A_n, its power of two satisfies

\[
v_2(g_n/n!)\le v_2(\ell_n)=s_n-n.
\]

Discarding negative odd-prime valuations only increases a positive rational number. Therefore

\[
\boxed{\frac{g_n}{n!}\le2^{s_n-n}\widehat G_n,}
\tag{12}
\]

and the strict analytic floor (1) gives the fully explicit certificate

\[
\boxed{Q_n>
\frac{2\,2^{2(n-s_n)}}{(2n+1)\widehat G_n^2}.}
\tag{13}
\]

The weaker version with G_n or with only the backward factor in E_n is also valid. No formula for s_n, approximate value of pi, computation of B_n, or assumption about generic prime behavior is required for (13). When its right-hand side exceeds a proposed quality threshold, that index is excluded rigorously. A finite range of exclusions is not an all-n result.

## 5. What this attack does not settle

Odd excess really occurs: at n=72,p=17 both A_n and B_n have valuation five, whereas v_17(72!)=four. Thus deleting exceptional odd cancellation from the argument would be false.

The lcm constraints, exact neighboring poles, and sieve improve what can be checked. They do not yet prove that the product of all surviving exceptional odd-prime factors stays below the exponential scale needed in (2). In particular, sparsely occurring large simultaneous gains are not excluded by a local pole argument.

The companion standard-library verifier checks the polynomial differential identities exactly, compares the factorization-free sieve against the actual reduced gcd through a finite range, and records strict exclusions obtained from (13). The all-n validity of the sieve is proved above; no extrapolation from its finite output is used.
