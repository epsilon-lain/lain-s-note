# Primitive-divisor barrier for elliptic formal-group cancellation

Date: 2026-11-20  
Status: **proved asymptotic denominator lower bound for every fixed finite integer combination**, using elliptic divisibility sequences (EDS). This explains why moving cancellation weight to larger multiplication indices has a genuine arithmetic cost.

## 1. Setup

Take
\[
E:\ y^2=x^3-2,\qquad P=(3,5),\qquad \xi=-x/y.
\]
The point \(P\) is non-torsion: the discriminant is \(-1728\), while \(25=y(P)^2\) does not divide \(1728\), so Lutz–Nagell excludes torsion.

For \(n\ge1\), write
\[
nP=\left(\frac{A_n}{B_n^2},\frac{C_n}{B_n^3}\right),
\qquad B_n>0,\quad \gcd(B_n,A_nC_n)=1.
\]
The sequence \(B_n\) is the elliptic divisibility sequence associated to \((E,P)\). In particular, \(B_m\mid B_n\) if \(m\mid n\), and for good primes its primitive divisors record the exact order of the reduction of \(P\).

Let a finite set \(S\subset\{2,3,\ldots\}\) and nonzero integer weights \(\lambda_j\) be fixed. Let
\[
J=\max S,\qquad \lambda_J\ne0,
\]
and define
\[
\mathcal B_{\lambda,s}
=\sum_{j\in S}\lambda_j
\frac{\xi([js]P)}{\xi([s]P)}\in\mathbb Q.
\tag{1}
\]

## 2. Primitive primes of \(B_{2Js}\) force poles in the last summand

Let \(N=2Js\), and let \(p\) be a primitive prime divisor of \(B_N\) of good reduction, with \(p>3\). Then the reduction \(\bar P\in E(\mathbb F_p)\) has exact order \(N\).

Put \(\bar Q=[s]\bar P\). Its order is exactly \(2J\).

For every \(j<J\), the order of \([j]\bar Q\) is
\[
\frac{2J}{\gcd(2J,j)}>2.
\]
Hence \([j]\bar Q\) is not a nontrivial 2-torsion point. The rational function \(\xi=-x/y\) has poles precisely at nontrivial 2-torsion points, so
\[
v_p\!\left(\frac{\xi([js]P)}{\xi([s]P)}\right)\ge0
\qquad(j<J).
\tag{2}
\]

For \(j=J\), the point \([J]\bar Q\) has exact order 2. Thus \(\xi\) has a simple pole there, while \(\xi(\bar Q)\) is a unit. Since \(p>3\), multiplication by 2 is étale at this 2-torsion point. A local parameter at the 2-torsion point therefore has the same \(p\)-adic valuation as the formal parameter \(\xi([N]P)\) at the identity. From
\[
\xi([N]P)=-\frac{A_NB_N}{C_N}
\]
and \(\gcd(B_N,A_NC_N)=1\), that valuation is \(v_p(B_N)\). Consequently
\[
\boxed{
v_p\!\left(\frac{\xi([Js]P)}{\xi([s]P)}\right)
=-v_p(B_N)<0.
}
\tag{3}
\]

Equations (2)–(3) show that the \(J\)-term is the unique term of negative \(p\)-adic valuation in (1). If \(p\nmid\lambda_J\), no cancellation is possible:
\[
v_p(\mathcal B_{\lambda,s})=-v_p(B_N).
\tag{4}
\]

For fixed \(\lambda_J\), all its prime divisors form a finite set. A fixed good prime can be primitive for at most one index \(N\), because the order of \(\bar P\) modulo that prime is fixed. The bad primes are also a finite set. Silverman's primitive-divisor theorem therefore implies that, for all sufficiently large \(s\), every prime in the primitive part \(D_N\) of \(B_N\) is good, exceeds 3, and does not divide \(\lambda_J\). It follows from (4) that
\[
\boxed{D_{2Js}\mid \operatorname{den}(\mathcal B_{\lambda,s})}
\qquad\text{for all sufficiently large }s.
\tag{5}
\]

Here \(D_N\) is the product of the full prime powers in \(B_N\) whose prime divisors are primitive at index \(N\).

## 3. Quantitative height obstruction

For the \(\mathbb Z\)-indexed EDS of an elliptic curve with complex multiplication over a quadratic extension, Streng proves the sharper primitive-part asymptotic
\[
\log D_N
=
\widehat h(P)
\left(\prod_{\ell\mid N}(1-\ell^{-2})\right)N^2
+O_\epsilon(N^\epsilon)
\qquad(\epsilon>0).
\tag{6}
\]
Our curve has \(j(E)=0\), with CM field \(\mathbb Q(\sqrt{-3})\), so this result applies to its rational-point EDS. Since
\[
\prod_{\ell\mid N}(1-\ell^{-2})
\ge \prod_{\ell\ {\rm prime}}(1-\ell^{-2})
=\frac1{\zeta(2)}>0.6079,
\]
equations (5)–(6), with \(N=2Js\), give
\[
\boxed{
\liminf_{s\to\infty}
\frac{h(\mathcal B_{\lambda,s})}{s^2}
\ge
\frac{4J^2\widehat h(P)}{\zeta(2)}.
}
\tag{7}
\]
Here \(h(a/b)=\log\max(|a|,b)\) for a reduced rational \(a/b\), \(b>0\).

**Interpretation:** for every fixed finite integer combination, its reduced denominator is forced to contain a primitive EDS part at index \(2Js\). Its unavoidable quadratic height cost grows at least like \(J^2s^2\), where \(J\) is the largest multiplication index.

This is stronger than the heuristic proxy \(\sum_j|\lambda_j|j^2\): the lower bound comes from actual prime divisors that cannot cancel in the rational sum.

## 4. Consequence for order-\(m\) formal-group cancellation

To obtain
\[
\sum_j\lambda_j f(j)=f(1)
\qquad\text{for every polynomial }\deg f\le m,
\tag{8}
\]
the support must contain at least \(m+1\) distinct indices. If all indices were at most \(m+2\), the support would have to be exactly \(\{2,\ldots,m+2\}\), and the Vandermonde system has the unique solution
\[
\lambda_j=(-1)^j\binom{m+1}{j-1}.
\]
Therefore every genuinely different integer moment scheme has
\[
J\ge m+3.
\]
Combining with (7):

- the binomial scheme has the smallest possible \(J=m+2\);
- every alternative scheme has an asymptotic primitive-denominator lower-bound coefficient at least
  \[
  \left(\frac{m+3}{m+2}\right)^2
  \]
  times the binomial lower-bound coefficient.

This does **not** prove that binomial weights minimize the full rational height among every imaginable construction. It proves that, within finite integer combinations of multiplication maps satisfying the moment conditions, moving the largest multiplier upward necessarily raises the primitive-divisor height barrier.

## 5. A second EDS phenomenon: coprime multipliers have disjoint good denominator primes

For a good prime \(p>3\), if \(p\) divides the reduced denominator of
\[
R_j(s)=\frac{\xi([js]P)}{\xi([s]P)},
\]
then the reduction of \([s]P\) is neither the identity, nor 2-torsion, nor a 3-torsion zero of \(\xi\); moreover \([j][s]P\) is 2-torsion. If the same \(p\) divides the denominator of \(R_k(s)\), the order \(r\) of \([s]\bar P\) must satisfy
\[
r\mid 2j,\qquad r\mid 2k,\qquad r>2.
\]
Thus \(r\mid2\gcd(j,k)\), which is impossible if \(\gcd(j,k)=1\). Hence
\[
\boxed{
\gcd(j,k)=1
\implies
\gcd\bigl(\operatorname{den}R_j(s),\operatorname{den}R_k(s)\bigr)
\text{ has no prime factor }p>3.
}
\tag{9}
\]
For our curve, the only excluded primes are the bad primes 2 and 3.

This is the opposite of the hoped-for “shared denominator” shortcut: coprime multipliers tend to bring disjoint good prime supports. EDS divisibility can create overlap when multipliers share factors, but the moment condition \(\sum_j\lambda_jj=1\) forces the gcd of all support indices to be 1.

## 6. Exact finite experiment

For \(E:y^2=x^3-2\), \(P=(3,5)\), the signed division-polynomial sequence starts
\[
W_1=1,\quad W_2=10,\quad W_3=171,\quad W_4=-7660,
\]
and satisfies
\[
W_{m+2}W_{m-2}
=W_{m+1}W_{m-1}W_2^2-W_3W_1W_m^2.
\]
The division-polynomial identities give
\[
\phi_n=3W_n^2-W_{n+1}W_{n-1},\qquad
\xi([n]P)=-\frac{2\phi_nW_n^2}{W_{2n}}.
\tag{10}
\]
For \(1\le n\le12\), the experiment verifies exactly that \(|W_n|=B_n\) and that (10) agrees with direct rational elliptic group arithmetic. It also verifies strong divisibility on a finite range and compares the exact denominators of the terms \(R_j(15)\).

For \(s=15\), the exact denominator bit lengths of \(R_j(s)\), \(2\le j\le12\), are:
\[
2629,\ 5257,\ 10514,\ 15771,\ 23657,\ 31542,\ 42056,\ 52570,\ 65712,\ 78855,\ 94625.
\]
In this finite sample, every pair with \(\gcd(j,k)=1\) has coprime reduced denominators. The exact denominators of the binomial and alternative sums also reproduce the height-cost gap seen in the previous experiment.

## 7. What this does and does not tell us about \(\pi\)

This is a genuine obstruction theorem for the elliptic formal-group route: high-order cancellation cannot make the output rational height arbitrarily cheap, because the largest multiplier forces a large primitive EDS denominator.

It does **not** prove that \(\pi\) is badly approximable. The EDS lower bound controls the height of our constructed rational combinations; it does not yet relate their real distance to \(\pi\) strongly enough to yield a positive lower bound for \(q\|q\pi\|\).

The useful next step is to compare the exact error \(|\mathcal B_{\lambda,s}-1|\) with the forced lower bound (7), and determine whether any choice of \(s,m\) can beat the height cost. This will tell us whether the elliptic construction can ever produce a useful Diophantine approximation input for the continuous \(\pi\)-lattice flow.

## References

- M. Streng, *Divisibility sequences for elliptic curves with complex multiplication*, Algebra & Number Theory 2 (2008), Proposition 1.3, especially the primitive-part asymptotic for the \(\mathbb Z\)-indexed sequence, pp. 184–185.
- K. Eisenstädter and G. Everest, *Descent on elliptic curves and Hilbert's tenth problem*, arXiv:0707.1485, for the local-order interpretation of primitive EDS divisors and the H10 use of elliptic divisibility sequences.
- M. Nowroozi and S. Siksek, *Perfect powers in elliptic divisibility sequences*, Bull. Lond. Math. Soc. (2024), for strong divisibility and local valuation properties.


For the application of Streng's theorem, base-change \(E\) from \(\mathbb Q\) to its CM field \(K=\mathbb Q(\sqrt{-3})\), where the endomorphisms are defined. The point \(P\) remains rational over \(K\), and the denominator ideal of \(x(nP)\) is generated by the rational integer \(B_n\). Its normalized ideal norm is \(B_n\), so the primitive-part size estimate is exactly the stated estimate for the rational EDS (not its square).


## 8. CM-sparse cancellation improves the multiplier bound

The later note CM_SPARSE_ORDER12_CANCELLATION.md exploits the order-six automorphism of \(y^2=x^3-2\). It uses the three weights
\[
(\lambda_2,\lambda_3,\lambda_4)=(-364,291,-36)
\]
to cancel the only non-linear formal coefficient below degree 13. The resulting error is \(O(\xi(V)^{12})\) while the largest multiplier is only \(J=4\). Thus the primitive-denominator lower bound becomes
\[
\liminf_{s\to\infty}\frac{h(\mathcal C([s]P))}{s^2}
\ge \frac{64\widehat h(P)}{\zeta(2)},
\]
rather than the \(324\widehat h(P)/\zeta(2)\) floor for the generic binomial order-seven filter. The exact \(s=41\) comparison is recorded in cm_sparse_cancellation_experiment.json.
