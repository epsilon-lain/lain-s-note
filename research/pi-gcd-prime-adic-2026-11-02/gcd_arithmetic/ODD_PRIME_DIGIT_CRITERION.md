# Exact odd-prime digit criterion for maximal denominator loss

2026-11-02. Same recurrence and initial data as TWO_ADIC_VALUATION.md. Write a_n=A_n/n!, b_n=B_n/n!, and g_n=gcd(A_n,4B_n). All congruences of rational numbers below use denominators coprime to the odd prime p after the indicated multiplication.

## 1. Convolution formula

Let Q(z)=1-2z-z² and F(z)=sum a_n z^n=Q(z)^(-1/2). The normalized recurrence gives

Q(z)F'(z)=(1+z)F(z),
Q(z)B'(z)=(1+z)B(z)+1,

where B(z)=sum b_n z^n and B(0)=0. Therefore

B(z)=F(z) integral_0^z F(t)dt,
b_n=sum_{k=1}^n a_{n-k}a_{k-1}/k.

The binomial expression for a_n shows a_n belongs to Z[1/2]. Consequently, for every odd prime p,

v_p(b_n)>=-floor(log_p n),
v_p(g_n)>=v_p(n!)-floor(log_p n),

whenever p<=n. The lower bound is nonnegative automatically. This is a baseline divisibility statement; it is not the upper bound needed to exclude rare unusually large gcd values.

## 2. Lucas congruence for a_n

In F_p[[z]],

F(z)=Q(z)^{(p-1)/2}F(z^p).

The polynomial Q(z)^{(p-1)/2} has degree p-1, with coefficients a_0,...,a_{p-1} modulo p. Thus, if n=sum_{j=0}^e d_jp^j is the base-p expansion,

a_n = product_{j=0}^e a_{d_j} mod p.

Also a_{p-1}=chi_p modp, where chi_p=(-1)^{(p-1)/2}, by taking the highest coefficient of Q(z)^{(p-1)/2}.

## 3. Maximal-pole residue

Let p<=n be odd, e=floor(log_p n)>=1, and n=sum_{j=0}^e d_jp^j with 1<=d_e<p. Then

p^e b_n = chi_p^e b_{d_e} product_{j=0}^{e-1}a_{d_j} modp.

Proof: in the convolution for p^e b_n, all k not divisible by p^e vanish modulo p. The remaining k=mp^e have 1<=m<=d_e<p. Lucas gives

a_{n-mp^e}=a_{d_e-m} product_{j<e}a_{d_j},
a_{mp^e-1}=a_{m-1} a_{p-1}^e.

The sum over m is exactly b_{d_e}, using the same convolution at d_e. This proves the formula.

Whenever the displayed right side is nonzero, it follows exactly that

v_p(b_n)=-e,
v_p(g_n)=v_p(n!)-e.

There is no cancellation caveat in this nonzero case, since a_n is p-integral. Conversely, a zero right side means only that the maximum possible pole has canceled; it does not by itself imply positive gcd surplus.

For p<=n<2p, n=p+r, this simplifies to

p b_n=chi_p a_r modp.

Thus any prime in this range with a_r nonzero modulo p is completely absent from g_n. In particular, when n is an odd prime, n does not divide g_n.

## 4. Exact counterexamples to an attractive but false shortcut

It is false that the odd part of g_n always divides n!. Exact recurrence computation gives odd valuation surplus at least at (n,p)=(72,17), (107,73), (166,131), (170,41), (176,23). These are finite counterexamples, not predictions about distribution. In particular, no argument may replace the missing odd-prime upper bound by that divisibility assertion.

The digit theorem isolates where maximal denominator loss can fail. Bounding the total positive valuation surplus at those exceptional primes remains necessary for the global estimate described in TWO_ADIC_VALUATION.md.

## 5. A uniform finite cap on every lifting exponent

There is also a general upper bound, independent of the digit residue. Let n>=2 and p be any odd prime. Put s_p(n)=min(v_p(a_n),v_p(b_n))=v_p(g_n)-v_p(n!). The normalized exact Wronskian is

a_n b_{n-1}-a_{n-1}b_n=(-1)^n/n.

Write e_p=floor(log_p(n-1)), with e_p=0 when p>n-1. Section1 gives v_p(a_{n-1})>=0 and v_p(b_{n-1})>=-e_p. Each term on the left therefore has valuation at least s_p(n)-e_p. Taking valuations of the identity proves

s_p(n)<=floor(log_p(n-1))-v_p(n).

In particular, positive odd-prime surplus is capped at max(0,e_p-v_p(n)). With L_{n-1}=lcm(1,...,n-1),

E_n divides the odd part of L_{n-1}/gcd(L_{n-1},n).

This also shows that every surplus prime is below n, that primes above sqrt(n-1) can occur in E_n only to the first power, and that a prime p with n=p^k cannot contribute positive p-surplus. These are exact all-index statements. The total logarithmic bound supplied by the lcm is still of order n; by itself it does not overcome the (log2)n/2 loss in the 2-adic theorem.
