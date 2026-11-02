# A quantitative limitation of generic recurrence-window witnesses

2026-11-02. This note examines a specific proposed route to the moving-prime bound: build a common integer witness from finitely many nearby recurrence entries, show that every exceptional prime divides it, and bound its absolute value. The result below explains a limitation of witnesses whose vanishing uses only the central zero condition. It does **not** prove the desired global S(n) estimate, and does not exclude witnesses that exploit the actual initial values in a further arithmetic way.

## 1. Fixed sequence and the intended prime product

Let A_0=A_1=1 and A_(n+1)=(2n+1)A_n+n²A_(n−1). Put

\[
U_n=2^n A_n/n!,\qquad C_n=\operatorname{odd}(U_n).
\]

The U_n are positive integers and satisfy

\[
(n+1)U_{n+1}=2(2n+1)U_n+4nU_{n-1}.
\tag{1}
\]

Here odd(m)=m/2^(v_2(m)) for a positive integer m; for a nonzero signed integer use oddabs(m)=odd(|m|).

The upper-half exceptional prime product P(n)=exp(S(n)) divides C_n, by the existing Schur congruence. The sharp asymptotic of A_n/n! and the exact binary valuation give

\[
\ln C_n=\ln(2+\sqrt2)n+O(\ln n).
\tag{2}
\]

Its exponential coefficient is about 1.227947, while the sufficient target for ln P(n) has coefficient below ln2/4, about 0.173287. A height argument using C_n itself is therefore too large.

## 2. Formal recurrence window

Use an indeterminate N and two independent variables x,y. Set V_0=x, V_(-1)=y, and propagate the same recurrence:

\[
(N+j+1)V_{j+1}=2(2N+2j+1)V_j+4(N+j)V_{j-1}
\quad(j\ge0).
\]

For h>=0 define D_h(N)=∏_(t=1)^h(N+t), with D_0=1. Each V_j for −1<=j<=h is a linear polynomial in x,y whose coefficients lie in D_h(N)^(-1) Z[N]. This follows by induction; the denominator of V_j divides D_j and D_j divides D_h.

Thus the following theorem has an explicit denominator rather than an unspecified elimination cost.

## 3. Common-witness theorem

Let h,d be nonnegative integers. Let

\[
F\in\mathbb Z[N,Z_{-1},Z_0,\ldots,Z_h]
\]

have total degree at most d in the Z variables; its degree and coefficients in N can be arbitrary. Suppose its vanishing is a consequence of the central zero alone, in the precise sense

\[
F(N,V_{-1}(N,0,y),\ldots,V_h(N,0,y))=0
\quad\text{in }\mathbb Q(N)[y].
\tag{3}
\]

For an actual integer n>=1, evaluate F_n=F(n,U_(n−1),U_n,...,U_(n+h)). Then, whenever F_n is nonzero,

\[
\boxed{
\frac{C_n}{\gcd(C_n,D_h(n)^d)}\ \mid\ \operatorname{oddabs}(F_n).
}
\tag{4}
\]

Proof. Clearing the common window denominator gives

\[
D_h(N)^d F(N,V_{-1},\ldots,V_h)\in\mathbb Z[N,x,y].
\]

By (3), its coefficient of x^0 is identically zero. Therefore it equals x H(N,x,y) for an integer polynomial H. Evaluate at N=n,x=U_n,y=U_(n−1):

\[
D_h(n)^d F_n=U_n H(n,U_n,U_{n-1}).
\]

Hence U_n divides D_h(n)^d F_n. At each odd prime, cancellation of at most the factors in D_h(n)^d gives (4). This proof uses no root-count hypothesis, random model, or asymptotic estimate. ∎

The exceptional primes which divide D_h(n) can be discarded separately when applying a modular witness. For n/2<p<n, p|n+t with 1<=t<=h implies n+t=mp for an integer m; no claim that all such primes automatically survive denominator clearing is needed for (4).

## 3a. Linear witnesses retain the entire central value

For linear witnesses the denominator allowance can be removed altogether, for every window length. Suppose F has degree at most one in its window variables and satisfies (3). Then there is a polynomial G(N) in Z[N] such that

\[
\boxed{F(N,V_{-1},\ldots,V_h)=G(N)x.}
\tag{4b}
\]

In particular, F_n=G(n)U_n exactly, and every nonzero such witness has oddabs(F_n)>=C_n. Increasing the window length does not change this conclusion for this class of generic linear witnesses.

Here is a proof of the integral-polynomial assertion, beyond the rational factorization. Write

\[
V_j=(a_j(N)x+b_j(N)y)/D_j(N).
\]

The initial pairs are (a_0,b_0)=(1,0), a_1=2(2N+1), b_1=4N. For j>=1 both numerator polynomials satisfy

\[
c_{j+1}=2(2N+2j+1)c_j+4(N+j)^2c_{j-1}.
\]

All possible poles of the rational coefficients are simple and lie at N=-t, 1<=t<=h. At N=-t, the recurrence first becomes singular in V_t. Its y-residue is nonzero: b_1(-t)=-4t, and the displayed recurrence proves (-1)^j b_j(-t)>0 for 1<=j<=t. Both summands have the asserted sign, since 2(-2t+2j+1)<0 for j<t.

Every residue vector of V_j, j>=t, is a scalar multiple of the residue vector of V_t. Indeed V_(t-1) has zero residue, and subsequent propagation uses a nonsingular denominator N+j+1 for j>=t. Thus, in any polynomial-coefficient linear combination, cancelling the y-residue also cancels the x-residue. Hypothesis (3) cancels the full y coefficient, so the surviving x coefficient G(N) has no finite pole. Hence G is in Q[N]. Since D_h G is in Z[N] and D_h is monic, ordinary polynomial division shows G is actually in Z[N]. The constant term of F is zero by (3), completing the proof.

This conclusion is stronger than a numerical height estimate, but its hypothesis is essential. It concerns identities vanishing for arbitrary y when x=0, not extra congruences known only at the actual neighboring sequence value.

## 4. Quantitative cost

For a nonzero witness, the positive integer on the left of (4) cannot exceed oddabs(F_n). Since D_h(n)<=(n+h)^h, one obtains

\[
\boxed{
\ln\operatorname{oddabs}(F_n)
\ge\ln C_n-dh\ln(n+h)
=\ln(2+\sqrt2)n-dh\ln(n+h)-O(\ln n).
}
\tag{5}
\]

This remains a valid estimate when h,d, and the choice of F depend on n, provided each chosen polynomial satisfies its own formal identity (3). In particular, if dh ln(n+h)=o(n), every nonzero witness in this class still has odd height at least exp((ln(2+sqrt2)+o(1))n). Fixed h,d are a special case.

Consequently, such a witness cannot by itself provide the desired much smaller common-integer upper bound on P(n). It merely retains essentially the entire odd numerator C_n, including prime factors outside the exceptional interval. This statement is about the small-height-witness route: it does not say that every use of a large witness, together with separately proved factor removals, must fail.

For example, the simplest two-neighbor eliminant is

\[
(n+1)U_{n+1}-4nU_{n-1}=2(2n+1)U_n.
\]

All exceptional primes divide it, but the exact right side shows that no smaller height has been obtained. Higher-degree universal eliminants have the same central-factor obstruction, with the explicit denominator allowance in (5).

If F_n=0, it is divisible by every prime and supplies no upper bound by absolute size. The theorem deliberately excludes that case rather than interpreting zero as an advantageous small nonzero witness.

## 5. What is outside this limitation

The argument does not cover an identity whose reduction at x=0 only vanishes for the particular neighbor value selected by the initial data A_0=A_1=1. Such an identity would add arithmetic information beyond the generic recurrence window. Nor does the small-height obstruction for nonlinear witnesses cover prime-dependent polynomials of unaccounted degree or windows whose denominator allowance is linear or greater in n. The stronger linear conclusion in Section 3a applies to every finite window length. Separately justified removal of factors from a large witness is outside both height limitations.

Those are exactly the kinds of additional structure that a useful cross-modulus product bound would have to supply. The existing moving-diagonal estimate S(n)<=n/6 remains unproved. The theorem above should not be presented as a general impossibility theorem for resultants, for the recurrence, or for rational approximation to pi.

## Inputs and verification

The only analytic input to (5) is (2), derived from the previously proved coefficient asymptotic and binary valuation. The divisibility theorem (4) is elementary polynomial algebra and holds for any particular solution of (1) consisting of integers.

The companion standard-library verifier constructs the formal linear window at exact integer base indices, checks both coefficient recurrences against the actual U sequence, and checks divisibility and odd-part cancellation for linear and quadratic sample witnesses. It also verifies 80 nonzero y-residue directions and 4,840 collinearity identities at the possible negative-integer poles. Those finite checks audit the algebra; the universal assertions follow from the proofs in Sections 3 and 3a.
