# Independent audit of the multiple-band certificate and witness obstruction

2026-11-02. This review checks the final multiple-band, generic-window and midpoint notes against their stated fixed recurrence. **The claims pass within their stated scope. The eventual bound on S(n), the required primorial-gcd bound, and exclusion of all successful subsequences remain unproved.** This is an independent mathematical and computational check within the project, not external peer review.

## 1. Exact multiple-band certificate

Write ell_n=A_n/n!, w_n=B_n/n!, R_n=g_n/n!=2^{d_n}E_n/D_n, and C_n=odd(2^n ell_n), with the notation of `multi-band-budget/MULTI_BAND_DEFICIT_CERTIFICATE.md`.

For a prime p dividing H but not P F, write n=bp+r. The cutoff ensures b<=B and b<p. Lucas gives ell_n=ell_b ell_r modulo p. Since p does not divide C_n, ell_r is nonzero. Also p does not divide B_b, and b! is a p-unit. The leading-pole congruence therefore gives a nonzero residue for p w_n. Hence v_p(w_n)=-1, while ell_n is a p-unit, so e_p(n)=-1 exactly. Squarefreeness now proves

\[
M=H/\gcd(H,PF)\mid D_n,
\qquad H\mid D_nPF.
\]

This includes r=0 and p=n. Retaining the overlap with gcd(H,PF), instead of dividing numerically by PF, is valid and sharper.

The refined gain divisor Ghat is also valid. Positive excess requires the entire alpha=v_p(C_n) to satisfy both neighboring pole caps. If alpha exceeds either cap, the prime contributes no positive excess; deleting its whole contribution is justified. Otherwise Ghat retains p^alpha, which contains the possible positive excess. Consequently E_n divides Ghat_n. Since d_n<=s_n-n and D_n>=M,

\[
R_n\le 2^{s_n-n}\widehat G_n/M,
\qquad
Q_n>\frac{2\,2^{2(n-s_n)}M^2}
{(2n+1)\widehat G_n^2}.
\]

This is an unconditional improvement by M squared over the earlier gain-only floor. It is not itself a global asymptotic exclusion.

An independently reproduced example is n=21, B=2:

| Object | Exact value |
|---|---:|
| H | 46189 |
| P | 17 |
| F | 1 |
| M | 2717 |
| Ghat | 17 |
| s_n | 13 |
| Earlier floor | 131072/12427 |
| Improved floor | 967585169408/12427 |
| Improvement factor | 7382089 |

Both floors already exceed one, so this example demonstrates the stronger certificate without claiming a newly excluded index.

## 2. Asymptotic concentration condition

The bound B_b<=b! (1+sqrt(2))^b gives ln F=O(B^2 ln B). With B=floor(n^(1/3)), this is o(n), and the cutoff t is asymptotic to n^(2/3). The prime number theorem gives ln H=n+o(n). Every selected prime can contribute positive exponent at most one and must then divide C_n; thus E_high divides P. The remaining positive weight is at most floor(t) ln(n-1)=o(n). These facts prove

\[
J_n\le 2\ln P-n+o(n),
\qquad
\ln Q_n\ge (2+\ln2)n-4\ln P+o(n).
\]

Therefore bounded Q_n along unbounded indices requires

\[
\ln P\ge (1/2+\ln2/4)n-o(n).
\]

There is no proved uniform upper bound for this P. In particular, its trivial height bound has a larger exponential coefficient and cannot close the argument.

One minor fixed-B proof wording was corrected during review: the low-prime contribution is bounded by theta(t)+O(sqrt(n) ln n), not literally by psi(t). Prime powers between t and n can contribute. The correction leaves the stated asymptotic fixed-B inequality unchanged.

## 3. Generic recurrence-window witnesses

For the formal recurrence starting from V_0=x,V_{-1}=y, clearing D_h^d places a degree-d witness substitution in Z[N,x,y]. Generic vanishing at x=0 means the cleared polynomial is divisible by x. Evaluation at the integer sequence proves

\[
C_n/\gcd(C_n,D_h(n)^d)\mid\operatorname{oddabs}(F_n)
\]

whenever F_n is nonzero. The stated quantitative height consequence follows from D_h(n)<=(n+h)^h and the established asymptotic for C_n.

The stronger linear result also passes. At each possible pole N=-t, the numerator b_t(-t) is nonzero by the sign induction in the note. All subsequent residue vectors lie in the same one-dimensional direction with nonzero y component. Cancelling the y coefficient therefore cancels every x pole. The resulting coefficient G is polynomial over Q; because D_h G is integral and D_h is monic, polynomial division gives G in Z[N]. Thus F_n=G(n)U_n exactly.

The arbitrary-y hypothesis is essential. These results do not exclude additional identities using actual initial data, justified removal of large factors, or all possible elimination methods. A zero witness does not give an upper bound by absolute size.

## 4. Edge capacity and midpoint arithmetic

The folded index k=min(n-p,2p-n-1) has at most two possible primes for each k. Hence the edge weight is at most 2(K+1) ln n. The endpoint-interval formulation and the prime number theorem further give o(n) for any K=o(n); this does not require a prime theorem with error o(K) in short intervals.

For the curve y^2=x(x^2+2x+2), Euler's criterion and coefficient extraction give a_p congruent to U_((p-1)/2) modulo p. Hasse's bound makes vanishing equivalent to a_p=0 for p>=5. Rational two-torsion excludes p=1 modulo 4; a rational point of order four then forces 2 to be a square, giving the necessary condition p=7 modulo 8. The argument proves neither the observed stronger residue pattern nor a moving-index zero-density bound. A fixed n has at most one exact midpoint prime.

## 5. Independent finite verification

`verify_crossprime_independent.py` uses exact integers and rational fractions and imports no author verifier. Its recorded run checks n through 500, using several band counts at each index:

| Check | Count |
|---|---:|
| Multiple-band divisibility and overlap certificates | 2431 |
| Exact forced prime deficits | 82305 |
| High-gain support and low-prime capacity | 2431 |
| Computable improved quality floors | 2431 |
| Linear and quadratic window witnesses | 2520 |
| Elliptic midpoint point counts | 94 |
| Edge affine-capacity checks | 1468 |

All passed. Linear witnesses also satisfy exact U_n divisibility. The midpoint zero primes in this finite range are 31, 223, 383 and 479. These checks audit identities and implementation; no extrapolation from finite samples supplies the missing global theorem.
