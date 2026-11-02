# Exact first-layer prime budget and audit of the proposed global reduction

2026-11-02. This note independently audits the reduction to the upper-half exceptional-prime sum. The reduction is valid. It can be expressed as an exact integer divisibility before using the prime number theorem. No uniform bound S(n)<=n/6 is proved here.

## 1. Objects

Let A_0=A_1=1, B_0=0, B_1=1, and let each integer sequence obey

\[
X_{n+1}=(2n+1)X_n+n^2X_{n-1}\qquad(n\ge1).
\]

For n>=2 define g_n=gcd(A_n,4B_n), R_n=g_n/n!, and

\[
e_p(n)=v_p(g_n)-v_p(n!)\quad(p\text{ odd prime}),
\]
\[
E_n=\prod_{p\text{ odd}}p^{\max(e_p(n),0)},\qquad
D_n=\prod_{p\text{ odd}}p^{\max(-e_p(n),0)}.
\]

Thus E_n,D_n are coprime odd integers, and the odd part of the positive rational R_n is E_n/D_n. Write odd(m) for the positive odd part of an integer m.

Let

\[
H_n=\prod_{n/2<p<n\atop p\text{ odd prime}}p,
\qquad
P_n=\prod_{n/2<p<n\atop p\mid A_{n-p}}p,
\]
\[
S(n)=\ln P_n,\qquad
L_n=\operatorname{odd}(\operatorname{lcm}(1,\ldots,n-1)).
\]

All products are finite, and an empty product is one. The symbol P_n here denotes this exceptional-prime product, not the endpoint polynomial in older geometric notes.

The existing prime-adic results used below are:

1. For every odd prime p and n>=2,
   \(e_p(n)\le\lfloor\log_p(n-1)\rfloor-v_p(n)\).
2. For n/2<p<n, write r=n-p, so 0<r<p. Then v_p(g_n)=0 if p does not divide A_r; otherwise v_p(g_n) is one or two. In this range v_p(n!)=1.

The first bound follows from the rational Wronskian and the coefficient denominator bound for B_n/n!. The second follows from the mod-p² shifted transfer criterion. They are established proof inputs, not density assumptions.

## 2. Exact integer certificate

For every n>=2,

\[
\boxed{E_n\operatorname{odd}(n)H_n^2
\mid D_nP_n^2L_n.}
\tag{1}
\]

To prove this, compare the p-adic exponents for each odd prime p.

- If p<=n/2, it contributes nothing to H_n or P_n. The required inequality is exactly e_p(n)+v_p(n)<=floor(log_p(n-1)). This includes the possible endpoint p=n/2 when n=2p; no upper-half classification is used at that endpoint.
- If n/2<p<n, it does not divide n and it occurs once in L_n. If p does not divide A_(n-p), its excess is -1, and the right-minus-left exponent in (1) is zero. If it is exceptional, its excess is at most one, and the right-minus-left exponent is at least zero. These two cases are e_p(n)<=-1+2·1_(p|A_(n-p)).
- If p>=n, it occurs in none of H_n,P_n,L_n. The general bound gives e_p(n)<=-v_p(n), which proves the required nonnegative difference. This handles p=n when n is an odd prime.

Every integer in (1) is odd, so there is no further 2-adic condition. The quotient in (1) is therefore a positive integer.

In particular, logarithms give the fully finite bound

\[
\boxed{\ln(E_n/D_n)\le
2S(n)+\psi_{\rm odd}(n-1)-2\ln H_n-\ln\operatorname{odd}(n),}
\tag{2}
\]

where psi_odd(m)=ln odd(lcm(1,...,m))=sum_(p odd,k>=1,p^k<=m) ln p. There is no asymptotic error in (1) or (2).

## 3. Asymptotic cancellation and the user's reduction

The prime number theorem gives psi_odd(n-1)=n+o(n) and

\[
\ln H_n=\vartheta(n-1)-\vartheta(n/2)=n/2+o(n).
\]

The omitted prime 2 changes this by at most a bounded amount at small indices. Since ln odd(n)=O(ln n), (2) proves

\[
\boxed{\ln(E_n/D_n)\le2S(n)+o(n).}
\tag{3}
\]

This is an all-index upper bound. It does not require a subsequence to have any density. The lower-half allowance can be written as the first-prime weight up to n/2 plus the prime-power tail; the latter is O(sqrt(n)(ln n)^2)=o(n) by a direct count, so no unexplained linear contribution from higher powers is present.

It is important that the subtracted upper-half nonexceptional deficits are retained: bounding all positive excess and dropping every deficit would give only a linear upper bound too large for (3).

## 4. Consequence for the quadratic approximation quality

Let Q_n=q_n²|pi-p_n/q_n| for the fully reduced approximant p_n/q_n=4B_n/A_n. The preceding checkpoint proves

\[
Q_n\sim\frac{\sqrt2}{n}R_n^{-2}.
\]

There is also a strict finite certificate, using the actual d_n=v_2(R_n). The all-index positive-integral bound Q_n>2R_n^(-2)/(2n+1), followed by (1), gives

\[
\boxed{Q_n>
\frac{2^{1-2d_n}\operatorname{odd}(n)^2H_n^4}
{(2n+1)P_n^4L_n^2}.}
\tag{4a}
\]

No pi digits, asymptotic remainder, or prime number theorem is needed to evaluate this rational lower bound at a specified index.

The exact binary saturation theorem gives R_n=2^(-n/2+beta_n) E_n/D_n for n>=4, with beta_n=O(ln n). Hence (2) implies

\[
\ln Q_n\ge n\ln2-4S(n)
-2\psi_{\rm odd}(n-1)+4\ln H_n
+2\ln\operatorname{odd}(n)-O(\ln n).
\tag{4}
\]

Applying PNT only to the explicit middle terms yields

\[
\boxed{\ln Q_n\ge n\ln2-4S(n)+o(n).}
\tag{5}
\]

Therefore the proposed eventual inequality S(n)<=n/6 would indeed prove

\[
Q_n\ge\exp((\ln2-2/3+o(1))n)\longrightarrow\infty.
\]

The coefficient is positive. A slightly weaker sufficient global target is available: it would be enough to prove limsup S(n)/n<ln2/4=0.1732867951..., rather than specifically 1/6. The exact equality case is not settled by this exponential comparison.

These implications exclude this constructed sequence if their global S-hypothesis is proved. Neither the hypothesis nor its negation follows from the present budget calculation, PNT, or a finite scan. Excluding the sequence would not decide pi's bad approximability over all reduced rationals.

A bound merely on average would not suffice without further structure. For example, the abstract sequence f(n)=n/4 at powers of two and f(n)=0 elsewhere has the Cesaro average of f(n)/n tending to zero, but limsup f(n)/n=1/4. This is not asserted to be realizable by the present S(n); it illustrates why a density-one or mean-value conclusion alone does not remove arbitrarily sparse exceptional indices. The proposed eventual pointwise lemma would remove them.

## 5. A sharper budget retaining the actual second lift

Optionally let P_n^(2) be the product of upper-half odd primes p with p² dividing g_n. It divides P_n by the local classification. At each such prime the exact identity is

\[
e_p(n)=-1+\mathbf1_{p\mid A_{n-p}}+\mathbf1_{p^2\mid g_n}.
\]

Consequently the stronger divisibility holds:

\[
\boxed{E_n\operatorname{odd}(n)H_n^2
\mid D_nP_nP_n^{(2)}L_n.}
\tag{6}
\]

Its logarithmic version has S(n)+ln P_n^(2) in place of 2S(n). Formula (1) follows on using P_n^(2)|P_n. This refinement separates first-layer anomalies from actual extra-factor gains; it is not needed to validate the proposed simpler target.

In particular, writing T(n)=ln P_n^(2), one has the stronger all-index asymptotic upper bound

\[
\boxed{\ln(E_n/D_n)\le S(n)+T(n)+o(n),\qquad 0\le T(n)\le S(n).}
\tag{7}
\]

The finite quality lower bound (4a) also strengthens on replacing its P_n^4 denominator by P_n²(P_n^(2))².

## Verification

The standard-library companion reconstructs A,B,n!,g and every upper-half exceptional set through a finite specified limit. It checks the local valuations, both divisibility certificates, and the exact lower-half budgets without floating-point logarithms. These finite checks audit arithmetic only; the universal reduction follows from Sections 2–4.

The default run through n=1000 passed 999 checks of each aggregate divisibility and 39,849 exact upper-half classifications: 39,708 nonexceptional deficits, 139 first-layer exceptions without positive excess, and two true extra-factor gains. It also independently certified S(n)<n/6 and S(n)/n<=ln(17)/21 throughout that finite range, with equality in the latter at n=21. These are narrower independent checks, not a claim to rerun the submitted n<=10000 scan.

Proof inputs are the [odd-prime Wronskian cap](https://github.com/epsilon-lain/lain-s-note/blob/db233082d7ee811cccc8d2cbd376a3daec6b89fb/research/pi-gcd-blast-2026-11-02/gcd-prime-blocks/ODD_PRIME_BLOCKS_AND_EXCESS.md), [prime-square classification](https://github.com/epsilon-lain/lain-s-note/blob/db233082d7ee811cccc8d2cbd376a3daec6b89fb/research/pi-gcd-prime-adic-2026-11-02/gcd_blocks/PRIME_BLOCKS_AND_LIFTING.md), and [binary saturation with exponent normalization](https://github.com/epsilon-lain/lain-s-note/blob/db233082d7ee811cccc8d2cbd376a3daec6b89fb/research/pi-error-exponent-2026-11-02/README.md). The classical PNT input is summarized in [NIST DLMF §27.12](https://dlmf.nist.gov/27.12). No claim of a new prime number theorem or density theorem is made.
