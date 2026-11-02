# Integer quarter-period coordinates and a shared-denominator control

2026-11-02. This is a constructive transfer from the supplied elliptic bridge. It retains the symmetric-square period and derivative pairing and measures its actual reduced denominator. The bounded control is not proposed for enlargement: sharing a denominator helps relative costs at two of the tested budgets, but it is not enough to produce useful quadratic-quality approximants.

## 1. A smaller integer second-order coordinate

Clausen's formula, with a=b=1/4, gives

\[
h(u)={}_2F_1(1/4,1/4;1;u),\qquad
H(u)=h(u)^2.
\]

Set t=u/64, and write F(t)=h(64t). Then

\[
F(t)=\sum_{n\ge0}v_nt^n,\qquad
v_n=b_n^2,\qquad
b_n=\frac{2^n\prod_{j=0}^{n-1}(4j+1)}{n!}\in\mathbb Z.
\]

For every odd prime p, each complete p^e-block of j contributes a factor divisible by p^e in 4j+1, cancelling v_p(n!). For p=2 the valuation of b_n is n-v_2(n!)=s_2(n). This proves integrality and the exact recurrence

\[
\boxed{(n+1)^2v_{n+1}=4(4n+1)^2v_n,\qquad v_0=1.}
\]

The first values are 1,4,100,3600,152100. The convolution identity is

\[
\sum_{j=0}^{n}v_jv_{n-j}=\binom{2n}{n}^3.
\]

For theta=t d/dt, the differential equation and first-order system are

\[
[(1-64t)\theta^2-32t\theta-4t]F=0,
\qquad
\theta\binom F{\theta F}=
\begin{pmatrix}0&1\\4t/(1-64t)&32t/(1-64t)\end{pmatrix}
\binom F{\theta F}.
\]

At t0=1/4096, the classical Ramanujan identity supplied by the user becomes the explicit bilinear pairing

\[
\boxed{F(t_0)[5F(t_0)+84\theta F(t_0)]=16/\pi.}
\]

This uses an analytic period F; a fundamental matrix with a logarithmic companion still needs its cusp normalization. An integer monodromy matrix does not normalize that companion or remove an unaccounted 2 pi i factor.

## 2. A common-denominator construction with a fixed source budget

Fix s>=1. Choose Q,P0,P1 of degree at most 2s so that

\[
QF-P_0=O(t^{3s+1}),\qquad
Q\theta F-P_1=O(t^{3s+1}),\qquad Q(0)=1.
\]

There are 2s rational linear equations for the nonconstant coefficients of Q. The accompanying program solves these exactly when the system is nonsingular, reconstructs both numerators, and clears the content of all three polynomials together. Only v_0,...,v_(3s) are used for generation. No value of pi or target rational is supplied.

For a common degree r=2s put P_(j,*)=4096^r P_j(1/4096), Q_*=4096^r Q(1/4096), after the integer polynomial normalization. All are integers. The generated rational is

\[
\boxed{\frac pq=\frac{16Q_*^2}{P_{0,*}(5P_{0,*}+84P_{1,*})}}
\]

followed by complete gcd reduction. In particular its reduced denominator is exactly the absolute denominator displayed above divided by the gcd with 16Q_*^2. The homogenization powers cancel before reduction; they are not separately charged a second time.

## 3. A fully rational achieved-quality certificate

For n>=0, the positive terms of F at t0 have successive ratio less than 1/64. For theta F, beyond a positive index, the successive ratio is less than 1/32. Thus

\[
\sum_{n>k}v_nt_0^n\le\frac{64}{63}v_{k+1}t_0^{k+1},\qquad
\sum_{n>k}nv_nt_0^n\le\frac{32}{31}(k+1)v_{k+1}t_0^{k+1}
\quad(k\ge1).
\]

Let q_j be the integer coefficients of Q, N=3s, and r_i=P_i(t0)/Q(t0). The two error radii are bounded by

\[
e_i=\frac1{|Q(t_0)|}\sum_{j=0}^{2s}|q_j|t_0^j\,B_i(N-j),
\]

where B0,B1 are the two rational tail bounds above. Then put

\[
C=5r_0^2+84r_0r_1,
\quad E=10|r_0|e_0+5e_0^2+
84(|r_0|e_1+|r_1|e_0+e_0e_1).
\]

If C>E>0, the candidate p/q=16/C has the joint upper certificate

\[
\boxed{q^2|\pi-p/q|\le q^2\frac{16E}{C(C-E)}.}
\]

All quantities are rational and include the actual reduced q. One additional recurrence coefficient is used to certify the tails; it does not enter the generation equations. The source budgets below count construction coefficients, as with the parent's post-construction error certification.

## 4. Shared Q is not exact differential compatibility

The polynomial

\[
\boxed{D=Q\theta P_0-P_0\theta Q-P_1Q}
\]

is divisible by t^(3s+1), and D/t^(3s+1) has degree at most s-1. This follows by differentiating the first approximation identity and substituting the second. The program records this exact low-degree residual.

It is important to keep this defect: P1/Q is not identically theta(P0/Q). If P0,Q are coprime and Q(0)!=0, exact compatibility with the same polynomial denominator would require Q to divide theta Q, hence theta Q=kQ. Evaluation at zero gives k=0, so Q must be constant. Nonconstant rational periods therefore need a higher-power derivative denominator or an explicit residual such as D. This elementary restriction is specific to the proposed shared-Q parametrization; it is not a general limitation on period-pair constructions.

## 5. Bounded fair controls

The program compares four fixed source budgets with scalar Padé of the Ramanujan series S(t)=sum(42n+5)binom(2n,n)^3t^n, using exactly the same number of source coefficients for each generated candidate.

| Source coefficients | Shared pair Q quality (display) | Scalar Padé Q quality (display) |
|---:|---:|---:|
| 4 | 7.89e17 | 9.16e10 |
| 7 | 1.56e42 | 3.30e29 |
| 10 | 5.05e56 | 1.31e64 |
| 13 | 4.17e77 | 4.30e105 |

The shared pair reduces the coefficient/denominator cost in the last two controls but remains a poor approximant. The exact data retain both normalized coefficient heights, both actual reduced fractions, Ramanujan-tail quality intervals and the separate period-tail upper certificate. No larger generic Hermite-Padé scan is recommended from these controls.

What should be retained is the smaller integer second-order system, the perfect-square coefficient normalization, the bilinear Ramanujan pairing, and the explicit low-degree compatibility residual. A further construction must control that residual together with the common denominator, rather than treating common Q as full differential preservation.

## Reproduction and sources

Run `python period_pair_control.py`. Python standard library only. The exact code checks the coefficient convolution, approximation defects, compatibility residual divisibility, complete gcd reduction, and rational period-tail quality certificate in the four stated cases.

- NIST DLMF16.12.2, Clausen's formula: https://dlmf.nist.gov/16.12 . Used with a=b=1/4.
- Jesus Guillera, *A method for proving Ramanujan series for 1/pi*, arXiv:1807.07394, formula(2): https://arxiv.org/pdf/1807.07394 . Supplies the same classical special-value identity as the parent bridge.
- The parent repository's pi-s6-period-bridge-2026-11-02/periods/MODULAR_SUBSYSTEM.md supplies the modular and differential normalization context. This note neither uses nor asserts a new theorem about the global topology of S6.
