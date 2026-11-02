# The 2-adic part of the gcd is fully saturated, except at index three

2026-11-02. For the fixed recurrence

\[
A_0=A_1=1,\quad B_0=0,\ B_1=1,\qquad
X_{n+1}=(2n+1)X_n+n^2X_{n-1},
\]

put g_n=gcd(A_n,4B_n), r_n=4B_n/A_n=p_n/q_n in lowest terms with q_n>0. The conjectural extra 2-adic loss in the supplied research package can be removed completely:

\[
\boxed{v_2(g_n)=v_2(A_n)\quad(n\ne3),
\qquad v_2(g_3)=v_2(A_3)-1.}
\tag{1}
\]

Equivalently, q_n is odd at every index except n=3, where q_3=6. The statement includes n=0. No exact formula for v_2(B_n) is needed.

## 1. Adjacent rational differences have an exact valuation

The established integer Wronskian gives, for n>=1,

\[
r_n-r_{n-1}
=\frac{4(-1)^{n+1}((n-1)!)^2}{A_nA_{n-1}}.
\]

Write f(n)=v_2(A_n), and let s_2(m) count the ones in the binary expansion of m. Legendre's factorial identity yields

\[
d_n:=v_2(r_n-r_{n-1})
=2+2\bigl(n-1-s_2(n-1)\bigr)-f(n)-f(n-1).
\tag{2}
\]

Use the all-index theorem already proved in the input package:

\[
f(4m)=2m,\quad f(4m+1)=2m,\quad f(4m+2)=2m+2,
\]
\[
f(4m+3)=2m+3+v_2(m+1).
\]

The binary digit identities s_2(2m)=s_2(m), s_2(2m+1)=s_2(m)+1, and
s_2(m-1)=s_2(m)-1+v_2(m) give

\[
\boxed{
\begin{aligned}
d_{4m}&=4m-3-2s_2(m)-3v_2(m),&&m\ge1,\\
d_{4m+1}&=4m+2-2s_2(m),&&m\ge0,\\
d_{4m+2}&=4m-2s_2(m),&&m\ge0,\\
d_{4m+3}&=4m-1-2s_2(m)-v_2(m+1),&&m\ge0.
\end{aligned}}
\tag{3}
\]

In particular d_3=d_4=-1.

## 2. All later adjacent differences are 2-adic integers

We verify d_n>=0 for every n>=5 directly from (3).

For the second and third classes, s_2(m)<=m suffices. For the fourth class m>=1; the elementary inequality v_2(m+1)<=m gives

\[
d_{4m+3}\ge4m-1-2m-m=m-1\ge0.
\]

For the first class m>=2, write m=2^t u with u odd. Then s_2(m)=s_2(u)<=u, so

\[
d_{4m}\ge u(4\cdot2^t-2)-3(t+1).
\]

If t=0, then u=m>=3 and this is at least 2u-3>0. If t>=1, it is at least
4*2^t-3t-5. This is zero at t=1 and strictly increases thereafter. Hence the first class is nonnegative as well.

This is an all-index sign proof, not an extrapolation from the finite valuation table.

## 3. A single base value proves saturation

The first values are

\[
r_0=0,\quad r_1=4,\quad r_2=3,\quad
r_3=19/6,\quad r_4=160/51.
\]

In particular r_4 has odd denominator. If r_(n-1) has odd denominator and d_n>=0, the sum r_n=r_(n-1)+(r_n-r_(n-1)) has odd denominator after reduction. This is just closure under addition of the ring of rational numbers with odd denominator.

Starting at n=4 and using Section 2 proves that every r_n for n>=4 has odd denominator. The indices 0,1,2 are immediate; index 3 has exactly one power of two in its denominator. Since q_n=A_n/g_n, this is precisely (1).

Thus an extra-loss variable t_n defined by t_n=v_2(A_n)-v_2(g_n) satisfies

\[
\boxed{t_n=0\ (n\ne3),\qquad t_3=1.}
\tag{4}
\]

## 4. Exact use in the exponent budget

For n>=1 define

\[
E_n=\prod_{p\text{ odd}}p^{\max(0,v_p(g_n)-v_p(n!))},
\quad
D_n=\prod_{p\text{ odd}}p^{\max(0,v_p(n!)-v_p(g_n))}.
\]

These are positive coprime odd integers. At every n>=4,

\[
\boxed{\frac{g_n}{n!}
=2^{f(n)-n+s_2(n)}\frac{E_n}{D_n}.}
\tag{5}
\]

The power of two is now known exactly, rather than merely bounded above. Its exponent is -n/2+O(log n), by the explicit formula for f(n). Consequently there is no unknown linear-scale 2-adic loss to include when translating approximation error into a power of the reduced denominator. The remaining unknown exponential contribution lies in the net odd-prime ratio E_n/D_n, not in E_n alone.

This theorem does not bound that odd-prime ratio. It therefore removes one uncertainty from the exponent budget without producing a successful approximation exponent or deciding pi's bad approximability.

## Verification

`verify_prime_lift_and_saturation.py` checks the four difference-valuation formulas, the precise exceptions, and actual gcd saturation against exact integer recurrences through n=5000. It also independently audits the new prime-square lifting rule. All mathematical checks run by default without the original package. Passing `--source` with its exact directory additionally validates `GCD_RESEARCH_MANIFEST.json`; `--output` selects a separate result file. The universal saturation theorem follows from the sign argument and the base value above.
