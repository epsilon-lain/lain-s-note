# Rational-basis truncation without the Taylor factor 3^K

2026-11-02. This records a constructive kernel control suggested by the parent agent and checked here. It gives a real arithmetic change and exact cancellation certificates; the bounded controls below do not recommend enlarging this version of the search.

## 1. Kernel and elementary geometric parameter

For 0<=x<=1,

\[
1+x^2=(1+x)(1-z(x)),\qquad
z(x)=\frac{x(1-x)}{1+x}\in[0,3-2\sqrt2].
\]

The rational upper bound z<=7/40 is sufficient for exact certificates: the polynomial x^2-(33/40)x+7/40 has negative discriminant and positive leading coefficient.

For integers m>=1,c=4m+1,K>=0 put

\[
w_m(x)=x^{2m}(1-x)^{2m},\qquad
S_K(x)=\sum_{k=0}^{K}\binom{c+k-1}{k}
\frac{x^k(1-x)^k}{(1+x)^{c+k}}.
\]

This is a positive rational truncation of (1+x^2)^(-c) on [0,1]. It introduces no coefficient denominator 3^K. Its single real pole is at x=-1.

## 2. Its moments contain rational and log(2) parts

For j=0,1,2 define

\[
\mu_j=\int_0^1 w_m(x)x^jS_K(x)\,dx=a_j+b_j\ln2,
\qquad a_j\in\mathbb Q,\quad b_j\in\mathbb Z.
\]

Each pair is computed exactly by u=1+x. For an individual term with A=2m+j+k,B=2m+k,D=c+k, expand

\[
(u-1)^A(2-u)^B=\sum_l n_lu^l.
\]

Its rational integral is sum over l!=D-1 of
n_l(2^(l-D+1)-1)/(l-D+1), and its log(2) coefficient is n_(D-1).

In particular a common denominator of all a_j,j<=2 divides

\[
\boxed{2^{c+K-1}\operatorname{lcm}(1,\ldots,c+K-1).}
\]

The b_j are integers. The disappearance of 3^K is therefore exact, rather than a change of overall scaling. However the rational and log coefficients can be large with opposing signs even when the actual positive moment is small; both coefficient sizes must be charged.

Take the cross product h=a cross b, clear denominators, and divide its integer content. If h!=0, the quadratic H(x)=h_0+h_1x+h_2x^2 satisfies

\[
\boxed{\int_0^1w_m(x)H(x)S_K(x)\,dx=0}
\]

because its rational and log(2) coefficients vanish separately. This is an explicit candidate generator without pi or a target rational as input.

## 3. The final quadratic-pole integral has no log term under the degree budget

Set

\[
I=\int_0^1\frac{w_m(x)H(x)}{(1+x^2)^c}\,dx.
\]

Although the numerator is not even, its degree is at most4m+2<=2c-2. Its odd part is xP(x^2), with deg P<=c-2. Under t=x^2,u=1+t, all exponents in the odd integral are at most -2. Thus the odd part has a rational primitive and contributes no log(2). The even part has rational and pi periods. Consequently

\[
I=\beta+\alpha\pi/4,\qquad\alpha,\beta\in\mathbb Q.
\]

This degree restriction is essential: arbitrary higher-degree odd numerators would introduce a final log(2) coefficient. The implementation checks the degree and each potential logarithmic exponent explicitly.

Whenever alpha!=0, extract the fully reduced p/q=-4beta/alpha. No general all-scale nonzero-alpha theorem is claimed for the cross-product rule.

## 4. Signed remainder and actual denominator certificate

Put delta=7/40 and suppose

\[
r=\delta\frac{c+K+1}{K+2}<1.
\]

Then the omitted positive tail is at most

\[
\varepsilon=\binom{c+K}{K+1}\frac{\delta^{K+1}}{1-r},
\]

because its successive term ratios decrease and (1+x)^(-c)<=1. With h=sum|h_j| and

\[
M=\int_0^1w_m(x)\,dx=\frac{((2m)!)^2}{(4m+1)!},
\]

the exact truncation cancellation yields

\[
\boxed{|I|\le hM\varepsilon,\qquad
q^2|\pi-p/q|\le\frac{4q^2hM\varepsilon}{|\alpha|}.}
\]

This retains the actual gcd-reduced denominator q. The finite implementation also computes a separate Machin interval only after generation, to distinguish a loose bound from bad actual quality.

## 5. Bounded controls and a higher-pole rational-only variant

The quadratic cross-product control completed m=1,2,3 and K in {m,2m,c,2c}:12 cases, all with nonzero pi coefficient. Its smallest actual quality was about184.1367 at m=1,K=1. The other completed cases were worse. These observations are not a global exclusion theorem and do not justify a larger version of the same scan.

There is a variant where every truncated moment is rational. Set c=4m+K+3 and use j=0,1. Every individual single-pole rational integrand then has numerator degree at most its denominator degree minus two, so its log coefficient is identically zero. The linear multiplier H=mu_1-mu_0 x cancels the entire truncation. This does remove the extra log-cancellation condition, but raises the pole order. Nine controls m=1,2,3,K in {2m,4m,8m} were also poor; their smallest actual quality was about4.75e8.

`rational_basis_truncation.py` generates these controls and exact certificates. The rational-only variant is selected by --rational-only. `free_rational_basis.py` additionally checks a small growing-degree saturated integer-kernel control using SymPy Smith transforms and exact integer LLL. Most short vectors had zero pi period; no tested free control had Q<1. This illustrates why short-vector existence alone cannot guarantee a short nonzero period. Those scripts are bounded controls, not proposed new main scans.

The useful retained mechanism is the exact denominator replacement and two-period cancellation. A next redesign would have to control the separately large rational/log moment coefficients, retain useful multiplier freedom, or change the geometry again. No conclusion about pi bad or nonbad approximability follows from these controls.

## 6. Period-aware saturated-kernel follow-up

A separate, predeclared feature control addressed the zero-period short-vector issue. It used m=1,2,3, K=2m, d=3m+1; saturated moment-kernel dimensions were3,6,9. For each scale the feature map evaluated the final integrand at M=max(kernel_dimension+2,c+2,9) equally spaced interior rational points. It added alpha and beta with five fixed weights2^(-k), k in {0,8m,12m,16m,20m}.

For each feature block all denominators were cleared by one common integer D, and the original integer H coefficients were appended unscaled as guard coordinates. Every reduced LLL row was reconstructed from those guards; all rational samples, period coordinates, and both moment constraints were then checked exactly. The LLL nearest-integer operation used exact numerator//denominator, with SymPy's exact final LLL assertions retained.

The complete bounded run took about6.3seconds. There were24 nonzero-pi output rows and66 zero-pi rows across the fifteen feature blocks. The best observed qualities at m=1,2,3 were respectively about6694.11,179752.94,154626732.31. None had Q<1. This is a different selection objective from plain coefficient LLL, but the completed control still gives no reason to enlarge this particular search.

Files: `period_aware_rational_basis.py` and its JSON. It imports the exact helpers in `free_rational_basis.py` and `rational_basis_truncation.py`; include those dependencies if preserving the control.
