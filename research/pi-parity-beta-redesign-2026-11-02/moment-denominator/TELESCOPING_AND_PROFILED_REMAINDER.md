# One differential residual replaces the full Taylor error

2026-11-02. This continues the prescribed complete-truncation construction. It derives an exact telescoping identity for the two-moment denominator combination, and a stronger signed-integral certificate which retains where the Taylor remainder is concentrated. The generated fractions and their reduced denominators stay the same; the new gain is in the proved error certificate, with no new multiplier cost.

## 1. Exact differential residual

Let c>=1, K>=0, f(t)=(1+t)^(-c), and

\[
S_K(t)=\left(\frac23\right)^c\sum_{k=0}^{K}
(-1)^k\binom{c+k-1}{k}\left(\frac{2t-1}{3}\right)^k.
\]

Define the rational scalar

\[
\boxed{\gamma_K=
\frac{c\,2^c(-1)^K\binom{c+K}{K}}{3^{c+K}}.}
\]

The binomial recurrence cancels every coefficient except the highest one in

\[
\boxed{(1+t)S_K'(t)+cS_K(t)=\gamma_K(2t-1)^K.}
\tag{1}
\]

Since S_K(1/2)=f(1/2), solving this first-order equation gives the exact error, not a tail estimate:

\[
\boxed{S_K(t)-f(t)=\frac{\gamma_K}{(1+t)^c}
\int_{1/2}^{t}(1+s)^{c-1}(2s-1)^K\,ds.}
\tag{2}
\]

## 2. The remainder retains an extra high power

For 0<=t<=1, put u=|2t-1|. After changing variables in (2), the absolute integral is of the form

\[
2^{-c}\int_0^u v^K(3\pm v)^{c-1}\,dv.
\]

The derivative of v^(K+1)(3-v)^(c-1) is at least
(K+1-(c-1)/2)v^K(3-v)^(c-1) on [0,1]. For the plus sign it is at least (K+1)v^K(3+v)^(c-1). Both lower bounds are positive whenever 2K-c+3>0. Consequently, under that weaker condition,

\[
\boxed{|S_K(t)-f(t)|\le
\frac{|\gamma_K|}{2K-c+3}\,
\frac{|2t-1|^{K+1}}{1+t}.}
\tag{3}
\]

The previous uniform epsilon_K bound discards the factor |2t-1|^(K+1). In the present construction the endpoint weight is small precisely where that factor is largest, so retaining it is useful.

## 3. A fully rational weighted certificate

Return to a=m,b=2m,c=4m+1,K=13m and w(x)=x^(2m)(1-x^2)^(2m). Let mu_j be the exact truncation moments from the preceding checkpoint and H(t)=mu_1-mu_0t. Then integral w(x)H(x^2)S_K(x^2) dx=0.

Set E=2 floor((K+1)/2), an even integer satisfying K<=E<=K+1, and define positive rational numbers

\[
Z_j=\int_0^1x^{2m+2j}(1-x^2)^{2m}(2x^2-1)^E\,dx,
\qquad j=0,1.
\]

They are explicitly computable by

\[
Z_j=\sum_{l=0}^{E}(-1)^{E-l}\binom El2^l M_{j+l},
\quad M_r=\frac{2^{2m}(2m)!}{\prod_{s=0}^{2m}(2m+2r+2s+1)}.
\]

Because |2t-1|^(K+1)<=|2t-1|^E and 1/(1+t)<=1, (3) proves

\[
\boxed{|I_m(H)|\le\frac{|\gamma_K|}{2K-c+3}
\bigl(|\mu_1|Z_0+|\mu_0|Z_1\bigr).}
\tag{4}
\]

Thus, for the same actual reduced denominator q and pi coefficient alpha as before,

\[
\boxed{Q\le\frac{4q^2}{|\alpha|}
\frac{|\gamma_K|}{2K-c+3}
\bigl(|\mu_1|Z_0+|\mu_0|Z_1\bigr).}
\tag{5}
\]

Every quantity on the right is rational. It costs no new coefficient denominator: only the certification uses the rational Z_j.

For K>=c, an exact comparison with the previous uniform epsilon_K factor is

\[
\frac{|\gamma_K|}{(2K-c+3)\varepsilon_K}
=1+\frac{c-1}{(K+2)(2K-c+3)}<2.
\tag{6}
\]

All finite comparisons below use actual Z_j, not the coarse all-scale bound in the next section.

## 4. A stronger all-scale absolute-error bound

For the fixed schedule K=13m, the rational moments obey Z_j<=64^(-m). To prove this, it suffices to bound

\[
t(1-t)^2|2t-1|^{13}\le\frac1{64},\qquad0\le t\le1.
\]

For t<=1/2, the left side is at most t(1-2t)^13, whose maximum is (1/28)(13/14)^13<1/64. For t>=1/2, set v=1-t and bound by v^2(1-2v)^13<=1/225<1/64. Multiplication by t^j<=1 and integration on an interval of length one establish the claim.

Using |mu_j|<=(1+epsilon_K)M_j, M_j<=(4/27)^m, equation (6), and epsilon_K<1 gives

\[
|I_m(H)|\le8\varepsilon_K(4/27)^m64^{-m}.
\]

The preceding closed-period estimates give |alpha|>=|Delta|/2 and
|Delta|>=256^(-m)/(4m^2-1), while
epsilon_(13m)<=(4/3)(2^21/3^17)^m. Hence

\[
\boxed{\left|\pi-\frac pq\right|\le
\frac{256}{3}(4m^2-1)
\left(\frac{2^{25}}{3^{20}}\right)^m.}
\tag{7}
\]

The exponential base is approximately 0.009623. This improves the previous base 2^33/3^23 by the factor 27/256. The sharper rational quality certificate remains (5), retaining the actual q rather than substituting a large denominator envelope.

## 5. Telescoping the denominator combination itself

For any rational primitive identity with R(0)=0 (in particular for the odd R in the current construction),

\[
N(x)=(1+x^2)R'(x)-2(c-1)xR(x)
+\alpha(1+x^2)^{c-1},
\]

integration by parts using (1) proves

\[
\int_0^1N(x)S_K(x^2)\,dx
=2R(1)S_K(1)
-2\gamma_K\int_0^1xR(x)(2x^2-1)^K\,dx
+\alpha\int_0^1(1+x^2)^{c-1}S_K(x^2)\,dx.
\tag{8}
\]

For N=w[(6m+1)-(2m-1)x^2], the pi coefficient is zero by the closed base-period relation. Its primitive satisfies R(1)=2^(c-1)/binom(4m,2m). Therefore

\[
\boxed{(6m+1)\mu_0-(2m-1)\mu_1
=\frac{2^cS_K(1)}{\binom{4m}{2m}}
-2\gamma_K\int_0^1xR(x)(2x^2-1)^K\,dx.}
\tag{9}
\]

Thus the combination entering the reduced denominator is a rational boundary term plus one highest-order residual moment. There is no need to treat its K+1 Taylor summands as arithmetically unrelated. Equation (9) is an exact arithmetic handle for a redesigned truncation with controlled boundary values. By itself it does not guarantee an additional gcd; any proposed shared factor must still divide both actual integers.

## 6. A positive three-term recurrence for moment generation

The exact moments can also be generated without expanding S_K. Put

\[
W_{j,k}=\int_0^1x^{2m+2j}(1-x^2)^{2m}(2x^2-1)^k\,dx,
\qquad j=0,1.
\]

Integration by parts for the transformed beta weight gives

\[
\boxed{(2k+6m+2j+3)W_{j,k+1}
=-(2m-2j+1)W_{j,k}+2kW_{j,k-1}.}
\tag{10}
\]

Its initial data are W_(j,0)=M_j and
W_(j,1)=-(2m-2j+1)M_j/(6m+2j+3). Hence (-1)^k W_(j,k)>0 for every k. Define

\[
F_{j,k}=\left(\frac23\right)^c
\frac{(-1)^k\binom{c+k-1}{k}}{3^k}W_{j,k}>0.
\]

Then mu_j(K)=sum_(k=0)^K F_(j,k), and the increment recurrence is

\[
\boxed{9(k+1)(2k+6m+2j+3)F_{j,k+1}
=3(2m-2j+1)(c+k)F_{j,k}
+2(c+k)(c+k-1)F_{j,k-1}.}
\tag{11}
\]

All coefficients and all increments are positive. Uniform convergence of the original Taylor series on [0,1] shows that mu_j(K) increases to J_j as K tends to infinity at fixed m. This is a constructive recurrence for the two rational moments and makes positive-series acceleration available. It does not assert that their ratio is monotone: the individual centered-moment ratios can oscillate.

The same argument applied to the weight x^(2m)(1-x^2)^(2m+1) shows that every Taylor increment for mu_0-mu_1 is positive: its transformed beta-weight drift coefficient is -(2m+3). Hence

\[
\boxed{0<\mu_1(K)<\mu_0(K)\qquad(m\ge1,\ K\ge0).}
\tag{12}
\]

The closed base periods then imply the simpler all-K nondegeneracy statement

\[
\boxed{\alpha_H=\alpha_0\left(\mu_1-
\frac{6m+1}{2m-1}\mu_0\right)<0\quad\text{for every }K\ge0.}
\tag{13}
\]

Thus K=13m is not needed to secure a nonzero pi period. This frees the truncation length to be chosen by its arithmetic/error cost. It also removes the moment normalization from a convenient upper bound: using (12), equation (5) yields

\[
\boxed{Q\le\frac{4q^2|\gamma_K|}{
(2K-c+3)\alpha_0\bigl((6m+1)/(2m-1)-1\bigr)}(Z_0+Z_1).}
\tag{14}
\]

Here Z_j uses the even E=2 floor((K+1)/2), and 2K-c+3>0 is still required. No uniform epsilon_K bound is used in (14).

## 7. Ordinary convergence already holds at K=2m

For c=4m+1 and K=2m, the denominator 2K-c+3 equals two, so the profiled estimate applies even though K<c. Here E=2m.

For 0<=t<=1 the elementary bound

\[
t(1-t)^2(1-2t)^2\le49/864
\]

follows as follows. On [0,1/8], the derivative has the sign of 1-9t+10t^2>0, so the maximum there is 441/8192<49/864. On [1/8,1/2], use (1-t)^2<=49/64 and t(1-2t)^2<=2/27. On [1/2,1], setting v=1-t gives a bound v^2(1-2v)^2<=1/64<49/864. Thus Z_j<=(49/864)^m.

The weighted binomial theorem gives

\[
\binom{6m}{2m}(1/3)^{2m}(2/3)^{4m}\le1,
\qquad
|\gamma_{2m}|\le\frac23(6m+1)<4m+1.
\]

Also alpha_0((6m+1)/(2m-1)-1)>=2*16^(-m)/(2m-1). Substituting into (14) and dropping q^2 proves

\[
\boxed{\left|\pi-\frac pq\right|\le
2(4m+1)(2m-1)(49/54)^m\longrightarrow0.}
\tag{15}
\]

This lower truncation schedule spends fewer coefficient denominators. It is not being selected as a good-quality candidate merely from ordinary convergence; the exact q in (14) remains decisive. Its constructive role is to eliminate the conservative nonzero-period restriction on K, while leaving (5) and (14) available for cost-based choices of K.

## Verification

The companion standard-library script checks (1), (9), all closed-period fractions used, and (3) at 21 rational test points for each of the two schedules K=13m and K=2m, m=1,...,6. It checks the positive three-term increments, all prefix nonzero periods, and the K=2m coarse moment bound. It compares the two K=13m rational quality certificates exactly. The improved-to-old factors are approximately 63.99,777.61,10274.55,112821.13,1440635.23,15703397.30, respectively. These are improvements in upper certificates for the same generated fractions, not newly generated approximation records.

The script reads utilities from the prior pi-prescribed-cancellation checkpoint when run inside the repository. All universal claims (1)--(9) have the direct proofs above; rational-point tests are implementation audits only.
