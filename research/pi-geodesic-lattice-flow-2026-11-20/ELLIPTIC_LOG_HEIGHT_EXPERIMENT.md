# Experiment: elliptic logarithm recurrence versus high-order formal cancellation height

Date: 2026-11-20  
Status: **exact rational experiment on one explicit elliptic curve, with high-precision real uniformization diagnostics**. This tests the height bottleneck; it is not a theorem about pi.

## 1. Curve and non-torsion point

Take
\[
E:\ y^2=x^3-2,\qquad P=(3,5).
\]
The discriminant is
\[
\Delta_E=-16\cdot27\cdot(-2)^2=-1728.
\]
By Lutz–Nagell, a rational torsion point on this integral short Weierstrass model has integer coordinates and either \(y=0\) or \(y^2\mid\Delta_E\). Here \(y^2=25\nmid1728\), so \(P\) is non-torsion.

Use the standard formal parameter at the identity
\[
\xi=-x/y.
\]

## 2. Real elliptic logarithm and an explicit recurrence

For the Weierstrass uniformization \(x=\wp(u)\), \(y=\wp'(u)/2\), the invariant differential is
\[
du=\frac{dx}{2y}.
\]
Let \(\Omega\) be the positive real period and \(u(P)\) the real elliptic logarithm, with sign chosen so that the distance-to-period calculation is positive. Numerically,
\[
\Omega=2.163681749013751009472202825697\ldots,
\]
\[
u(P)=0.580500474815398977631230884312\ldots,
\]
and
\[
\theta=\frac{u(P)}{\Omega}
=0.268292910951438486839781900129\ldots.
\]
At \(s=41\),
\[
41\theta=11.0000093490089779604310579053\ldots.
\]
Thus
\[
|41u(P)-11\Omega|
=2.0228280096978686235174119362\cdot10^{-5}.
\]
Exact rational elliptic arithmetic gives
\[
\boxed{
|\xi([41]P)|
=2.0228280096978686235174119479\cdot10^{-5}.
}
\]
The ratio of these two numerical quantities is \(0.9999999999999999999999999942\ldots\), as expected from the local expansion of the uniformization parameter near the identity.

This is a concrete instance of the linear-form-in-elliptic-logarithms problem: the small quantity is \(41u(P)-11\Omega\). David–Hirata-Kohno prove effective lower bounds for such linear forms, with a simultaneous version depending on the heights of the curve and points and the field degree; their abstract notes that the numerical constant is not fully explicit in a practically sharp form. The present experiment does not instantiate their theorem's constant.

## 3. Exact rational construction of B_m([s]P)

For each \(s\), define
\[
B_m([s]P)=
\sum_{j=2}^{m+2}
(-1)^j\binom{m+1}{j-1}
\frac{\xi([js]P)}{\xi([s]P)}.
\]
Every point multiple and every B_m was computed using exact rational arithmetic on the curve. The decimal logarithms below are diagnostics of those exact rational outputs.

### Three recurrence indices

| \(s\) | \(|\xi([s]P)|\) |
|---:|---:|
| 15 | \(5.2780127154\times10^{-2}\) |
| 26 | \(5.2759898871\times10^{-2}\) |
| 41 | \(2.0228280097\times10^{-5}\) |

Among \(1\le s\le100\), \(s=41\) gives the smallest observed \(|\xi([s]P)|\). This is a finite search, not a claim about the best recurrence.

### High-order cancellation at \(s=41\)

Here \(h(B_m)=\log\max(|\operatorname{num}B_m|,\operatorname{den}B_m)\) is the logarithmic rational height after reduction.

| \(m\) | \(\log_{10}|B_m-1|\) | \(h(B_m)\) | numerator/denominator bit length |
|---:|---:|---:|---:|
| 1 | \(-24.9452\) | \(40,835.5\) | 58,914 |
| 2 | \(-24.2223\) | \(95,282.8\) | 137,465 |
| 3 | \(-23.8298\) | \(176,953.8\) | 255,291 |
| 4 | \(-23.7271\) | \(285,846.9\) | 412,390 |
| 5 | \(-23.9267\) | \(449,188.8\) | 648,043 |
| 6 | \(-24.5288\) | \(666,978.2\) | 962,247 |

The weight l1 norms are \(3,7,15,31,63,127\), respectively.

**Observation:** for this fixed V=[41]P, increasing m from 1 to 6 raises the rational height by a factor of about 16, while \(|B_m-1|\) stays around \(10^{-24}\) and is not monotone decreasing. The formal estimate \(B_m=1+O_m(\xi(V)^m)\) is valid for each fixed m, but its implied constant depends on m; this one-point experiment shows why that dependence cannot be ignored.

For s=15 and s=26, the same exact experiment with m<=4 also shows \(|B_m-1|\) increasing as m grows. The high-order identity is algebraically correct, but it does not automatically provide a useful height/error tradeoff.

## 4. Canonical-height scale

The naive logarithmic heights satisfy
\[
\frac{h(x([41]P))}{2\cdot41^2}
\approx0.6747884178399647,
\]
\[
\frac{h(x([82]P))}{2\cdot82^2}
\approx0.6747884178399647.
\]
These are finite approximations to \(\widehat h_E(P)\), not a substitute for a full canonical-height computation.

The output heights h(B_m) are already enormous at m<=6. This is the quantitative bottleneck in our formal-group strategy.

## 5. What this tells us about David–Hirata-Kohno

The elliptic-logarithm theorem controls how small a nonzero form
\[
s\,u(P)-r\,\Omega
\]
can be in terms of arithmetic heights and coefficients. It explains why the recurrence at s=41 is exceptional but cannot be arbitrarily strong without arithmetic cost.

However, the theorem bounds the input small parameter \(\xi([s]P)\); it does not directly bound the height of the rational combination B_m([s]P). Our exact experiment isolates this missing step:
\[
\text{small elliptic logarithm}
\quad\not\Rightarrow\quad
\text{small-height high-order rational combination}.
\]

A plausible next improvement is to replace the binomial finite-difference weights by a redundant integer-weight system optimized for small coefficient norm, while preserving the moment conditions. Any such optimization must be proved to retain the p-adic congruences.

## 6. Reproducibility and limits

- elliptic_formal_height_experiment.py performs the elliptic group law and all B_m calculations with exact Fraction arithmetic.
- elliptic_formal_height_experiment.json stores the output diagnostics.
- verify_elliptic_high_order_weights.py independently verifies the integer interpolation identities for m=1,...,30.

The elliptic-period and logarithm values are high-precision numerical integrations, not interval-certified. The rational point arithmetic and the reduced rational B_m values are exact. This experiment is about the height bottleneck in our elliptic construction; it does not prove or disprove bad approximability of pi.
