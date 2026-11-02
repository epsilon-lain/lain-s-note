# Critical algebraic endpoints and residue-normalized trace

This is a constructive fixed-number-field variant of the rational endpoint design. It was suggested by the parent investigator to avoid paying for increasingly large rational Pell coordinates. The exact construction works, but the tested quartic-kernel variants do not yet improve approximation quality.

Let K=Q(sqrt2), v=sqrt2−1, and

\[
R(x)=\frac{((x-1)^2+v^2)^{2a}(x-1)^{2b}}
{x^{c+1}(2-x)^{c+1}},\qquad a\ge1,b,c\ge0.
\]

Integrate from 1−iv to 1+iv along the vertical segment and divide by i. The partial-fraction symmetry x↦2−x gives

\[
I_1=\frac\pi2\rho+B,
\qquad \rho,B\in K,
\]

where rho is the simple-pole coefficient. Let sigma be the nontrivial field embedding. Its signed endpoint parameter is v'=−sqrt2−1; keeping that sign is essential. The conjugate contour, from 1−iv' to 1+iv', is oriented downward. Continuous Log x−Log(2−x) along it changes by −3iπ/2, hence

\[
I_2=-\frac{3\pi}{2}\sigma(\rho)+\sigma(B).
\]

If rho≠0,

\[
\boxed{\frac{I_1}\rho+\frac{I_2}{\sigma(\rho)}
=-\pi+\operatorname{Tr}_{K/\mathbb Q}(B/\rho).}
\]

Thus the rational candidate is the trace of B/rho. It is derived from the chosen kernel; no target convergent is supplied.

The arithmetic charge is explicit. Multiply the entire coefficient pair rho,B by one common positive integer D, and write

\[
D\rho=u+w\sqrt2,\quad DB=r+s\sqrt2,
\qquad u,w,r,s\in\mathbb Z.
\]

Then

\[
N=u^2-2w^2,\quad M=2(ru-2sw),\qquad
\operatorname{Tr}(B/\rho)=M/N.
\]

For g=gcd(|N|,|M|), the reduced denominator is |N|/g and

\[
\boxed{Q=\frac{|N|}{g^2}|N\pi-M|.}
\]

Taking the trace avoids an increasing rational Pell endpoint denominator, but residue normalization introduces the norm N. Its size and the two conjugate integral errors must both be controlled. A construction with small I1 alone is insufficient. This is an arithmetic accounting identity, not a success theorem.

## Executed tests

`quadratic_trace_contour.py` computes all coefficients exactly in K, then evaluates the rational candidates using exact Machin alternating-series error intervals. Among a≤4,b≤4,c≤10,a+b+c≤14, 114 positive rational candidates were evaluable. No good approximation emerged; the smallest observed quality in that finite set was the poor candidate 11, with Q≈7.8584. This experiment is not a reason to promote the unmodified quartic trace kernel as a promising candidate.

`verify_quadratic_trace.py` independently integrates both signed conjugate contours at 100 decimal digits for six parameter choices and checks the two different pi coefficients and the normalized trace identity. These numerical checks supplement the exact field derivation.

The parent investigator is separately testing simpler arctangent Padé trace constructions. The exact norm formula above is available to evaluate that broader idea; the negative quartic control should not be mistaken for a theorem excluding all algebraic-trace designs.

## Coupled modification actually executed

To control both embeddings together, the final experiment replaced the square of one endpoint factor by independent powers

\[
R(x)=\frac{Q_-(x)^aQ_+(x)^b(x-1)^{2d}}
{x^{c+1}(2-x)^{c+1}},
\qquad Q_\pm(x)=(x-1)^2+(\sqrt2\pm1)^2.
\]

Here a,b≥0,a+b>0,d≥0,c≥1. Under sigma, the two factors swap. On the wide conjugate contour the small factor has interior zeros, so choosing its power can affect the conjugate cancellation rather than only scale the original integral. Properness is precisely c≥a+b+d; polynomial parts otherwise remain explicit. The residue-normalized trace and the integer norm-cost formula above remain unchanged.

`coupled_trace.py` scanned a,b≤4,d≤2,c≤8: 387 positive evaluable candidates. The best over this finite set was 3, with Q≈0.14159265; it occurs at (a,b,d,c)=(0,2,1,2), when the large endpoint factor is used alone. At the same pole exponent, the original small-factor square with d=1 gives 11. This improves that poor baseline, but is not a new good pi approximation. Requiring both a,b>0 is recorded separately in the JSON so the one-factor winner is not misrepresented as evidence for simultaneous-factor success.

`verify_coupled_trace.py` independently checks six parameter choices, using twelve 100-decimal signed-conjugate contour quadratures. Both different pi branches and the final rational trace identity pass. Zero-residue cases are excluded before attempting normalization.

This final coupled experiment does not yield a quality-zero sequence or a promising new record. A next orthogonality-based choice of the two factor weights would need to improve the actual norm-cost formula, not only the small embedding's integral.
