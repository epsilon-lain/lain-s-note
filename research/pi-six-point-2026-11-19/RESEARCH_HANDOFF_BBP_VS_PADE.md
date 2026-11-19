# Research handoff: do not duplicate the Padé line; connect the BBP sparse-exception result

Date: 2026-11-19  
Purpose: consolidate the newly supplied Yohaku relay with the existing Padé work and set the next research priority.

## 1. The newly supplied result is a different construction

The report [SPARSE-EXCEPTION-RELAY-RESULTS.zh.md](https://github.com/epsilon-lain/yohaku/blob/aa7188a8ba945cfbb31b047111fbf3ee0668e9f5/research/pi-equivalence-native-2026-11-18/SPARSE-EXCEPTION-RELAY-RESULTS.zh.md) studies the specific BBP–finite-integral rational family \(X_{N,M}\), its reduced denominator \(Q_{N,M}\), and
\[
Y_{N,M}=16^N X_{N,M}.
\]
Its main mechanism is 2-adic: a 12-term exact finite-difference formula, a terminal sum modulo \(2^r\), and a sparse-exception argument. It is **not** a Padé construction and does not repeat the two-arctangent Taylor-tail determinant.

The reported results are:

- A's finite-difference/terminal-congruence lemma reduces \(v_2(Y_{N,4N})\) to nonvanishing of a specific terminal sum \(U_r(N)\) modulo \(2^r\), for even \(N\ge4\).
- B proves an infinite sparse subfamily: if \(N=2^s+14\), \(s\ge13\), and \(M=4N\), then \(v_2(Y_{N,4N})=0\), hence \(v_2(Q_{N,4N})=4N\).
- The report explicitly leaves open the terminal nonvanishing for all even \(N\), the full class \(M=4N+o(N)\), and all rational approximants.

The linked follow-up explains why this is a meaningful but incomplete bridge: the two-dimensional spacing theorem makes linearly large 2-adic losses sparse, but \(O(\log B)\) exceptional points up to scale \(B\) can still be infinite. The global problem is to rule out or control those isolated exceptions on the same parameter point.

## 2. Padé work already exists; do not restart it

There are at least three distinct earlier lines in the notes:

### 2.1 Niven-integral Padé compensation

See [pi-niven-pade-compensation-2026-11-18/README.md](../pi-niven-pade-compensation-2026-11-18/README.md).

Already proved there:

- For the critical Niven-integral scale, diagonal \([m/m]\) Padé with \(1\le m\le\sqrt n\) cannot make the absolute Padé remainder smaller than the relevant rational lattice spacing after accounting for its reduced denominator.
- The joint endpoint-displacement/Padé-remainder comparison rules out cancellation for all nonconstant rectangular Padé choices when the threshold is chosen at the actual error scale.
- The \([0/b]\) inverse-series construction can genuinely remove the original denominator's 2-adic factors, including all powers of 2 when \(q\) is a power of 2; however, the joint estimate still shows that at the natural threshold its remainder is smaller than the endpoint displacement, so their sum does not become small.

Therefore do not spend another cycle re-deriving Padé denominator cancellation or testing generic rectangular Padé degrees. A genuinely new Padé proposal must change the endpoint functional or prove a new joint estimate not covered by that note.

### 2.2 Exponential-function diagonal Padé detector

See [pi-critical-wave-2026-11-19/PADE.md](../pi-critical-wave-2026-11-19/PADE.md).

This is the separate \(e^{i\pi}=-1\) detector. Its exact rational denominator, positive periodic remainder, and derivative scale are known. The current note explains why separately forcing both the periodic remainder and the displacement below the \(q^{-n}\) grid is much too expensive; the \(p/q>\pi\) case still has a possible sign cancellation not resolved by separate bounds.

### 2.3 OpenAI interpolation determinant

See [openai-pi-exact-determinant-2026-11-18/README.md](../openai-pi-exact-determinant-2026-11-18/README.md).

This is a different determinant route based on the OpenAI interpolation proof and a fixed nonzero subminor. It has already improved the arithmetic denominator and analytic estimates for the same minor, but the two bounds have not been shown to cross.

These three lines are not the same as the newly supplied BBP finite-integral construction. Keep their assumptions and target objects separate.

## 3. Relation to the current six-point tail-column line

The six-point construction based on
\[
\arctan(2/5)+\arctan(3/7)=\pi/4
\]
uses the signed linear forms
\[
f_{m,n}=\pi/4-T_m(2/5)-T_n(3/7).
\]
Its exact tail-column determinant is a legitimate new interface, but the raw Taylor denominator budget is not currently enough. This is a separate obstacle from the Padé compensation theorem above.

Do not silently substitute a Padé denominator estimate from the Niven-integral problem into the two-arctangent tail determinant: the rational functions, endpoint, denominators, and error terms are different. Any proposed transfer must be proved for the exact matrix and exact rational entries at hand.

## 4. Recommended next work

Priority should move to the new BBP relay, not to another generic Padé experiment:

1. Independently reconstruct the exact rational \(Y_{N,4N}\) or the terminal residue for a set of even \(N\), and verify the A–B terminal-sum identity against the direct fraction.
2. Attack the all-even nonvanishing problem for \(U_r(N)\), beginning with the known exceptional sample \(N=110\) where the terminal residue has valuation 14 and the direct \(Y\) has valuation 4. Any proposed universal low-bit rule must survive this counterexample.
3. Prove or refute a larger explicit residue class than \(N=2^s+14\). Keep the exact seed-shift hypotheses visible: the present result is a fixed-seed theorem, not a theorem for all even \(N\).
4. Only after the slope-4 case is controlled, extend to \(M=4N+b\) for fixed \(b\), then to \(M=4N+o(N)\). The relay's own task sheet makes this quantifier escalation explicit.
5. Connect the BBP-family exclusion to the global problem only through a separate theorem. Excluding this one approximation family does not prove \(\inf_q q\|q\pi\|>0\) over all rational denominators.

The first research artifact should be a short independent verification note with the exact hypotheses and counterexamples, not another Padé manuscript.

## 5. Scope discipline

- Do not call the BBP relay a Padé result.
- Do not call a sparse exceptional set finite merely because its count is \(O(\log B)\).
- Do not promote \(N=2^s+14\) to all even \(N\).
- Do not promote \(M=4N\) to \(M=4N+o(N)\) without a uniform error proof.
- Do not promote exclusion of the BBP family to bad approximability of \(\pi\).

This handoff is a routing and non-duplication record, not a new proof of the global conjecture.
