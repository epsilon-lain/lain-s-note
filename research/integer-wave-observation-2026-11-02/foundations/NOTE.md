# Integer-wave foundations: stable topological readings and compatible refinement

Date: 2026-11-02. Status: a mathematically specified model and a finite-resolution theorem; not a physical claim, not a new construction of the number systems, and not a theorem on π's approximation quality.

## 1. Preserve the user's mechanism, improve the mathematical object

The earlier divisor kernel Δ_d(x)=d^-2 |sum_{r=0}^{d-1} exp(2π i r x/d)|² has exact divisor readings on integer x. Its exactness alone does not give a dynamical law, a finite-noise guarantee, or a new source of cancellations. Nor is a conjectural p-adic energy-transfer rule established by drawing interacting waves.

A more substantial next object is a continuous complex wave trajectory, together with an explicit observation operation. The observation should be invariant under many continuous deformations. This is stronger than merely writing an integer as an eigenvalue of an operator whose boundary conditions already contain integer periodicity.

Take a continuous closed path γ:[0,1]→C\{0}, with γ(0)=γ(1). Its normalized direction admits a continuous lifted angle θ:[0,1]→R satisfying γ(t)/|γ(t)|=exp(i θ(t)). Define

    W(γ) = [θ(1)-θ(0)]/(2π).

Closure forces W(γ)∈Z. This is the standard winding number, not a newly discovered fact. γ can be a highly irregular or modulated continuous wave, not a pure sinusoid. The integer is a stable topological reading of the whole path, rather than its instantaneous amplitude. Changing the integer requires crossing the excluded zero or breaking closure.

Underlying assumptions are explicit: a reference unit of time, continuity, a nonzero complex wave, and closure. Real/complex arithmetic is used to formulate the model. Thus this does not prove that arithmetic, or physical reality, originates from wave collapse.

## 2. Algebra of readings

For nonzero closed complex waves γ,η,

    W(γ η)=W(γ)+W(η),
    W(1/γ)=-W(γ),
    W(γ∘η)=W(γ)W(η)

in the last identity γ and η are interpreted as maps S¹→S¹ after normalization. The first two follow by adding/subtracting lifted phases; the third is multiplication of map degrees. The constant wave has reading 0, and one positively oriented loop has reading 1. These operations realize integer arithmetic through operations on continuous objects. They do not, by themselves, select a useful approximation to π.

A concrete phase-slip family is γ_λ(t)=exp(2π i k t)-λ with positive integer k and real λ≥0. W=k for λ<1 and W=0 for λ>1. At λ=1 the wave passes through zero k times; that is precisely where the discrete reading can change. For a general varying wave, the relevant new law must govern the appearance and transport of such zeros, not assert unexplained energy transfer from a picture.

## 3. Finite-resolution/noise theorem (proved here in explicit form)

Assume γ is closed and L-Lipschitz in the time variable:

    |γ(s)-γ(t)| ≤ L |s-t|,
    |γ(t)| ≥ m > 0.

Sample at t_j=j/N, j=0,…,N, with measured values y_j satisfying |y_j-γ(t_j)|≤η and y_N=y_0. Let P be the closed polygonal interpolation through y_j. Then

    ||P-γ||_∞ ≤ η + L/(2N).

If η+L/(2N)<m, P never meets zero and W(P)=W(γ) exactly.

Proof. For t=(j+s)/N, 0≤s≤1, let P₀ be interpolation through exact samples. By the triangle inequality,

    |P₀(t)-γ(t)|
    ≤(1-s)L s/N+sL(1-s)/N
    =2L s(1-s)/N≤L/(2N).

Interpolation of the measurement errors has norm≤η. If δ=η+L/(2N)<m, the straight-line homotopy H_u=(1-u)γ+uP satisfies |H_u(t)|≥m-δ>0. Winding is constant during this homotopy. QED.

This is a usable observation certificate: a continuous wave produces an exact integer despite finite sampling and bounded measurement noise, provided one knows a nonzero margin and a speed bound. The theorem is elementary; no claim of novelty is made.

The bound is not an aesthetic convention. Without a regularity/speed bound, finite samples cannot determine a winding: for any sample count N, all paths γ_ell(t)=exp(2π i ell N t) have the same samples γ_ell(j/N)=1 but winding ell N, for arbitrary integer ell. Thus wave→reading models must declare their information budget.

The winding can be computed from the polygon as the sum of principal arguments of y_{j+1}/y_j divided by 2π, when no edge passes through zero. In an exact implementation, orientation/crossing tests plus interval arithmetic can replace floating arguments.

## 3A. Data-driven upgrade: no unknown true margin is required

The true margin m may be unavailable. It can be eliminated from the certificate. Suppose only the Lipschitz bound L and sample-error bound η are known. Compute the measured polygon P and its exact distance d_P from the origin. Put δ=η+L/(2N). If

    d_P² > δ²,

then the unknown continuous wave γ is nonzero, and W(γ)=W(P). This is sufficient even though no lower bound on |γ| was supplied beforehand.

Proof. The interpolation estimate gives ||P-γ||∞≤δ. Every point of the straight homotopy P+u(γ-P), 0≤u≤1, has norm≥d_P-δ>0. In particular γ never vanishes and the winding is invariant. QED.

The test is exact for rational sample coordinates. For an edge a→b, let v=b-a. If v≠0, the nearest segment parameter is s=clamp(-a·v/(v·v),0,1), and the squared distance is |a+s v|². If v=0, it is |a|². Take the minimum over edges. All operations use rational arithmetic; no square root, angle integration or numerical zero search is needed.

### Application: certified root counts of auxiliary polynomials

Let A(z)=sum_{j=0}^d a_j z^j with rational real coefficients. Traverse the positively oriented square [-1,1]+i[-1,1] at constant speed over t∈[0,1]. Its speed is 8 and |z(t)|≤2, so γ(t)=A(z(t)) has the rigorous bound

    L = 8 sum_{j=1}^d j |a_j| 2^(j-1).

Contour samples at rational times have rational real and imaginary parts. Therefore the data-driven test above gives an exact certificate that A is nonzero on the entire contour. Its polygon winding then equals the number of zeros inside the square, counted with multiplicity, by the argument principle. This does not require knowing the roots or a precomputed boundary margin. Rational complex coefficients work as well if |a_j| is replaced by the rational upper bound |Re a_j|+|Im a_j|.

`verify_polynomial_contour.py` implements segment-distance projection, adaptive dyadic sampling, exact ray-crossing winding, and the argument-principle certificate with `Fraction` arithmetic. Its six examples (three polynomials, each with noiseless and bounded-noise samples) all pass:

| Polynomial | Noise bound | Certified sample count | Zeros inside |
| --- | ---: | ---: | ---: |
| z²+1/4 | 0 and 1/1000 | 32 | 2 |
| z²+4 | 0 and 1/1000 | 8 | 0 |
| (z-1/2)³ | 0 and 1/1000 | 1024 | 3 |

The polynomial z-1 has a zero on the contour and is correctly rejected rather than assigned an interior count. A failed margin test only means “not yet certified”; unless an exact zero sample is found in the noiseless case, it does not prove a boundary zero. Adaptive finite sampling terminates for a genuinely nonzero boundary in the noiseless case: d(P,0) tends to min|γ|>0 while δ tends to zero. With noise bound η, termination is assured if the true boundary margin is greater than 2η.

The same construction works on a translated or scaled rational square enclosing a particular near-root region; its path speed and radius bound must be updated. This is a concrete bridge back to the earlier auxiliary-polynomial near-root work. It supplies zero-count certification, while the separate rational-height/denominator problem remains open.

## 4. Rational and real readings from the same dynamical viewpoint

An integer reading counts net turns over a closed observation. A rational p/q can represent p turns after q reference-clock steps. For the rigid wave exp(2π i α t), integer-time recurrence occurs exactly when α is rational; if α=p/q is reduced, q is the smallest positive integer observation time giving closure. This is standard commensurability.

For genuinely variable local motion, use an increasing continuous map F:R→R with F(x+1)=F(x)+1. This is a lift of an orientation-preserving circle homeomorphism. Define its translation number

    ρ(F)=lim_{n→∞}(F^n(x)-x)/n.

The limit exists and is independent of x. Rational ρ=p/q is equivalent to a periodic orbit on the circle. Irrational ρ has no periodic orbit. These are standard Poincaré theory, not hypotheses about π.

There is a particularly clean compatible-observation theorem. Let

    m_n = floor(F^n(x)-x),
    J_N = intersection_{1≤n≤N} [m_n/n,(m_n+1)/n].

Then ρ(F)∈J_N, J_{N+1}⊂J_N, and length(J_N)≤1/N. Consequently the observations have exactly one compatible real limit. No continued-fraction digit is an input.

Proof of the interval assertion. Put G=F^n-m_n; it is increasing and commutes with unit translations. G(x)≥x implies G^k(x)≥x, hence nρ(F)-m_n≥0. Similarly F^n-(m_n+1) sends x below x, so nρ(F)-(m_n+1)≤0. The other assertions follow by intersection and the n=N interval. Existence and x-independence follow from the standard fact that displacement of each iterate varies by less than 1 over a unit interval, followed by sub-/superadditivity.

This supplies the previously missing consistency/refinement rule for a specified class of wave dynamics. Its universal guarantee is 1/N, not a hidden q^-2 improvement. The dynamics F must still be provided independently. Choosing F(x)=x+π puts the target in by hand and is only a diagnostic representation; it constructs no new information about π.

## 5. Complex extension records both phase and amplitude

Let X:[0,T]→C\{0} be an observed trajectory with an unwrapped continuous phase θ(t). Define the complex rate

    z_T = [θ(T)-θ(0)]/(2πT)
          - i [ln|X(T)|-ln|X(0)|]/(2πT).

For X(t)=X(0)exp(2π i z t), z=a+ib, this reading is exactly z: a is phase advance per unit time and b is negative logarithmic amplitude growth per 2π unit time. If only the endpoint multiplier X(T)/X(0) is retained, z_T is ambiguous by integer multiples of 1/T in its real part. Continuous path history or a winding certificate fixes that ambiguity. This is the logarithm/Floquet-exponent structure, not a mysterious new complex multiplication law.

This phase-plus-growth complexification must not be confused with a two-phase construction. For z=a+ib with b≠0, X(t)=exp(2π i zt) has nonconstant amplitude, so X(q)=X(0) is impossible for any positive integer q, even if both a and b are rational. It is false that Q+iQ produces finite closure in this one-complex-wave model.

A separate two-phase model is Y_z(t)=(exp(2π i a t),exp(2π i b t)) in S¹×S¹. Here integer-time closure exists exactly when a,b are both rational, using a common integer denominator. This construction represents C as two real phase rates and uses no exponential growth. Complex multiplication acts by the matrix [[a,-b],[b,a]] on two real coordinates; it is extra structure, not an automatic consequence of closure.

These models also induce different approximation tasks. Restricting to a common positive integer denominator for two independent real coordinates is simultaneous approximation, whose generic Dirichlet exponent is 3/2. Complex rational approximation p/q with Gaussian integers p,q∈Z[i], q≠0 instead has the natural critical scale |q|^-2. Neither exponent should silently be transferred to the other model. No complex approximation theorem is being claimed here.

The hierarchy is therefore: integer=stable closed-path charge; rational=commensurate charge/clock rate; real=compatible long-time phase rate. Complex numbers admit at least two different extensions—phase plus growth, or two phases—whose observations and closure rules must remain distinct.

## 6. What this contributes to the constructive π project

The immediate reusable tool is a finite-data winding/zero-count certificate. For analytic auxiliary functions, contour winding can certify root counts without locating every root. For continuous deformations, it detects the exact event at which the integer reading changes: a zero crosses the observation contour. This is relevant to the previous nearby-complex-root work and to choosing new analytic kernels.

It does not yet give a new rational approximant. The constructive next law must couple (a) a trajectory or analytic function defined independently of π's digits, (b) its observable integer charge/period, and (c) a quantitative small remainder with the fully reduced denominator accounted for. Merely replacing p,q by winding numbers leaves that task open.

A useful experimental module would evolve a finite integer-coefficient wave/kernel, certify its boundary margin and winding, and keep both the analytic residual and the coefficient/denominator budget. The other team's exact period Gram construction can supply the arithmetic component. Topological stability protects the reading; it cannot be silently substituted for a favorable approximation rate.

## 7. Validation and primary sources

`verify_rational_polygon.py` uses exact `fractions.Fraction` arithmetic only. It passes 80 noisy, margin-certified examples, an origin-crossing rejection, orientation reversal and additivity. The underlying closed wave traverses four rationally parametrized quarter circles: ((1-u²)/(1+u²),2u/(1+u²)), 0≤u≤1, rotated into each quadrant. With unit radius, each whole turn in unit time has speed at most 8. For k turns and horizontal shift c≥0, L=8|k| and m=|c-1| are rigorous rational bounds; the script checks η+L/(2N)<m exactly before reading a winding. Thus the margin theorem is exercised with exact inputs, not floating arguments. The numerical script is an additional independent consistency check.

`verify_winding.py` passes 120 deterministic noisy examples of shifted circular waves, including positive/negative/zero winding and both sides of a phase-slip threshold, plus 15 aliasing examples. The script uses floating arithmetic only as consistency checks; the theorem above is proved analytically.

- R. Melrose, *Lectures 3 and 4: The Circle*, MIT: circle spectral/Fourier framework. https://math.mit.edu/~rbm/18.158-S13/Lectures3-4.pdf . Retrieved reference `turn79view0` (parent must reopen before citing).
- R. Sjamaar, *Manifolds and Differential Forms*, §4, Corollary 4.14 and (4.13): integer winding as an integral of the angle form. https://pi.math.cornell.edu/~sjamaar/manifolds/manifold.pdf . Retrieved reference `turn79view1`.
- J. Milnor, *Dynamics Notes*, Chapter IV §14, Lemma 14.5, Theorem 14.10 and formula (14:6): translation numbers, rational periodicity, and compatible finite observations. https://www.math.stonybrook.edu/~jack/DYNOTES/dn14.pdf . Retrieved reference `turn81view0`.

Source-backed material is standard. The quantitative sampling proof and organization into an observation protocol are worked out here; no originality priority is asserted.
