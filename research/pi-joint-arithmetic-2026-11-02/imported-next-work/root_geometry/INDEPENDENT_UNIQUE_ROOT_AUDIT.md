# Independent audit of the unique near-period root argument

2026-11-02. This note audits `UNIQUE_NEAR_PERIOD_ROOT.md` by checking its all-dimension estimates rather than extrapolating finite root computations. It is an independent agent review, not an external peer review.

## Result of the audit

No gap was found. Subject to the cited inherited period nonvanishing and endpoint-trace asymptotic, the argument proves that, for all sufficiently large n, the disk centered at omega=2*pi*i of radius 48/K_n contains exactly one zero of P_n, counted with multiplicity. That zero is simple, purely imaginary, and the nearest zero to omega. Its displacement is

    zeta_n-omega = -P_n(omega)/P_n'(omega) * (1+O(n^(-2))).

It also proves exactly one simple root in the corresponding independent endpoint disk; that endpoint root lies on the unit circle. These two roots belong to different specializations of the two-variable determinant and are not identified with each other.

## Checks that matter for the proof

1. **Mixed determinant normalization.** Choosing k endpoint functionals and m=n-k actual interval functionals gives the prefactor 1/(k!m!), because binomial(n,k)/n!=1/(k!m!). Integrating the interval Vandermonde majorant gives m!*Delta_m. The endpoint Vandermonde contributes |Z|^(-k(k-1)), the cross terms contribute |Z|^(-2mk), and the endpoint measures contribute |Z|^(-k). The combined power is exactly k(2n-k).

2. **Laguerre integral.** After shifting v to y=v+|Z| and enlarging the nonnegative majorant domain, the remaining integral divided by k! is the factorial-moment Hankel determinant. Its value is A_k times the product of (2m+j)! for j=0,...,k-1. There is no missing k! or additional factorial from the shift.

3. **All-k recurrence and edges.** Direct factorial division gives

       B_(n,k+1)/B_(n,k)=k! K_(n-k)/(2n-2k)_k.

   This includes k=0 and k=n-1. In particular B_(n,1)=K_n and B_(n,n)=1/U_n. K_m is not increasing from m=1: K_2<K_1 at t=2*pi. The draft correctly uses the eventual maximum instead; K_m increases for m>=2 and K_3>K_1, so the maximum among m<n is K_(n-1) for n>=4.

4. **Uniformity in the moving period variable.** On the radius-1/n disk, the factors arising from |Z| versus t are bounded uniformly in k by exp(k). The endpoint coefficients are analytic there. The radius-1/(2n) Cauchy circle around any point in the inner half-disk gives the derivative bound 12*e*c*n*K_n. This controls variation of the linear coefficient without assuming it is constant in Z.

5. **Nonlinear tail.** On |w|<=L/K_n the endpoint perturbation eta=1-exp(w) has magnitude at most 2|w|. The ratio c*K_(n-1)*|eta| tends to zero, and the complete k>=2 tail is O_L(K_(n-1)/K_n)=O_L(n^(-2)). The omitted linear-coefficient and actual-integral errors are O_L(n/K_n), which are smaller. No unproved termwise asymptotic in k is used.

6. **Rouche and location.** The linear comparison 1-w*T_n has boundary magnitude at least 3 on |w|=48/K_n once |T_n|>=K_n/12. The remainder tends to zero uniformly. The factor Z^(n^2)/A_n^2 has no zero in that disk. Exactly one zero counted with multiplicity implies simplicity. The reflection Z->-conjugate(Z) preserves both the polynomial zero set and this disk, forcing its unique zero onto the imaginary axis.

7. **Nearest root and endpoint circle.** The zero displacement is at most 24/K_n eventually, while all other zeros are outside 48/K_n. For the independent endpoint root, reciprocal conjugation maps a root at distance at most 24/K_n from 1 into distance at most 25/K_n eventually; uniqueness therefore forces modulus one. This extra inner-distance step is necessary because reciprocal conjugation does not preserve an arbitrary disk about 1 exactly.

## Scope

The logarithmic-distance refinement and exact imaginary-axis location are analytic advances for the same family. The proof supplies no bound for the irreducible degree or height of the selected algebraic ordinate and no rational approximation with a controlled denominator. It does not establish bad approximability, its negation, or a new irrationality exponent for pi.
