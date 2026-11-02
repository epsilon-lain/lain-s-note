# A unique, simple, purely imaginary root at the Newton scale

Draft research result, 2026-11-02. This is a proof from the pinned repository's
period nonvanishing, endpoint moments and endpoint trace asymptotic. It does
not assert an arithmetic approximation theorem for pi.

## Objects and inherited input

Put t=2*pi, omega=it, A_j=product_{a=0}^{j-1} a!,
Delta_j=A_j^4/A_{2j}, with A_0=Delta_0=1. Let

H_n(Z)_{ij}=integral_0^1 u^(i+j) exp(Zu) du,

V_n(Z)_{ij}=Z^(-1) integral_0^infinity exp(-v)(1-v/Z)^(i+j) dv.

Let D_n(Z,eta)=det(H_n(Z)+eta V_n(Z)). The formal-period polynomial satisfies

P_n(Z)=Z^(n^2)/A_n^2 * D_n(Z,1-exp(Z)).

The established period estimates are

|det H_n(omega)|>Delta_n/6,

T_n=tr(H_n(omega)^(-1)V_n(omega))
   =i(-1)^n K_n chi_{n-1}(t)/chi_n(t) (1+O(n^(-1/2))),

K_n=[(2n-2)!]^2(2n-1)!/[(n-1)!]^4 t^(2n-1),

where 1/6<chi_j(t)<=1 for j>=2. In particular, eventually
K_n/12<=|T_n|<=12K_n. Also K_n grows faster than every fixed power of n.

## Mixed endpoint insertions

For 0<=k<=n define m=n-k and

B_{n,k}(t)=(Delta_m/Delta_n) A_k
             product_{j=0}^{k-1}(2m+j)! * t^(-k(2n-k)).

Empty products are 1, so B_{n,0}=1 and B_{n,1}=K_n.
If C_{n,k}(Z) is the coefficient of eta^k in
D_n(Z,eta)/det H_n(omega), then for Z!=0,

|C_{n,k}(Z)|
 <=6 exp(m max(Re Z,0)+k|Z|) B_{n,k}(|Z|).                 (1)

Proof: determinant multilinearity and the Vandermonde-square identity give
the mixed integral with factor 1/(k!m!), m variables u in [0,1] and k
endpoint variables x_a=1-v_a/Z, weighted exp(-v_a)/Z. The physical
Vandermonde integral is m!Delta_m. Bound each cross term by
|x_a-u_j|<=(v_a+|Z|)/|Z|. The endpoint Vandermonde contributes
product_{a<b}|v_a-v_b|^2/|Z|^(k(k-1)), and the k endpoint measures
contribute |Z|^(-k). Thus the total denominator exponent is
2mk+k(k-1)+k=k(2n-k). Shift y_a=v_a+|Z|. The Vandermonde is unchanged;
the exponential incurs exp(k|Z|); enlarging the y domain to [0,infinity)^k
is valid because the majorant is nonnegative. The Laguerre integral equals

1/k! integral exp(-sum y_a) product y_a^(2m)
                      product_{a<b}(y_b-y_a)^2 dy
 = A_k product_{j=0}^{k-1}(2m+j)!.

For completeness, the monic Rodrigues polynomial

L_j(y)=(-1)^j e^y y^(-2m) (d/dy)^j[e^(-y)y^(j+2m)]

is orthogonal to degrees below j by j integrations by parts; its squared
norm is j!(j+2m)!. The monic change of basis has determinant one, so the
factorial-moment Hankel determinant is the product of those norms. Divide
by |det H_n(omega)|>Delta_n/6. The formula includes m=0 and k=0.

For 0<=k<n, exact cancellation gives

B_{n,k+1}/B_{n,k}=k! K_{n-k}/(2n-2k)_k,                    (2)

where (a)_k=a(a+1)...(a+k-1), with (a)_0=1.
Since (2m)_k>=k! for m>=1, and eventually
max_{1<=m<n} K_m=K_{n-1},

B_{n,k}<=K_n K_{n-1}^{k-1}  (1<=k<=n).                   (3)

To justify the maximum, direct division gives

K_{m+1}/K_m=8(2m-1)^2(2m+1)/(m t^2).

This is >1 for m>=2, while K_n eventually exceeds K_1. Moreover

K_{n-1}/K_n=(n-1)t^2/[8(2n-3)^2(2n-1)]=O(n^(-2)).       (4)

## Uniform control as the period variable moves

On |Z-omega|<=1/n (n>=2), inequality (1) yields

|C_{n,k}(Z)|<=6e c^k B_{n,k}(t), c=exp(t+1).              (5)

Indeed m max(Re Z,0)<=1, |Z|<=t+1/n, and
log(t/|Z|)<=1/(n(t-1/n))<=1/(5n). Multiplication by
k(2n-k)<=2nk costs at most 2k/5 in the exponential; together with
k/n<=k/2 this is bounded by k. This estimate is uniform in k.

Let w=Z-omega and fix any constant L>0. Eventually |w|<=L/K_n is
contained in |w|<=1/(2n). For the k=0 term, the physical integral gives

|C_{n,0}(omega+w)-1|<=6n|w|exp(n|w|).                    (6)

For the k=1 term, the Cauchy derivative bound in the radius-1/n disk gives

|C_{n,1}(omega+w)-T_n|<=12e c n K_n |w|                  (7)

on its inner half. Also eta=1-exp(w) satisfies
|eta|<=2|w| and |eta+w|<=|w|^2 for all sufficiently small w.

Using (3),(5), uniformly on |w|<=L/K_n,

sum_{k=2}^n |C_{n,k}(omega+w)eta^k|
 <=6e c^2 K_n K_{n-1}|eta|^2/(1-cK_{n-1}|eta|)
 =O_L(K_{n-1}/K_n)=O_L(n^(-2)).                         (8)

Equations (6)--(8) prove the locally uniform normal form

D_n(omega+w,1-exp(w))/det H_n(omega)
       =1-wT_n+O_L(n^(-2))                              (9)

on |w|<=L/K_n. The additional errors O_L(n/K_n) are smaller than n^-2.

## Root theorem

Take L=48. On |w|=48/K_n the linear polynomial has
|1-wT_n|>=3; the error in (9) tends to zero. Rouché's theorem therefore
gives exactly one zero (counted with algebraic multiplicity) of D_n in this
disk, hence exactly one of P_n since the prefactor is nonzero there.
It is simple. Applying (9) at that zero w_n gives

w_n=1/T_n * (1+O(n^(-2))).                               (10)

The real-coefficient parity P_n(-Z)=(-1)^nP_n(Z) makes its zero set invariant
under Z -> -conjugate(Z). That map fixes omega and preserves the disk.
Uniqueness forces the root zeta_n=omega+w_n to be purely imaginary.

The inherited chain-rule identity says

P_n'(omega)/P_n(omega)=F_n'(omega)/F_n(omega)-T_n,

with F_n'/F_n=O(n^2). Consequently

zeta_n-omega=-P_n(omega)/P_n'(omega) * (1+O(n^(-2))).       (11)

In particular,

zeta_n-omega=-i(-1)^n chi_n(t)/(K_n chi_{n-1}(t))
                            * (1+O(n^(-1/2))),

and

log|zeta_n-omega|=-2n log n-[log(16/pi^2)-2]n+O(log n).    (12)

This unique disk root is the nearest root, since its distance is at most
24/K_n eventually whereas every other root lies outside 48/K_n.

For even n the root lies below 2pi i; for odd n it lies above it, for all
sufficiently large n.

## Independent endpoint variable

At fixed Z=omega, the same coefficient bounds give exactly one simple root
Y_n in |Y-1|<=48/K_n and

Y_n-1=-1/T_n*(1+O(n^-2)).

The reciprocal-conjugate symmetry maps roots by Y -> 1/conjugate Y. To use
uniqueness rigorously, first note |Y_n-1|<=24/K_n eventually, so its image
has distance |Y_n-1|/|Y_n|<=25/K_n eventually and remains in the unique disk.
Hence Y_n=1/conjugate Y_n and |Y_n|=1.

These axis and unit-circle results add analytic location, uniqueness and a
two-term logarithmic scale. They give no degree, height or rational-denominator
improvement. The real algebraic ordinate zeta_n/i is not generally rational.

## Finite checks

`verify_unique_root.py` passed 1377 exact B recurrences (including all k=n−1
last transitions), 81 exact k=n endpoint identities, 1377 geometric majorants,
49 Laguerre factorial determinants, and 27 mixed coefficient inequalities.
The mixed coefficient checks deliberately use synthetic physical measure du
and rational imaginary endpoint Z=6i, so every operation is Gaussian rational.
They audit the integral algebra and majorant, not the transcendental period
value. Optional 180-digit mpmath examples for n=4,6,8,10,12 find imaginary
roots with K_n times distance 1.477,1.318,1.211,1.157,1.125. These are numerical
illustrations, not certified root isolations or asymptotic proofs.
