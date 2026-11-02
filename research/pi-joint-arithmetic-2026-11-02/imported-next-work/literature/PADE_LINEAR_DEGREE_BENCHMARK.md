# Classical exponential Padé family: a linear-degree benchmark for the near-period root mechanism

Internal research note, 2026-11-02. The polynomials below are classical diagonal exponential Padé polynomials, not newly discovered objects and not asserted to be factors of the repository's P_n. The deductions here give an explicit comparison with its factorial near-root mechanism.

## 1. Exact normalization and remainder

For an integer m >= 1 set

    a_(m,k) = (2m-k)! m! / ((2m)! k! (m-k)!),  0 <= k <= m,
    p_m(Z) = sum a_(m,k) Z^k,
    q_m(Z) = p_m(-Z),
    g_m(Z) = p_m(Z)-q_m(Z).

The standard diagonal exponential Padé remainder is

    q_m(Z) exp(Z)-p_m(Z)
      = (-1)^m Z^(2m+1)/(2m)! * integral_0^1 t^m(1-t)^m exp(tZ) dt.       (1)

A source for the corresponding explicit Padé polynomials and beta-integral remainder (under a monic alternative normalization) is Driver–Temme, Quaestiones Mathematicae 22 (1999), formulas (1.2)–(1.5):
https://ir.cwi.nl/pub/1229/1229D.pdf

Equation (1) can also be proved without importing asymptotics: integrate the beta integral by parts 2m+1 times, using vanishing endpoint derivatives through order m−1. The remaining endpoint derivatives produce the displayed coefficients.

For ODD m >= 3 define the classical integer polynomial

    C_m(Z) = sum_(1 <= k <= m, k odd) (2m-k)!/[k!(m-k)!] * Z^(k-1)
           = (2m)!/(2 m! Z) * g_m(Z).                                  (2)

Every coefficient is integer since it is binom(2m-k,m) * (m!/k!). C_m is monic, even, and has degree m−1. Its nonintegral imaginary rational values have full degree denominator: for coprime integers a,b with b >= 2,

    den_reduced C_m(ia/b) = b^(m-1),   |C_m(ia/b)| >= b^(-(m-1)).         (3)

Indeed multiplying by b^(m−1), reduction modulo each prime divisor of b retains precisely the unit leading term i^(m−1)a^(m−1).

## 2. Two-sided actual-period value, with the linear term

Put omega=2 pi i and T=2 pi. Let X_m have the beta(m+1,m+1) density t^m(1−t)^m/B(m+1,m+1), and set

    chi_m = E cos(T(X_m−1/2)).

By symmetry, E exp(omega X_m) = −chi_m. The elementary inequality cos u >= 1−u^2/2 and Var(X_m)=1/[4(2m+3)] give

    0 < 1−pi^2/[2(2m+3)] <= chi_m <= 1,  m>=1.                         (4)

At omega, exp(omega)=1. Equations (1)–(2), together with B(m+1,m+1)=(m!)^2/(2m+1)!, give for odd m

    C_m(omega) = T^(2m) m!/[2(2m+1)!] * chi_m > 0.                    (5)

Consequently

    log |C_m(omega)|
      = −m log m + [1+2 log pi]m + O(log m).                            (6)

This is the arithmetic smallness of the monic integer polynomial. Its root error has twice the leading logarithmic decay, as shown next. Those quantities must not be conflated.

## 3. A unique simple imaginary root with the same factorial root scale

For each fixed radius R, coefficient comparison gives

    |p_m(Z)−exp(Z/2)| <= R^2 exp(R/2)/(8m),  |Z|<=R,
    |g_m(Z)−2sinh(Z/2)| <= R^2 exp(R/2)/(4m),
    |g_m'(Z)−cosh(Z/2)| <= exp(R/2)(R/2+R^2/8)/m.                      (7)

To prove the first bound, extend a_(m,k)=0 for k>m and use

    a_(m,k) = 1/(2^k k!) * product_(j=0)^(k-1) [1−j/(2m−j)]

for k<=m. The product deficit is at most k(k−1)/(2m); this bound also controls the omitted coefficients k>m. Sum the resulting exponential majorant and differentiate it to obtain (7).

The limiting function 2sinh(Z/2) has precisely one zero, simple, in |Z−omega|<1. Its boundary modulus is at least 2−cosh(1/2)>0. Rouché plus (7) therefore yields exactly one zero zeta_m of g_m, counted with multiplicity, in that disk, for sufficiently large m. Thus the zero is simple. The disk is invariant under Z -> −conjugate(Z), and the real odd polynomial g_m has this root symmetry. Uniqueness forces zeta_m to be purely imaginary. It is also the unique root of C_m in that disk.

For an explicit conservative threshold one can take every odd m>=10001: in |Z−omega|<=1 (hence |Z|<8), (7) implies |g_m'(Z)+1|<1/4. The remainder makes |g_m(omega)|<1/2. Rouché on the smaller disk |Z−omega|=2|g_m(omega)| then gives a unique root and

    (4/5)|g_m(omega)| <= |zeta_m−omega| <= (4/3)|g_m(omega)|.            (8)

Here (8) follows directly by integrating g_m' along the segment from omega to zeta_m. All other roots have distance >=1, so this root is nearest to omega.

In fact (7) and the small root distance refine this to

    zeta_m−omega = g_m(omega)[1+O(1/m)]
                 = i T^(2m+1)(m!)^2/[(2m)!(2m+1)!] * [1+O(1/m)],      (9)

where the last formula is for odd m. Thus zeta_m lies above omega on the imaginary axis for sufficiently large odd m, and

    log |zeta_m−omega|
      = −2m log m + [2+2 log(pi/2)]m + O(log m).                        (10)

The near-period root scale from the large determinant therefore already occurs in a classical family of degree m−1, with a unique simple imaginary root. This is a benchmark, not a contradiction: different constructions can produce the same leading root rate at very different degrees.

## 4. Exact height and why linear degree still does not yield q^(−2)

The positive coefficients c_k=(2m−k)!/[k!(m−k)!] for odd k decrease:

    c_(k+2)/c_k
      = (m−k)(m−k−1)/[(k+1)(k+2)(2m−k)(2m−k−1)] < 1.

Thus the naive coefficient height is exactly

    H(C_m) = (2m−1)!/(m−1)!,
    log H(C_m) = m log m + (2 log 2−1)m + O(log m).                     (11)

From (6), the baseline-value criterion |C_m(omega)| <= 1/[2 b^(m−1)] has asymptotic threshold

    m/b = e pi^2.                                                       (12)

More precisely: along m/b -> c>0, odd m -> infinity, its log left side after multiplying by b^(m−1), divided by m, tends log(e pi^2/c). Hence any fixed c>e pi^2 suffices eventually; c<e pi^2 fails eventually. No impossibility claim for other uses of the family is intended.

A direct derivative majorant from (11) on |Z|<=8 is

    |C_m'(Z)| <= m^2 H(C_m) 8^(m−1).

For p/q reduced, reduce 2p/q=a/b, b=q/gcd(2,q). Under the baseline criterion and |ia/b−omega|<=1, equations (3) and the mean-value integral yield

    |pi−p/q| >= 1/[4 b^(m−1) m^2 H(C_m) 8^(m−1)].                     (13)

The non-near case is automatic. Taking m~c b makes the logarithm of this bound

    −2c b log b + O_c(b),

still far weaker than a fixed positive constant times q^(−2). The benchmark shows that reducing the auxiliary degree from quadratic to linear is a useful simplification, but cannot alone bridge the arithmetic gap. The factorial coefficient height consumes the benefit.

## 5. Best matching primary literature and how to use it

1. Driver–Temme (1999), Zero and pole distribution of diagonal Padé approximants to the exponential function.
   https://ir.cwi.nl/pub/1229/1229D.pdf
   Exact coefficients, beta-integral remainder, Bessel/Hankel identities and root asymptotics. Formulas (1.2)–(1.5) directly support the control family above. Their P_n notation differs from the repository's determinant P_n.

2. Cullinan–Scheel (2020 preprint; 2022 publication), On the arithmetic of Padé approximants to the exponential function.
   https://arxiv.org/abs/2007.01329
   https://arxiv.org/pdf/2007.01329
   Explicit integer-normalized polynomials P(u,v,x)=sum_(j<=u) (u+v−j)!/v! * binom(u,j)x^j, represented by generalized Laguerre polynomials. Diagonal P(m,m,x) has full symmetric Galois group and is irreducible, by earlier Bessel-polynomial work. This irreducibility does NOT transfer automatically to the odd part divided by Z, nor to the determinant reflection factors F_±; use Newton polygons as a tool rather than cite a false identification.

3. Celsus–Deaño–Huybrechs–Iserles, The kissing polynomials and their Hankel determinants, arXiv:1504.07297 (2015 preprint; 2022 journal article).
   https://arxiv.org/abs/1504.07297
   https://doi.org/10.1093/imatrm/tnab005
   Exact oscillatory weight exp(iw x) on [-1,1]; Hankel determinants, complex zeros and a sigma-Painlevé V equation. The repository's genuine exponential moment determinant maps to this after u=(x+1)/2, w=−iZ/2. With h_(n−1)(w)=det integral_(-1)^1 x^(i+j)exp(iw x)dx,

       det integral_0^1 u^(i+j)exp(Zu)du = exp(nZ/2) 2^(−n^2) h_(n−1)(−iZ/2).

   However, globally this is the actual-endpoint determinant G_n(Z,exp Z), not the formal-endpoint polynomial G_n(Z,1). At Z=2pi i values agree but derivatives differ by ∂_Y G. Importing their zeros to P_n without this distinction would be incorrect. The paper is an exact mechanism match, not directly an arithmetic-cost theorem.

4. Basor–Chen–Ehrhardt (2009), Painlevé V and time dependent Jacobi polynomials.
   https://arxiv.org/abs/0905.2620
   https://doi.org/10.1088/1751-8113/43/1/015204
   Exponentially deformed Jacobi weight, Toda flow and Painlevé V. It supplies recurrence/differential structure for actual moments. Any formal-endpoint transfer must use D=∂_Z+Y∂_Y along Y=exp Z, not ∂_Z alone.

5. Zeilberger–Zudilin (2019 preprint; 2020 publication), The irrationality measure of pi is at most 7.103205334137...
   https://arxiv.org/abs/1912.06345
   https://sites.math.rutgers.edu/~zeilberg/mamarim/mamarimPDF/pimeas.pdf
   A proper arithmetic comparison: integral smallness is paired with integer linear forms and denominator reduction. It is a benchmark for the missing arithmetic output, not a statement about current best bounds. Search also found Sep2026 Bai arXiv:2609.11276 (7.101862832357) and the OpenAI Sept24 2026 repository claim mu(pi)=2; neither establishes bad approximability, which needs a fixed q^(−2) constant. The latter original proof has not been independently certified in this note.

## 6. Optional stronger classical root geometry, valid at every odd m>=3

All roots of C_m are purely imaginary and simple, without waiting for large m. This follows from two established classical facts, both stated and proved/used in Csordas–Charalambides–Waleffe, A new property of a class of Jacobi polynomials (Proc. AMS 133 (2005), 3551–3560): Corollary 3 on Bessel stability, and the Hermite–Biehler theorem on p.6 of the author manuscript.
https://people.math.wisc.edu/~fwaleffe/submitted/Jacobi_CCW2.pdf
https://doi.org/10.1090/S0002-9939-05-07898-6

Indeed, with y_m denoting their classical Bessel polynomial,

    B_m(Z):=(2m)!/m! * p_m(Z) = Z^m y_m(2/Z).

The zeros of y_m lie strictly in Re Z<0, so reciprocating shows that B_m is Hurwitz stable too. All zeros of B_m(ix), as a polynomial in x, lie in Im x>0. Multiplying by i^(−m) if necessary makes its leading coefficient real. The Hermite–Biehler theorem says that its real and imaginary parts have only real simple interlacing roots. Up to nonzero real scalar factors these are the even and odd parts of B_m(ix). The latter is

    x C_m(i x).

This is a real polynomial because C_m is even; equivalently, writing
C_m(Z)=R_m(Z^2), it is x R_m(−x^2).

Its zeros are all real and simple, with the zero at x=0 already simple because C_m(0)>0. Thus C_m has only nonzero, purely imaginary simple zeros. This theorem proves reality/simplicity but does not by itself determine the smallest positive zero or give the factorial convergence rate to 2pi; Section 3 supplies the local rate.
